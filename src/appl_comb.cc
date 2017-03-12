// applconv.cc
// Code to generate FastKernel-Like grids
// from applgrid files and evolution grids.
//
// Input: appl_comb parameter file
// Output: Combined FastKernel grid
//
// n.p.hartland@ed.ac.uk  03/12

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

using namespace std;

bool verifyGrid(APP::appl_param const& param, const appl::grid * g )
{
  if (param.ndata>g->Nobs())
  {
    cerr <<"Error: number of datapoints in settings: "<<param.ndata<<" > appl grid Nobs: "<<g->Nobs()<<endl;
    return false;
  }
  
  if (param.mask.size()!=g->Nobs())
  {
    cerr << "Error: mask size: "<<param.mask.size()<<" does not match number of bins in APPLgrid: "<<g->Nobs()<<endl;
    return false;
  }

  return true;
}

int main(int argc, char* argv[]) {
  
  if (argc!=3)
  {
    cout << "Usage: appl_comb <database id> <theory id>"<<endl;
    exit(1);
  }

  Splash();
  
  // APPLgrid and theory indices
  const int iDB = atoi(argv[1]);
  const int iTh = atoi(argv[2]);

  // Parse parameters
  APP::appl_param par;
  QCD::parse_input(iTh, par);
  APP::parse_input(iDB, par);
  par.datamap = APP::parse_map(par.common_subgrids, iDB);
  std::cout << "APPLgrid repository version: "<< applCommit() <<std::endl;

  const std::string commonfile = dataPath() + "commondata/DATA_" + par.setname + ".dat"; 
  const std::string sysfile    = dataPath() + "commondata/systypes/SYSTYPE_" + par.setname + "_DEFAULT.dat";  
  const NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(commonfile, sysfile);

  // Setup directory
  setupDir(iTh, par.setname);
  APP::grid sourceGrid(par);

  if (!verifyGrid(par, sourceGrid.g)) exit(-1);
  DisplayHR();

  // Initialise QCD
  QCD::initQCD(par, APP::getQ2max(sourceGrid.g));

  // Initialise truth x-grid
  // Need to get the absolute smallest x-grid value in APPLgrid here, not the value in par, hence the false
  // This is due to APPLgrids convolution routine requesting the smallest x-grid value from the PDF
  QCD::initTruthGrid(par, APP::getXmin(sourceGrid.g,false)); 

  if (par.ppbar == true && par.xiF != 1)
    std::cout << "WARNING: ppbar ROTATION NOT TESTED - APPLgrid does not support fac. scale variation with ppbar so I cannot cross-check" <<std::endl;

  // Compute with applgrid interface
  const int pto = (par.ptmin == 1) ? -1:(par.pto-1);
  vector<double> xsec;
  if (par.ppbar == true && par.xiF == 1)
    xsec = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::evolpdf_applgrid_pbar, QCD::alphas, pto, par.xiR, par.xiF );
  else
    xsec = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto, par.xiR, par.xiF);
  for (double& obs : xsec)
    obs *= par.nrmdat;
  
  cout << "  --  Compute FKTABLE "<<endl;

  // Re-init evol grid with correct (nonzero) x-grid limits
  const double targetXmin = APP::parse_xmin(par.common_subgrids);
  QCD::initEvolgrid(par.nx, targetXmin);   DisplayHR();
  std::cout << "Evolution initialised with xmin: " << targetXmin << " and nx: " << par.nx <<std::endl;

  // Setup FastKernel Header
  NNPDF::FKHeader FKhead;
  APP::set_params(par, FKhead);

  // Generate FK table
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );
  APP::computeFK(par,sourceGrid.g,FK);

  DisplayHR();
  cout << "  --  Verify FKTABLE"<<endl;

  APFELPDFSet apfelPDF;
  NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout << "<data>\t <FK> \t <APPLgrid> \t <Rel. Error.>"<<endl;

  for (int i=0; i<par.ndata; i++)
    {
      const int iappl = par.maskmap[i];
      const std::string point = to_string(i) + " (" + to_string(iappl) + ")";
      const double rel_err = abs((theory.GetObsCV(i)-xsec[iappl])/xsec[iappl]);
      const double targetPrec = APP::computeTargetPrecision(par.datamap[i], cd);
      cout << setw(10) << left <<point<< setw(10) << left<<theory.GetObsCV(i)<< setw(15) << left<<xsec[iappl]<< setw(15) << left<<rel_err<< setw(15) << left<<targetPrec<<endl;
    }

  DisplayHR();
  const std::string outFile = getOutputFilename(iTh, par.gridname);
  ofstream outFK;  outFK.open(outFile.c_str());
  FK->Print(outFK);
  outFK.close();

  DisplayHR();
  cout << "  --  Verify Export"<<endl;
  NNPDF::LHAPDFSet nn30("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);
  NNPDF::FKTable *impFK = new NNPDF::FKTable(outFile);

  NNPDF::ThPredictions rat = (NNPDF::ThPredictions(&nn30, impFK) - NNPDF::ThPredictions(&nn30, FK))
                           /  NNPDF::ThPredictions(&nn30, FK);

  for (int i=0; i<theory.GetNData(); i++)
  {
    const double rel_err = abs(rat.GetObsCV(i));
    if (rel_err > 1E-5)
    {
      cerr << "Error: FK Table Export Verification failed"<<endl;
      rat.Print(cout);
      remove(outFile.c_str());
      exit(1);
    }
  }

   // --  Cleanup ************************************************************

  delete FK;
  
  cout <<endl<< "--  APPLComb Complete "<<endl;
  DisplayHR();
  exit(0);
}

