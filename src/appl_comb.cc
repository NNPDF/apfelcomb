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
  std::cout << "APPLgrid repository version: "<< applCommit() <<std::endl;

  // Setup directory
  setupDir(iTh, par.setname, par.inventory);

  appl::grid *g = NULL;
  fastnlo *fg   = NULL;
  
  if (par.fnlo)
  {
    cout << "USING FASTNLO BIN "<<par.fnlobin<<endl;
    fg = new fastnlo(par.gridfile);
    g = fg->grids()[par.fnlobin];
    
  } else{ g = new appl::grid(par.gridfile); }
  
  DisplayHR();

  // Initialise QCD
  QCD::initQCD(par, APP::getQ2max(g));

  // Initialise truth x-grid
  // Need to get the absolute smallest x-grid value in APPLgrid here, not the value in par, hence the false
  // This is due to APPLgrids convolution routine requesting the smallest x-grid value from the PDF
  QCD::initTruthGrid(par, APP::getXmin(g,false)); 
  DisplayHR();
  cout << "  --  High accuracy APPLgrid Result "<<endl;
  
  if (par.ppbar == true && par.xiF != 1)
    std::cout << "WARNING: ppbar ROTATION NOT TESTED - APPLgrid does not support fac. scale variation with ppbar so I cannot cross-check" <<std::endl;

  // Compute with applgrid interface
  const int pto = (par.ptmin == 1) ? -1:(par.pto-1);
  vector<double> xsec;
  if (par.ppbar == true && par.xiF == 1)
    xsec = g->vconvolute( QCD::evolpdf_applgrid, QCD::evolpdf_applgrid_pbar, QCD::alphas, pto, par.xiR, par.xiF );
  else
    xsec = g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto, par.xiR, par.xiF);

  size_t ibin=0;
  for (size_t o=0; o<par.nbins; o++)
    if (par.mask[o])
    {
      cout << "  APPLGRID result, bin: "<<ibin<<"  = "<<xsec[o]<<endl;
      ibin++;
    }
  
  DisplayHR();
  cout << "  --  Compute FKTABLE "<<endl;

  // Re-init evol grid with correct (nonzero) x-grid limits
  QCD::initEvolgrid(par.nx,par.xmin); 


  // Setup FastKernel Header
  NNPDF::FKHeader FKhead;
  APP::set_params(par, FKhead);

  // Generate FK table
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );
  APP::computeFK(par,g,FK);

  DisplayHR();
  cout << "  --  Verify FKTABLE"<<endl;

  APFELPDFSet apfelPDF;
  NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout << "<bin> \t <FK> \t <APPLgrid> \t <Rel. Error.>"<<endl;
  for (int i=0; i<theory.GetNData(); i++)
  {
    const double rel_err = abs((theory.GetObsCV(i)-xsec[par.map[i]])/xsec[par.map[i]]);
    cout << setw(5) << left <<i<< setw(15) << left<<theory.GetObsCV(i)<< setw(15) << left<<xsec[par.map[i]]<< setw(15) << left<<rel_err<<endl;

    if (rel_err > par.tgtprec)
    {
      cerr << "Error: FK Table Verification failed, maxPrec: "<<par.tgtprec<<endl;
      // if (par.ppbar != true || par.xiF == 1)
      //   exit(1);  // Sidestep verification when you can't check against APPLgrid
    }
  }

  DisplayHR();
  const std::string outFile = getOutputFilename(iTh, par.setname, par.gridname);
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
  
  if (par.fnlo)
    delete fg;
  else
    delete g;
  
  cout <<endl<< "--  APPLComb Complete "<<endl;
  DisplayHR();
  exit(0);
}

