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
#include "apfelcomb/fk_grids.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

using namespace std;


int main(int argc, char* argv[]) {
  
  if (argc!=3)
  {
    cout << "Usage: "<<argv[0]<<" <database id> <theory id>"<<endl;
    exit(1);
  }

  // APPLgrid and theory indices
  const int iDB = atoi(argv[1]);
  const int iTh = atoi(argv[2]);
  setupDir(iTh);

  // Parse parameters
   Splash(); QCD::qcd_param par; QCD::parse_input(iTh, par);
  std::cout << "APPLgrid repository version: "<< applCommit() <<std::endl;

  NNPDF::SetVerbosity(0); appl::setVerbose(false);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", "app_subgrids");

  // Read grid information
  const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,iDB,"fktarget");
  const int target = NNPDF::dbmatch(grid_db, "name", fktarget)[0];
  FKTarget table(grid_db, target); table.ReadSubGrids(subgrid_db);

  DisplayHR();

  // // Initialise QCD
  QCD::initQCD(par, table.GetQ2max());
  QCD::initEvolgrid(table.GetNX(), table.GetXmin());  DisplayHR();

  // Initialise empty mFK table
  NNPDF::FKHeader FKhead; table.SetFKHeader(FKhead); QCD::set_params(par, FKhead);
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

  // Compute FK table
  table.Combine(par, FK);

  DisplayHR();
  cout << "  --  Verify FKTABLE"<<endl;

  APFELPDFSet apfelPDF;

  const vector<double>       xsec = table.Compute(par);
  const NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout << "<data>\t <FK> \t <APPLgrid> \t <Rel. Error.>"<<endl;

  for (int i=0; i<theory.GetNData(); i++)
    {
      const double applpred = xsec[i];
      const double FKpred  = theory.GetObsCV(i);
      const double rel_err = abs((FKpred-applpred)/applpred);
      cout  << setw(5)  << left << i
            << setw(10) << left << FKpred
            << setw(15) << left << applpred
            << setw(15) << left << rel_err
            << endl;
    }

  DisplayHR();


  // // Initialise truth x-grid
  // // Need to get the absolute smallest x-grid value in APPLgrid here, not the value in par, hence the false
  // // This is due to APPLgrids convolution routine requesting the smallest x-grid value from the PDF
  // QCD::initTruthGrid(par, APP::getXmin(sourceGrid.g,false)); 

  // if (par.ppbar == true && par.xiF != 1)
  //   std::cout << "WARNING: ppbar ROTATION NOT TESTED - APPLgrid does not support fac. scale variation with ppbar so I cannot cross-check" <<std::endl;

  // // Compute with applgrid interface
  // const int pto = (par.ptmin == 1) ? -1:(par.pto-1);
  // vector<double> xsec;
  // if (par.ppbar == true && par.xiF == 1)
  //   xsec = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::evolpdf_applgrid_pbar, QCD::alphas, pto, par.xiR, par.xiF );
  // else
  //   xsec = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto, par.xiR, par.xiF);
  // for (double& obs : xsec)
  //   obs *= par.nrmdat;
  
  // cout << "  --  Compute FKTABLE "<<endl;

  // // Re-init evol grid with correct (nonzero) x-grid limits
  // const double targetXmin = APP::parse_xmin(par.common_subgrids);
  // QCD::initEvolgrid(par.nx, targetXmin);   DisplayHR();
  // std::cout << "Evolution initialised with xmin: " << targetXmin << " and nx: " << par.nx <<std::endl;

  // // Setup FastKernel Header
  // NNPDF::FKHeader FKhead;
  // APP::set_params(par, FKhead);

  // // Generate FK table
  // std::stringstream IO; FKhead.Print(IO);
  // NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );
  // APP::computeFK(par,sourceGrid.g,FK);

  // DisplayHR();
  // cout << "  --  Verify FKTABLE"<<endl;

  // APFELPDFSet apfelPDF;
  // NNPDF::ThPredictions theory(&apfelPDF, FK);

  // cout << "<data>\t <FK> \t <APPLgrid> \t <Rel. Error.>"<<endl;

  // for (int i=0; i<par.ndata; i++)
  //   {
  //     const double applpred = xsec[par.maskmap[i]];
  //     const double FKpred  = theory.GetObsCV(par.datamap[i][0]);
  //     const double rel_err = abs((FKpred-applpred)/applpred);
  //     const double targetPrec = APP::computeTargetPrecision(par.datamap[i], cd);
  //     cout  << setw(5)  << left << i
  //           << setw(10) << left << FKpred
  //           << setw(15) << left << applpred
  //           << setw(15) << left << rel_err
  //           << setw(15) << left << targetPrec<<endl;
  //   }

  // DisplayHR();
  // const std::string outFile = getOutputFilename(iTh, par.gridname);
  // ofstream outFK;  outFK.open(outFile.c_str());
  // FK->Print(outFK);
  // outFK.close();

  // DisplayHR();
  // cout << "  --  Verify Export"<<endl;
  // NNPDF::LHAPDFSet nn30("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);
  // NNPDF::FKTable *impFK = new NNPDF::FKTable(outFile);

  // NNPDF::ThPredictions rat = (NNPDF::ThPredictions(&nn30, impFK) - NNPDF::ThPredictions(&nn30, FK))
  //                          /  NNPDF::ThPredictions(&nn30, FK);

  // for (int i=0; i<theory.GetNData(); i++)
  // {
  //   const double rel_err = abs(rat.GetObsCV(i));
  //   if (rel_err > 1E-5)
  //   {
  //     cerr << "Error: FK Table Export Verification failed"<<endl;
  //     rat.Print(cout);
  //     remove(outFile.c_str());
  //     exit(1);
  //   }
  // }

  //  // --  Cleanup ************************************************************

  // delete FK;
  
  cout <<endl<< "--  APPLComb Complete "<<endl;
  DisplayHR();
  exit(0);
}

