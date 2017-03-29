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
  
  if (argc!=4)
  {
    cout << "Usage: "<<argv[0]<<" <source=app/dis/dyp> <database id> <theory id>"<<endl;
    exit(1);
  }

  // APPLgrid and theory indices
  const string source = argv[1];
  const int iDB = atoi(argv[2]);
  const int iTh = atoi(argv[3]);
  setupDir(iTh);

  // Parse parameters
  Splash(); QCD::qcd_param par; QCD::parse_input(iTh, par);
  NNPDF::SetVerbosity(0); appl::setVerbose(false);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", source+"_subgrids");

  // Read grid information
  const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,iDB,"fktarget");
  const int target = NNPDF::dbmatch(grid_db, "name", fktarget)[0];
  FKTarget table(grid_db, target); table.ReadSubGrids(subgrid_db);

  // // Initialise QCD
  if (table.GetSource() == FKTarget::DYP) QCD::setFTDYmode(true);
  QCD::initQCD(par, table.GetPositivity(), table.GetQ2max());
  QCD::initEvolgrid(table.GetNX(), table.GetXmin());  DisplayHR();

  // Initialise empty mFK table
  NNPDF::FKHeader FKhead; table.SetFKHeader(FKhead); QCD::set_params(par, FKhead);
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

  // Compute FK table
  table.GetSubgrid(iDB)->Combine(par, FK);

  cout << "                        Verification                  "<<endl;
  DisplayHR();
  APFELPDFSet apfelPDF;
  const vector<double> xsec = table.Compute(par);
  const NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout  << setw(10) << left << "<idat>"
        << setw(15) << left << "<FK>"
        << setw(15) << left << "<Source>"
        << setw(15) << left << "<Rel.Err>"
        << endl;

  double max_relerr = 0;
  for (auto targets : table.GetSubgrid(iDB)->GetDataMap())
    for (int i : targets)
    {
      const double applpred = xsec[i];
      const double FKpred  = theory.GetObsCV(i);
      const double rel_err = abs((FKpred-applpred)/applpred);
      max_relerr = max(max_relerr, rel_err);
      cout  << setw(10) << left << i
            << setw(15) << left << FKpred
            << setw(15) << left << applpred
            << setw(15) << left << rel_err
            << endl;
    }

  if ( max_relerr > 1E-5 && table.GetSource() != FKTarget::DYP )
  {
    cerr << "Error: Relative error exceeds expectations - something has gone wrong in the combination"<<endl;
    exit(1);
  }

  // // Print to file
  const std::string outFile = resultsPath()+"theory_" + to_string(iTh) + "/subgrids/FK_"+table.GetTargetName()+"_"+to_string(iDB) + ".subgrid.dat";
  ofstream outFK;  outFK.open(outFile.c_str()); 
  FK->Print(outFK); outFK.close();

  DisplayHR();
  NNPDF::FKTable *impFK = new NNPDF::FKTable(outFile);
  NNPDF::ThPredictions rat = (NNPDF::ThPredictions(&apfelPDF, impFK) - NNPDF::ThPredictions(&apfelPDF, FK))
                           /  NNPDF::ThPredictions(&apfelPDF, FK);
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

  //  // --  Cleanup ************************************************************

  delete FK;
  
  cout << "                      APPLComb Complete "<<endl;
  DisplayHR();
  exit(0);
}

