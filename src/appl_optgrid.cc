// optgrid.cc
// Code to optimise x-grid nX
//
// nph  09/14

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <numeric>
#include <sys/time.h>

#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_qcd.h"

using namespace std;

double trgPrec(int const& i, NNPDF::CommonData const& cd)
{
  double targetPrec = std::numeric_limits<double>::infinity();
  if (abs(cd.GetUncE(i)/cd.GetData(i)) > 1E-5)
    targetPrec = std::min(targetPrec, abs(cd.GetUncE(i)/cd.GetData(i))/5.0);
  if (targetPrec == std::numeric_limits<double>::infinity())
    targetPrec = std::min(targetPrec, abs(cd.GetCorE(i)/cd.GetData(i))/10.0);
  if (targetPrec == std::numeric_limits<double>::infinity())
  {
    targetPrec = 0.001;
    std::cout << "WARNING: NO ERROR AVAILABLE, SETTING TO PERMILLE ACCURACY" <<std::endl;
  }
  return targetPrec;
}


int main(int argc, char* argv[]) {

  if (argc!=3)
  {
    cout << "Usage: "<<argv[0]<<" <grid id> <theory id>"<<endl;
    exit(1);
  }

  // APPLgrid and theory indices
  const int iDB = atoi(argv[1]);
  const int iTh = atoi(argv[2]);
  setupDir(iTh);

  // Parse parameters
  Splash(); QCD::qcd_param par; QCD::parse_input(iTh, par);

  NNPDF::SetVerbosity(0);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  const string source = "app";
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", source+"_subgrids");

  // Read grid information
  FKTarget table(grid_db, iDB, par.global_nx);
  table.ReadSubGrids(subgrid_db);

  // // Initialise QCD
  QCD::initQCD(par, table.GetPositivity(), table.GetQ2max());
  // QCD::initPDF("NNPDF30_nlo_as_0118.LHgrid",3);

  DisplayHR();
  cout << "                     Testing Predictions"<<endl;

  // Compute 'ideal' predictions
  double maxdiff = 1; int nx = 10;
  vector<double> xsec_prev(table.GetCommonData().GetNData(), std::numeric_limits<double>::infinity());
  while (maxdiff >= 1)
  {
    QCD::initEvolgrid(nx+=5, table.GetXmin());
    const vector<double> xsec_test = table.Compute(par);

    double maxdiff = 0;
    for (int i=0; i<xsec_prev.size(); i++)
      maxdiff = max(maxdiff, abs((xsec_prev[i] - xsec_test[i])/xsec_test[i])/trgPrec(i, table.GetCommonData()));
    xsec_prev = xsec_test;
    std::cout << "  - MAXDIFF: " <<maxdiff<<std::endl;
  }


  // cout << "                        Verification                  "<<endl;
  // DisplayHR();
  // APFELPDFSet apfelPDF;
  // const vector<double> xsec = table.Compute(par);
  // const NNPDF::ThPredictions theory(&apfelPDF, FK);

  // cout  << setw(10) << left << "<idat>"
  //       << setw(15) << left << "<FK>"
  //       << setw(15) << left << "<Source>"
  //       << setw(15) << left << "<Rel.Err>"
  //       << endl;

  // double max_relerr = 0;
  // for (auto targets : table.GetSubgrid(iDB)->GetDataMap())
  //   for (int i : targets)
  //   {
  //     const double applpred = xsec[i];
  //     const double FKpred  = theory.GetObsCV(i);
  //     const double rel_err = abs((FKpred-applpred)/applpred);
  //     max_relerr = max(max_relerr, rel_err);
  //     cout  << setw(10) << left << i
  //           << setw(15) << left << FKpred
  //           << setw(15) << left << applpred
  //           << setw(15) << left << rel_err
  //           << endl;
  //   }

  cout << "                      OptGrid Complete "<<endl;
  DisplayHR();

  exit(0);
}

