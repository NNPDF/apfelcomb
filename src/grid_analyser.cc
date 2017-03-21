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
#include "apfelcomb/fk_grids.h"

#include "NNPDF/nnpdfdb.h"
#include "NNPDF/commondata.h"

using namespace std;

int main(int argc, char* argv[]) {

  if (argc!=2)
  {
    cout << "Usage: grid_analyser <grid id>"<<endl;
    exit(1);
  }

  Splash(); NNPDF::SetVerbosity(0); appl::setVerbose(false);
  const int iDB = atoi(argv[1]);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", "app_subgrids");

  FKTarget table(grid_db, iDB);
  table.ReadSubGrids(subgrid_db);
  table.Splash(std::cout);

  QCD::qcd_param par;
  QCD::parse_input(52, par);
  QCD::initQCD(par, table.GetQ2max());
  QCD::initTruthGrid(par, table.GetComputeXmin()); 
  const vector<double> xsec = table.Compute(par);

  double chi2 = 0;
  for (size_t i=0; i<xsec.size(); i++)
  {
    chi2 += pow(xsec[i] - table.GetCommonData().GetData(i), 2)/pow(table.GetCommonData().GetUncE(i),2);
    std::cout << i <<"\t"<<xsec[i]<<"\t"<<table.GetCommonData().GetData(i)<<std::endl;
  }
  std::cout << "Uncorrelated chi2: " << chi2/(double)xsec.size()<<std::endl;
  exit(0);
}

