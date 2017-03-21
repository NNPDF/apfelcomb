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

  // Fetch number of entries
  const int entries =grid_db.GetNEntries();
  if (iDB < 0 || iDB > entries)
  {
    cerr << "Error: grid ID ("<<iDB<<") must be between 1 and "<<entries<<endl;
    exit(-1);
  }

  FKTarget table(grid_db, iDB);
  table.ReadSubGrids(subgrid_db);
  table.Splash(std::cout);


  QCD::qcd_param par;
  QCD::parse_input(52, par);
  QCD::initQCD(par, table.GetQ2max());
  QCD::initTruthGrid(par, table.GetComputeXmin()); 
  const vector<double> xsec = table.Compute(par);
  for (size_t i=0; i<xsec.size(); i++)
    std::cout << i <<"\t"<<xsec[i]<<std::endl;

  exit(0);
}

