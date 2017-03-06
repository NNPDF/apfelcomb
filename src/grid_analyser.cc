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

#include "NNPDF/nnpdfdb.h"


using namespace std;

int main(int argc, char* argv[]) {
  
  if (argc!=2)
  {
    cout << "Usage: grid_analyser <grid id>"<<endl;
    exit(1);
  }

  Splash();
  const int iDB = atoi(argv[1]);
  NNPDF::IndexDB grid_db(databasePath()+"applgrid.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"applgrid.db", "subgrids");

  // Fetch number of entries
  const int entries =grid_db.GetNEntries();
  if (iDB < 0 || iDB > entries)
  {
    cerr << "Error: grid ID ("<<iDB<<") must be between 1 and "<<entries<<endl;
    exit(-1);
  }

  const string name =  NNPDF::dbquery<string>(grid_db,iDB,"name");
  const string desc =  NNPDF::dbquery<string>(grid_db,iDB,"description");
  const int    nx   =  NNPDF::dbquery<int>(grid_db,iDB,"nx");

  std::cout << "Analysing grid: " << name <<std::endl;
  std::cout << desc <<std::endl;
  std::cout << "Target nx: " << nx <<std::endl;

  // Fetch matching subgrids
  const std::vector<int> subgridIDs = NNPDF::dbmatch(subgrid_db, "fktarget", name);
  std::cout << "Found " << subgridIDs.size() << " subgrids" <<std::endl;
  DisplayHR();

  std::vector<APP::grid*> subgrids;
  std::vector<APP::appl_param> subgrid_parameters;
  double xmin = 1.0;
  for (int i: subgridIDs)
  {
    APP::appl_param par;
    APP::parse_input(i, par);
    APP::grid* subgrid = new APP::grid(par);
    xmin = std::min(xmin, APP::getXmin(subgrid->g,true));

    subgrid_parameters.push_back(par);
    subgrids.push_back(subgrid);
  }

  std::cout << xmin <<std::endl;

  exit(0);
}

