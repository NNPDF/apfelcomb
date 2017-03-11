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
  const string setname =  NNPDF::dbquery<string>(grid_db,iDB,"setname");
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
  int lastNdat = 0;
  for (int i: subgridIDs)
  {
    APP::appl_param par;
    APP::parse_input(i, par, true);
    APP::grid* subgrid = new APP::grid(par);
    xmin = std::min(xmin, APP::getXmin(subgrid->g,true));

    subgrid_parameters.push_back(par);
    subgrids.push_back(subgrid);

    // Check NDat
    lastNdat += par.incdat+par.muldat*par.ndata;
  }

  // Compute target interpolation accuracy
  APP::appl_param& bpar = subgrid_parameters[0];
  const std::string commonfile = dataPath() + "commondata/DATA_" + setname + ".dat"; 
  const std::string sysfile    = dataPath() + "commondata/systypes/SYSTYPE_" + setname + "_DEFAULT.dat";  
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(commonfile, sysfile);

  if (lastNdat != cd.GetNData())
  {
    std::cout << xmin <<"  "<<lastNdat<<"  "<<cd.GetNData()<<std::endl;
    std::cerr << "ERROR: no match for " <<setname <<std::endl;
    exit(-1);
  }

  for (int i: subgridIDs)
  {
    const std::vector< std::vector<int> > data_map = APP::parse_map(subgridIDs, i);
    for (auto j : data_map)
      for (auto k : j)
        std::cout << k <<" "<< APP::computeTargetPrecision(j, cd) <<"  "<<cd.GetUncE(k)/cd.GetData(k)<<"  "<<cd.GetCorE(k)/cd.GetData(k)<<std::endl;
  }

  exit(0);
}

