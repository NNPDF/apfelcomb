// ftdy_comb.cc
// Code to generate FK tables from FTDY hard xsec tables
//
// nph  03/15

#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdio>
#include <vector>
#include <string>

#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_ftdy.h"

#include "APFEL/APFEL.h"

using namespace std;


int main(int argc, char* argv[]) {

  std::vector<std::string> setnames = {"POSDYC", "POSDYCBSB", "POSDYCSB", "POSDYUBD", "POSDYUD", "POSDYCBD", "POSDYCD", "POSDYD", "POSDYUBDB", "POSDYUDB", "POSDYCBDB", "POSDYCDB", "POSDYS", "POSDYUBS", "POSDYUS", "POSDYCBS", "POSDYCS", "POSDYU", "POSDYUBSB", "POSDYUSB"};
  
  int nsets = setnames.size();

    //  Currently out of action
  for (int i=0; i<nsets; i++)
  {
      cout << "Computing: "<<setnames[i]<< "  ("<<i<<"/"<<nsets<<")"<<endl;

      const std::string commonfile = dataPath() + "commondata/DATA_" + setnames[i] + ".dat";
      const std::string hcxfile = "./data/FTDY/"+setnames[i]+".hcx";

      APFEL::ComputeHardCrossSectionsDY(commonfile,hcxfile);
  }

  exit(0);
}

