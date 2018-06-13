// ftdy_comb.cc
// Code to generate FK tables from FTDY hard xsec tables
//
// nph  03/15

#include <iostream>
#include <cstdlib>
#include <string>
#include <cstdio>

#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_ftdy.h"

#include "APFEL/APFEL.h"

using namespace std;


int main(int argc, char* argv[]) {

    //  Currently out of action
//  for (int i=13; i<FTDY::nsets; i++)
//  {
//      cout << "Computing: "<<FTDY::setnames[i]<< "  ("<<i<<"/"<<FTDY::nsets<<")"<<endl;
//
//      const std::string commonfile = dataPath() + "commondata/DATA_" + FTDY::setnames[i] + ".dat";
//      const std::string hcxfile = "./data/FTDY/"+FTDY::setnames[i]+".hcx";
//
//      APFEL::ComputeHardCrossSectionsDY(commonfile,hcxfile);
//  }

  exit(0);
}

