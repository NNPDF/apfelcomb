// ftdy_comb.cc
// Code to generate FK tables from FTDY hard xsec tables
//
// nph  03/15

#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_ftdy.h"
#include "fk_pdf.h"

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>

#include "APFEL/APFEL.h"

using namespace std;


int main(int argc, char* argv[]) {
  
  for (int i=13; i<FTDY::nsets; i++)
  {
      cout << "Computing: "<<FTDY::setnames[i]<< "  ("<<i<<"/"<<FTDY::nsets<<")"<<endl;

      const std::string commonfile = dataPath() + "commondata/DATA_" + FTDY::setnames[i] + ".dat";
      const std::string hcxfile = "./data/FTDY/"+FTDY::setnames[i]+".hcx";

      APFEL::ComputeHardCrossSectionsDY(commonfile,hcxfile);
  }

  exit(0);
}

