#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"

#include "NNPDF/commondata.h"

using namespace std;

int main(int argc, char* argv[]) {

  if ( argc != 4 && argc != 5 )
  {
    std::cerr << "Usage: " <<argv[0] <<" <setID> <target ThID>" <<std::endl;
    exit(-1);
  }

  Splash();

  // target theoryIDs
  const int tTh = atoi(argv[1]);
  QCD::qcd_param  tPar;
  QCD::parse_input(tTh, tPar);

  exit(0);
}

