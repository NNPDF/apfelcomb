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

using namespace std;

void exportGrid(QCD::qcd_param const& par, NNPDF::CommonData const& cd, std::string const& gridname, std::string const& outFile)
{
  // Initalise FK table
  std::stringstream desc;
  desc  <<  "-----------------------------------------------------------"<<std::endl
        <<  " FK_"<<gridname<<".dat"<<std::endl
        <<  "-----------------------------------------------------------";

  // Setup FastKernel Header
  NNPDF::FKHeader FKhead;
  QCD::set_params(par, FKhead);

  FKhead.AddTag(NNPDF::FKHeader::BLOB, "GridDesc", desc.str());
  FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "SETNAME", cd.GetSetName());
  FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "NDATA", cd.GetNData());
  FKhead.AddTag(NNPDF::FKHeader::GRIDINFO, "HADRONIC", true );

  // Generate FK table
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );

  FTDY::processFK(par, cd, gridname, FK);

  APFELPDFSet apfelPDF;
  NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout << "<bin> \t <FK>"<<endl;
  for (int i=0; i<theory.GetNData(); i++)
    cout << setw(5) << left <<i<< setw(15) << left<<theory.GetObsCV(i)<<endl;

  // Export grid
  ofstream outFK;   outFK.open(outFile.c_str());
  FK->Print(outFK); outFK.close();

  delete FK;
}


int main(int argc, char* argv[]) {
  
  if (argc!=3)
  {
    cout << "Usage: ftdy_comb <dataset id> <theory id>"<<endl;
    cout << "        Available dataset IDs: "<<endl;
    for (int i=0; i<FTDY::nsets; i++)
      cout << "        "<<i+1<<" : "<<FTDY::setnames[i]<<endl;
    exit(1);
  }

  Splash();

  // Setup FTDY mode
  QCD::setFTDYmode(true);

  const int iDt = atoi(argv[1]);
  const int iTh = atoi(argv[2]);
  
  QCD::qcd_param par;
  QCD::parse_input(iTh, par);

  // Setup directory
  const std::string setname = FTDY::setnames[iDt-1];
  setupDir(iTh, setname);

  const std::string commonfile = dataPath() + "commondata/DATA_" + setname + ".dat"; //!< Path for the commondata file
  const std::string sysfile    = dataPath() + "commondata/systypes/SYSTYPE_" + setname + "_0.dat"; //!< Path for the SYSTYPE file
   
  // Read CommonData
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(commonfile, sysfile);

  const double nx = FTDY::nx[iDt-1];
  cout << "Using " << nx << " x-grid points, xmin = "<<FTDY::getXmin(cd)<<endl;

  // Initialise QCD
  QCD::initQCD(par, FTDY::getQ2max(cd));
  QCD::initEvolgrid(nx, FTDY::getXmin(cd)); // Need to be a little below xmin according to vb

  // Compute FK grids
  FTDY::computeGrid(par,cd);

  // Export the grid
  if (setname == "DYE886R")
  {
    exportGrid(par, cd, "DYE886R_P", getOutputFilename(iTh, setname, "DYE886R_P"));
    exportGrid(par, cd, "DYE886R_D", getOutputFilename(iTh, setname, "DYE886R_D"));
  }
  else
  {
    exportGrid(par, cd, setname, getOutputFilename(iTh, setname, setname));
  }

  cout <<endl<< "--  FTDYComb Complete **********************************"<<endl;
  
  exit(0);
}

