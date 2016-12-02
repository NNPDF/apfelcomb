// applconv.cc
// Code to generate FastKernel-Like grids
// from applgrid files and evolution grids.
//
// Input: appl_comb parameter file
// Output: Combined FastKernel grid
//
// n.p.hartland@ed.ac.uk  03/12

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

#include "NNPDF/fastkernel.h"
#include "NNPDF/thpredictions.h"
#include "NNPDF/lhapdfset.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  if ( argc<3 || argc > 5)
  {
    cout << "Usage: dis_comb <database id> <theory id> [database file] [timelike evol (0/1)]"<<endl;
    exit(1);
  }
  
  Splash();
  
  // APPLgrid and theory indices
  const int iDB = atoi(argv[1]);
  const int iTh = atoi(argv[2]);
  const bool timelike = ( argc == 5 ) ? ((atoi(argv[4]) == 1) ? true:false):false;

  // Set mode flags
  QCD::setDISmode(true);
  if (timelike) QCD::setSIAmode(true);

  // Parse parameters
  DIS::dis_param par;
  QCD::parse_input(iTh, par);

  // Parse DIS information
  if (argc>3)
    DIS::parse_input(iDB, par, std::string(argv[3]));
  else
    DIS::parse_input(iDB, par, std::string("dis.db"));

  // Setup directory
  setupDir(iTh, par.setname, par.inventory);

  // Read CommonData
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(par.commonfile, par.sysfile);
  DisplayHR();

  // Initialise QCD
  QCD::initQCD(par, DIS::getQ2max(cd));

  par.nx = 50;
  par.maxprec = 1E-7;
  par.ndata = cd.GetNData();
  par.xmin = DIS::getXmin(cd);

  DisplayHR();
  cout << "  --  APFEL Result"<<endl;
  QCD::initEvolgrid(par.nx,DIS::getXmin(cd));

  // Vector for truth values
  vector<double> xsec;

  // Timers
  timeval t1, t2;
  double elapsedTime;
  
  gettimeofday(&t1, NULL);
  for (int i=0; i<cd.GetNData(); i++)
  {
    const double x = cd.GetKinematics(i,0);
    const double Q = sqrt(cd.GetKinematics(i,1));
    const double y = cd.GetKinematics(i,2);

    xsec.push_back(QCD::disobs(par.process, x, Q, y));
  }
  gettimeofday(&t2, NULL);

  for (int o=0; o<cd.GetNData(); o++)
      cout << "APFEL result, bin: "<<o<<"  = "<<xsec[o]<<endl;

    // compute and print the elapsed time in millisec
  elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
  elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
  cout << Colour::FG_GREEN;
  cout << "Total time: "<<elapsedTime << " ms. ";
  cout << "Time per observable: "<<elapsedTime/((double) xsec.size()*(10)) << " ms.\n";
  cout << Colour::FG_DEFAULT;

  DisplayHR();
  cout << " --  Compute FKTABLE "<<endl;

  // Setup FastKernel Header
  NNPDF::FKHeader FKhead;
  DIS::set_params(par, FKhead);

  // Generate FK table
  std::stringstream IO; FKhead.Print(IO);
  NNPDF::FKGenerator* FK = new NNPDF::FKGenerator( IO );
  DIS::computeFK(par,cd,FK);

  DisplayHR();
  cout << "  --  Verify FKTABLE "<<endl;

  APFELPDFSet apfelPDF;
  NNPDF::ThPredictions theory(&apfelPDF, FK);

  cout << "<bin> \t <FK> \t <APFEL> \t <Rel. Error.> \t <Data>"<<endl;
  for (int i=0; i<theory.GetNData(); i++)
  {
    const double rel_err = abs((theory.GetObsCV(i)-xsec[i])/xsec[i]);
    cout << setw(5) << left <<i<< setw(15) << left<<theory.GetObsCV(i)<< setw(15) << left<< xsec[i]<< setw(15) << left<<rel_err<<setw(15) << left<< cd.GetData()[i]<<endl;

    if (rel_err > 1E-2)
    {
      cerr << "Error: FK Table Verification failed"<<endl;
      exit(1);
    }
  }

  DisplayHR();
  cout << "  --  Exporting FKTABLE"<<endl;

  const std::string outFile = getOutputFilename(iTh, par.setname, par.gridname);
  ofstream outFK;  outFK.open(outFile.c_str());
  
  FK->Print(outFK);
  outFK.close();

  DisplayHR();
  cout << "  --  Verifying Export"<<endl;

  NNPDF::LHAPDFSet nn30("NNPDF30_nlo_as_0118", NNPDF::PDFSet::ER_MC);
  NNPDF::FKTable *impFK = new NNPDF::FKTable(outFile);

  NNPDF::ThPredictions rat = (NNPDF::ThPredictions(&nn30, impFK) - NNPDF::ThPredictions(&nn30, FK))
                           /  NNPDF::ThPredictions(&nn30, FK);

  for (int i=0; i<theory.GetNData(); i++)
  {
    const double rel_err = abs(rat.GetObsCV(i));
    if (rel_err > 1E-5)
    {
      cerr << "Error: FK Table Export Verification failed"<<endl;
      rat.Print(cout);
      remove(outFile.c_str());
      exit(1);
    }
  }

   // --  Cleanup ************************************************************

  delete FK;

  DisplayHR();
  cout << "  --  DISComb Complete "<<endl;
  
  exit(0);
}

