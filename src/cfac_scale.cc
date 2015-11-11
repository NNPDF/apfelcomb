// cfac_scale.cc
// Code to perform scaling of c-factors according to alpha_S variations
// Scaling:  C(alphas_1) = 1.0 + (C(alphas_0)-1.0)*pow(alphas_1/alphas_0,2.0);
#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>

#include "fk_dis.h"
#include "fk_utils.h"
#include "fk_pdf.h"
#include "fk_qcd.h"

#include "NNPDF/commondata.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  Splash();

  // base and target theoryIDs
  const int bTh = 3;
  const int tTh = 8;

  QCD::qcd_param bPar, tPar;
  QCD::parse_input(bTh, bPar);
  QCD::parse_input(tTh, tPar);

  const std::string cFacName = argv[1];
  const std::string cFacPath = "../nnpdfcpp/data/NNLOCFAC_0119/CF_QCD_"+cFacName+".dat";
  const std::string cFacOut = "../nnpdfcpp/data/NNLOCFAC_0118/CF_QCD_"+cFacName+".dat";

  const std::string setName = argv[2];
  const std::string cDataPath = "../nnpdfcpp/data/commondata/DATA_"+setName+".dat";
  const std::string sysTypePath = "../nnpdfcpp/data/commondata/systypes/SYSTYPE_"+setName+"_0.dat";

  // Read CommonData 
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(cDataPath, sysTypePath);

  double* alphas_0 = new double[cd.GetNData()];
  double* alphas_1 = new double[cd.GetNData()];

  // Read Cfactors
  double* cFac = new double[cd.GetNData()];
  std::ifstream instream(cFacPath);
  std::ofstream outstream(cFacOut);

  for (int i=0; i< 9; i++)
  {
    std::string dum;
    getline(instream,dum);
    outstream << dum<<endl;
  }

  for (int i=0; i< cd.GetNData(); i++)
  {
    std::string dum; stringstream cvrt;
    getline(instream,dum);
    cvrt << dum;
    cvrt >> cFac[i];
  }
  // alphas_0
  QCD::initQCD(bPar, DIS::getQ2max(cd));
  for (int i=0; i< cd.GetNData(); i++)
    alphas_0[i] = QCD::alphas( std::sqrt(cd.GetKinematics(i,1)));

  // alphas_1
  QCD::initQCD(tPar, DIS::getQ2max(cd));
  for (int i=0; i< cd.GetNData(); i++)
    alphas_1[i] = QCD::alphas( std::sqrt(cd.GetKinematics(i,1)));

  for (int i=0; i< cd.GetNData(); i++)
  {
    const double fac1 = 1.0 + (cFac[i]-1.0)*pow(alphas_1[i]/alphas_0[i],2.0);
    const double perDiff = 100.0*fabs(cFac[i]-fac1)/cFac[i];

    std::cout << cFac[i]<<"  "<<fac1 <<"  "<<std::sqrt(cd.GetKinematics(i,1))<<"  "<<100.0*fabs(cFac[i]-fac1)/cFac[i]<<std::endl; 
    outstream << fac1 <<endl;
  }

  instream.close();
  outstream.close();

  exit(0);
}

