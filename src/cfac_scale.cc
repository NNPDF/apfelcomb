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

#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"

#include "NNPDF/commondata.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  if ( argc != 4 )
  {
    std::cerr << "Usage: " <<argv[0] <<" <target ThID> <cFactor_name> <set_name>" <<std::endl;
    exit(-1);
  }

  Splash();

  // base and target theoryIDs
  const int bTh = 3;
  const int tTh = atoi(argv[1]);

  QCD::qcd_param bPar, tPar;
  QCD::parse_input(bTh, bPar);
  QCD::parse_input(tTh, tPar);

  const std::string cFacName = argv[2];
  const std::string cFacPath = dataPath() + "/NNLOCFAC/CF_QCD_"+cFacName+".dat";
  const std::string cFacOut  = dataPath() + "/theory_" + std::to_string(tTh) + "/cfactors/CF_QCD_"+cFacName+".dat";

  const std::string setName = argv[3];
  const std::string cDataPath = dataPath() + "/commondata/DATA_"+setName+".dat";
  const std::string sysTypePath = dataPath() + "/commondata/systypes/SYSTYPE_"+setName+"_0.dat";

  // Read CommonData 
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(cDataPath, sysTypePath);

  double* alphas_0 = new double[cd.GetNData()];
  double* alphas_1 = new double[cd.GetNData()];

  // Read Cfactors
  double* cFac = new double[cd.GetNData()];
  std::ifstream instream(cFacPath.c_str());
  std::ofstream outstream(cFacOut.c_str());

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

  std::cout <<std::setw(8)<< "Scale"<<"\t"<< std::setw(8)<< "OldCFAC" << "\t"<<std::setw(8)<<"NewCFAC"<<std::endl;
  for (int i=0; i< cd.GetNData(); i++)
  {
    const double fac1 = 1.0 + (cFac[i]-1.0)*pow(alphas_1[i]/alphas_0[i],2.0);
    outstream << fac1 <<endl;
    std::cout <<std::setw(8)<< std::sqrt(cd.GetKinematics(i,1)) <<"\t"<< std::setw(8)<< cFac[i] << "\t"<<std::setw(8)<<fac1<<std::endl;
  }

  std::cout << "CFactors for " << setName << " converted" <<std::endl;

  instream.close();
  outstream.close();

  exit(0);
}

