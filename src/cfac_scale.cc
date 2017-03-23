// cfac_scale.cc
// Code to perform scaling of c-factors according to alpha_S variations
// Scaling:  C(alphas_1) = 1.0 + (C(alphas_0)-1.0)*pow(alphas_1/alphas_0,pto);
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
    std::cerr << "Usage: " <<argv[0] <<" <target ThID> <cFactor_name> <set_name> [MZ/MW/MT]" <<std::endl;
    exit(-1);
  }

  Splash();

  // target theoryIDs
  const int tTh = atoi(argv[1]);
  QCD::qcd_param  tPar;
  QCD::parse_input(tTh, tPar);

  const int bTh = 3;
  QCD::qcd_param  bPar;
  QCD::parse_input(bTh, bPar);

  const std::string scalearg = argc == 5 ? argv[4]:"";
  const bool bMW =  scalearg=="MW" ? true:false; const double MW = atof(tPar.thMap["MW"].c_str()); 
  const bool bMZ =  scalearg=="MZ" ? true:false; const double MZ = atof(tPar.thMap["MZ"].c_str()); 
  const bool bMT =  scalearg=="MT" ? true:false; const double MT = atof(tPar.thMap["mt"].c_str()); 
  const bool bfixed = bMW || bMZ || bMT; const double fixed = bMW ? MW:(bMZ ? MZ:(bMT? MT:0));
  if (bfixed) std::cout << "Using "<<scalearg<<" ("<<fixed<<") as scale" <<std::endl;

  if (tPar.evol_pto != 1 && tPar.evol_pto !=2) 
  {    
    std::cerr << "TheoryID must be either NLO or NNLO for scaling C-factors!" <<std::endl;
    exit(-1);
  }

  const std::string ptoString = (tPar.evol_pto == 1) ? "NLO":"NNLO";

  const std::string cFacName = argv[2];
  const std::string cFacPath = dataPath() + "/"+ptoString+"CFAC/CF_QCD_"+cFacName+".dat";
  const std::string cFacOut  = dataPath() + "/theory_" + to_string(tTh) + "/cfactor/CF_QCD_"+cFacName+".dat";

  const std::string setName = argv[3];
  const std::string cDataPath = dataPath() + "/commondata/DATA_"+setName+".dat";
  const std::string sysTypePath = dataPath() + "/commondata/systypes/SYSTYPE_"+setName+"_DEFAULT.dat";

  // Read CommonData 
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(cDataPath, sysTypePath);

  double* alphas_0 = new double[cd.GetNData()];
  double* alphas_1 = new double[cd.GetNData()];

  // Read Cfactors
  double* cFac = new double[cd.GetNData()];
  std::ifstream instream(cFacPath.c_str());
  std::ofstream outstream(cFacOut.c_str());


  if (instream.good() == false)
  {
    std::cerr << "Bad stream: " <<cFacPath.c_str()<<std::endl;
    exit(-1);
  }

  for (int i=0; i< 9; i++)
  {
    std::string dum;
    getline(instream,dum);
    if ( i == 5 )
      outstream << dum<<" -> Converted to alpha_s: " << tPar.thMap["alphas"]<< endl;
    else
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
    alphas_0[i] = QCD::alphas( bfixed ? fixed:std::sqrt(cd.GetKinematics(i,1)));

  // alphas_1
  QCD::initQCD(tPar, DIS::getQ2max(cd));
  for (int i=0; i< cd.GetNData(); i++)
    alphas_1[i] = QCD::alphas( bfixed ? fixed:std::sqrt(cd.GetKinematics(i,1)));

  std::cout <<std::setw(8)<< "Scale"<<"\t"<< std::setw(8)<< "OldCFAC" << "\t"<<std::setw(8)<<"NewCFAC"<<std::endl;
  for (int i=0; i< cd.GetNData(); i++)
  {
    const double fac1 = 1.0 + (cFac[i]-1.0)*pow(alphas_1[i]/alphas_0[i], tPar.pto);
    outstream << fac1 <<endl;

    // Print to screen
    std::cout << std::setw(8)
              << std::sqrt(cd.GetKinematics(i,1)) 
              <<"\t"<< std::setw(8)
              << cFac[i] 
              << "\t"<<std::setw(8)
              <<fac1
              <<std::endl;
  }

  std::cout << "CFactors for " << setName << " converted" <<std::endl;

  instream.close();
  outstream.close();

  exit(0);
}

