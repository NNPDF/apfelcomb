#include "LHAPDF/LHAPDF.h"

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <sys/time.h>
#include <cstdio>
#include <limits>

#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_pdf.h"
#include "apfelcomb/fk_qcd.h"

#include "NNPDF/commondata.h"

using namespace std;

int main(int argc, char* argv[]) {

  if ( argc != 6 )
  {
    std::cerr << "Usage: " <<argv[0] <<"<source=app/dis/dyp> <database id> <theory id> <set_name> <proton_baseline>" <<std::endl;
    exit(-1);
  }

  // Paths to CommonData
  // APPLgrid and theory indices
  const string source = argv[1];
  const int iDB = atoi(argv[2]);
  const int iTh = atoi(argv[3]);
  const string setName = argv[4];
  const string pdfset=argv[5];
  setupDir(iTh);

  const std::string cDataPath = dataPath() + "/commondata/DATA_"+setName+".dat";
  const std::string sysTypePath = dataPath() + "/commondata/systypes/SYSTYPE_"+setName+"_DEFAULT.dat";

  cout<<"dataPath() = "<<dataPath()<<endl;
  // Read CommonData
  NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(cDataPath, sysTypePath);


  Splash();
  // target theoryIDs
  QCD::qcd_param par;
  QCD::parse_input(iTh, par);

  //Read Grids and subgrids
  NNPDF::SetVerbosity(0);
  NNPDF::IndexDB grid_db(databasePath()+"apfelcomb.db", "grids");
  NNPDF::IndexDB subgrid_db(databasePath()+"apfelcomb.db", source+"_subgrids");

  // Read grid information
  const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,iDB,"fktarget");
  const std::vector<int> grid_matches = NNPDF::dbmatch(grid_db, "name", fktarget);

  if ( fktarget == "" ) // A bit clumsy, doing better would require modifying dbquery to throw an error when no match is found
  throw std::runtime_error("Cannot find subgrid "+to_string(iDB)+" in apfelcomb database");
  if ( grid_matches.size() == 0 )
  throw std::runtime_error("Cannot find FK target "+fktarget+" in apfelcomb database");

  const int target = grid_matches[0];
  FKTarget table(grid_db, target, par.global_nx); table.ReadSubGrids(subgrid_db);

  QCD::initQCD(par, table.GetPositivity(), table.GetQ2max());
  QCD::initEvolgrid(table.GetNX(), table.GetXmin());
  double F2_D;

  double x,Q,Q2,y;

  string command1 = "mkdir "+dataPath() + "/theory_" + to_string(iTh) + "/F2D_unnormalization/"+setName;
  system(command1.c_str());

  int Nrep=100;//100
  for(int nrep=0; nrep<1; nrep++)
  {
    std::string F2DPath = dataPath() + "/theory_" + to_string(iTh) + "/F2D_unnormalization/"+setName+"/F2D_"+setName+"_"+to_string(nrep)+".dat";
    std::ofstream outstream(F2DPath.c_str());
    outstream<<"*******************************************************************************************"<<endl;
    outstream<<"SetName: "<<setName<<endl;
    outstream<<"Author: Rabah Abdul Khalek rabah.khalek@gmail.com"<<endl;
    outstream<<"Date: 2018"<<endl;
    outstream<<"CodesUsed: None"<<endl;
    outstream<<"TheoryInput: theoryID "<<iTh<<endl;
    outstream<<"PDFset: "<<pdfset<<endl;
    outstream<<"Warnings: None"<<endl;
    outstream<<"Replica (LHAPDF): "<<nrep<<endl;
    outstream<<"Description: Unnormalization corrections for F2_A/F2_D nuclear ratios"<<endl;
    outstream<<"*******************************************************************************************"<<endl;

    bool FOR_APFELGRID=0;

    if(FOR_APFELGRID)
    outstream<<endl;



    for (int ndat=0; ndat<cd.GetNData(); ndat++)
    {

      if(setName == "nEMCC" || setName == "nEMCCA" || setName == "nEMCFE")
      {
        if(FOR_APFELGRID)
        outstream<<1.<<endl;//" "<<0<<endl;
        else
        outstream<<1.<<" "<<0<<endl;
      }
      else
      {
        x=cd.GetKinematics(ndat,0);
        Q2=cd.GetKinematics(ndat,1);
        Q=sqrt(Q2);
        y=cd.GetKinematics(ndat,2);

        QCD::initPDF(pdfset+".LHgrid",nrep);
        string obs;
        //if(setName== nE139AGD nE139ALD nE139AUD nE139BED nE139CAD nE139CD  nE139FED nE139HED
        //if(setName.find("nE139"))
        if(setName == "nE139AGD" || setName ==  "nE139ALD" || setName == "nE139AUD" ||
           setName == "nE139BED" || setName ==  "nE139CAD" || setName == "nE139CD"  ||
           setName == "nE139FED" || setName == "nE139HED")
        F2_D = QCD::disobs("DIS_NCE_D",x,Q,y);//Deuterium Cross section
        else
        F2_D = QCD::disobs("DIS_F2D",x,Q,y);//Deuterium structure function

        if(F2_D == 0)
        {
          if(FOR_APFELGRID)
          outstream<<std::numeric_limits<double>::signaling_NaN()<<endl;//"\t"<<std::numeric_limits<double>::signaling_NaN()<<endl;
          else
          outstream<<std::numeric_limits<double>::signaling_NaN()<<" "<<std::numeric_limits<double>::signaling_NaN()<<endl;

        }
        else
        {
          if(FOR_APFELGRID)
          outstream<<1./F2_D<<endl;//<<" "<<0<<endl; // or outstream<<2./F2_D<<" "<<0<<endl; ?!
          else
          outstream<<1./F2_D<<" "<<0<<endl;
        }
      }
    }
    outstream.close();


  }
  exit(0);
}
