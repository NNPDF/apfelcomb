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

/*
double getQ2Max_CommonData(NNPDF::CommonData const& cd)
{
double q2max = 0;
for (int i=0; i<cd.GetNData(); i++)
q2max = fmax(q2max, cd.GetKinematics(i,1));
return q2max;
}*/

int main(int argc, char* argv[]) {

  if ( argc != 5 )
  {
    std::cerr << "Usage: " <<argv[0] <<"<source=app/dis/dyp> <database id> <theory id> <set_name>" <<std::endl;
    exit(-1);
  }

  // Paths to CommonData
  // APPLgrid and theory indices
  const string source = argv[1];
  const int iDB = atoi(argv[2]);
  const int iTh = atoi(argv[3]);
  const std::string setName = argv[4];
  setupDir(iTh);

  const std::string cDataPath = dataPath() + "/commondata/DATA_"+setName+".dat";
  const std::string sysTypePath = dataPath() + "/commondata/systypes/SYSTYPE_"+setName+"_DEFAULT.dat";

  const std::string F2DPath = dataPath() + "/F2D_unnormalization/F2D_"+setName+".dat";

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

  /*
  //Check kinematics
  cout<<"double x[]={";
  for (int i=0; i<cd.GetNData(); i++)
  {    x=cd.GetKinematics(i,0);
  cout<<x<<", ";
  }
  cout<<"};"<<endl;

  cout<<"double Q[]={";
  for (int i=0; i<cd.GetNData(); i++)
  {
    Q2=cd.GetKinematics(i,1);
    Q=sqrt(Q2);
    cout<<Q<<", ";
  }
  cout<<"};"<<endl;*/

  std::ofstream outstream(F2DPath.c_str());
  int Nrep=100;
  for (int i=0; i<cd.GetNData(); i++)
  {
    x=cd.GetKinematics(i,0);
    Q2=cd.GetKinematics(i,1);
    Q=sqrt(Q2);
    y=cd.GetKinematics(i,2);

    QCD::initPDF("NNPDF31_nnlo_as_0118.LHgrid",0);
    F2_D = QCD::disobs("DIS_F2D",x,Q,y);//Deuterium structure function
    outstream<<1./F2_D<<"  ";

    for(int j=1; j<=Nrep; j++)
    {
      QCD::initPDF("NNPDF31_nnlo_as_0118.LHgrid",j);
      F2_D = QCD::disobs("DIS_F2D",x,Q,y);//Deuterium structure function
      outstream<<1./F2_D<<" ";
    }
    outstream<<endl;
  }
  outstream.close();
  exit(0);
}
