
#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include <sys/time.h>
#include <map>

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>
#include "NNPDF/nnpdfdb.h"


using namespace std;
using NNPDF::FKHeader;

namespace DIS
{

  // ******************* CommonData parsing ******************************

  void parse_input(int innum, dis_param& param, std::string dbfile)
  {
    // Setup db connection
    NNPDF::IndexDB grid_db(databasePath()+dbfile, "grids");
    NNPDF::IndexDB subgrid_db(databasePath()+dbfile, "dis_subgrids");

    // Fetch number of entries
    const int entries =subgrid_db.GetNEntries();
    if (innum < 0 || innum > entries)
    {
      cerr << "Error: DIS ID ("<<innum<<") must be between 1 and "<<entries<<endl;
      exit(-1);
    }

    // Read grid information
    const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,innum,"fktarget");
    const int target = NNPDF::dbmatch(grid_db, "name", fktarget)[0];
    param.nx      =  NNPDF::dbquery<int>(grid_db,target,"nx");
    param.setname =  NNPDF::dbquery<string>(grid_db,target,"setname");

    param.gridname    = fktarget + "_" + to_string(innum) + ".subgrid";
    param.process     = NNPDF::dbquery<string>(subgrid_db,innum,"process");
    param.positivity  = NNPDF::dbquery<bool>(subgrid_db,innum,"positivity");

    // Fix positivity observables to NLO and disable TMCs
    if (param.positivity)
    {
      param.thMap["TMC"] = '0'; 
      param.thMap["PTO"] = to_string(std::min(param.evol_pto,(size_t)1)); 
      param.evol_pto = std::min(param.evol_pto,(size_t)1);

      std::cout<< "****** POSITIVITY OBSERVABLE ******"<<std::endl;
      std::cout<< "Limiting PTO to NLO, disabling TMCs"<<std::endl;
      std::cout<< "***********************************"<<std::endl;
    }

    param.commonfile = dataPath() + "commondata/DATA_" + param.setname + ".dat"; //!< Path for the commondata file
    param.sysfile    = dataPath() + "commondata/systypes/SYSTYPE_" + param.setname + "_DEFAULT.dat";    //!< Path for the SYSTYPE file
    
    stringstream desc;
    desc << "-------------------------------" <<std::endl
         << " FK_"<<param.gridname<<".dat" <<std::endl
         << "-------------------------------";    
    param.desc = desc.str();      //!< FKTable description
    
    param.ndata = 0;  //!< Number of selected datapoints
    param.xmin  = 0;  //!< Minimum x-value to be used in interpolation

    cout << "    - Gridname: "<<param.gridname<<endl;
    cout << "    - SetName: "<<param.setname<<endl;
    DisplayHR();

    return;
  }

  void set_params(dis_param const& par, NNPDF::FKHeader& FK)
  {
    FK.AddTag(NNPDF::FKHeader::BLOB, "GridDesc", par.desc);
    FK.AddTag(NNPDF::FKHeader::GRIDINFO, "SETNAME", par.setname);
    FK.AddTag(NNPDF::FKHeader::GRIDINFO, "NDATA", par.ndata);
    FK.AddTag(NNPDF::FKHeader::GRIDINFO, "HADRONIC", false);

    // Full flavourmap
    stringstream fMapHeader;
    for (int i=0; i<14; i++)
        fMapHeader << "1 ";
    FK.AddTag(FKHeader::BLOB, "FlavourMap", fMapHeader.str());

    // Set QCD parameters
    QCD::set_params(par, FK);
  }



  // ******************* FK Table computation ****************************

  void computeFK(dis_param const& par, NNPDF::CommonData const& cd, NNPDF::FKGenerator* FK)
  {
    for (size_t d=0; d<par.ndata; d++)
    {
      std::cout << d<<"/"<<cd.GetNData() <<" datapoints computed"<<std::endl;
      const double x = cd.GetKinematics(d,0);
      const double Q = sqrt(cd.GetKinematics(d,1));
      const double y = cd.GetKinematics(d,2);

      // (For the moment) skip data with Q<Q0
      if (Q < par.Q0)
        continue;

      for(size_t ix=0; ix<FK->GetNx(); ix++) 
        for(int ifl=0; ifl<14; ifl++) 
          FK->Fill(d, ix, ifl, QCD::diskernel(par.process, x, Q, y, ifl, ix) );
    }

    return;
  }

  double getXmin(const NNPDF::CommonData& g)
  {
    // Proton mass for scaling variable
    const double M2proton = pow(APFEL::GetProtonMass(),2);

    double ximin = 1.0;
    for (int i=0; i<g.GetNData(); i++)
    {
      const double x = g.GetKinematics(i,0);
      const double Q2 = g.GetKinematics(i,1);

      const double tau = 1.0 + 4.0*(M2proton/Q2)*pow(x,2.0);
      const double xi = 2.0 * x / ( 1.0 + sqrt(tau) );

      ximin = fmin(ximin, xi);
    }
    
    return ximin;
  }

  // Get the maximum scale of an applgrid
  double getQ2max(const NNPDF::CommonData& g)
  {
    double q2max = 0;
    
    for (int i=0; i<g.GetNData(); i++)
     q2max = fmax(q2max, g.GetKinematics(i,1));
    
    return q2max;
  }

    // Get the maximum scale of an applgrid
  double getQ2min(const NNPDF::CommonData& g)
  {
    double q2min = 0;
    
    for (int i=0; i<g.GetNData(); i++)
     q2min = fmin(q2min, g.GetKinematics(i,1));
    
    return q2min;
  }

}
