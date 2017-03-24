
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
    NNPDF::IndexDB db(databasePath()+dbfile, "sets");

    // Fetch number of entries
    const int entries =db.GetNEntries();
    if (innum < 0 || innum > entries)
    {
      cerr << "Error: DIS ID ("<<innum<<") must be between 1 and "<<entries<<endl;
      exit(-1);
    }

    param.setname   = NNPDF::dbquery<string>(db,innum,"setname");
    param.gridname  = NNPDF::dbquery<string>(db,innum,"gridname");
    param.process   = NNPDF::dbquery<string>(db,innum,"process");

    param.positivity   = NNPDF::dbquery<bool>(db,innum,"positivity");

    // Get common grids
    vector<int> commonGrids = NNPDF::dbmatch(db, "setname", param.setname);
    for ( auto i : commonGrids)
      param.inventory.push_back(NNPDF::dbquery<string>(db,i,"gridname"));

    // Fix positivity observables to NLO and disable TMCs
    if (param.positivity)
    {
      param.thMap["TMC"] = '0'; 
      param.thMap["PTO"] = to_string(std::min(param.evol_pto,(size_t)1)); 
      param.evol_pto = std::min(param.evol_pto,(size_t)1);
      param.thMap["XIF"] = '1';
      param.thMap["XIR"] = '1';

      std::cout<< "****** POSITIVITY OBSERVABLE ******"<<std::endl;
      std::cout<< "Limiting PTO to NLO, disabling TMCs"<<std::endl;
      std::cout<< "Setting xi_F and xi_R to 1         "<<std::endl;
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
    param.nx    = 0;  //!< Number of interpolation grid x-points to be used
    param.xmin  = 0;  //!< Minimum x-value to be used in interpolation

    cout << "    - Gridname: "<<param.gridname<<endl;
    cout << "    - SetName: "<<param.setname<<endl;
    cout << "    - Common grids: "<<param.inventory.size()<<endl;
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
