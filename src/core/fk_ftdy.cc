
#include "apfelcomb/fk_ftdy.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>

#include "APFEL/APFEL.h"

using namespace std;
using NNPDF::FKHeader;

namespace FTDY
{

  size_t SubGrid::GetNdat() const 
  {
    return parent.GetCommonData().GetNData();
  }

  double SubGrid::GetQ2max() const
  {
    const NNPDF::CommonData& g = parent.GetCommonData();
    double q2max = 0;
    for (int i=0; i<g.GetNData(); i++)
     q2max = fmax(q2max, g.GetKinematics(i,1));
    return q2max;
  }

  double SubGrid::GetXmin() const
  {
    const NNPDF::CommonData& g = parent.GetCommonData();
    double xmin = 1.0;
    for (int i=0; i<g.GetNData(); i++)
    {
      const double y     = g.GetKinematics(i, 0);
      const double m2    = g.GetKinematics(i, 1);
      const double sshad = g.GetKinematics(i, 2);
      const double STAUdat = sqrt(m2)/sshad;
      xmin = min(xmin, STAUdat * exp(-y));
    }
    return xmin;
  }


  void SubGrid::Compute(QCD::qcd_param const& par, vector<double>& xsec) const
  {
   for (double& i : xsec) i = 0;
  }


  void SubGrid::Splash(ostream& o) const
  {
    FKSubGrid::Splash(o);
  }

  void SubGrid::Combine(QCD::qcd_param const& par, NNPDF::FKGenerator* FK) const
  {
    int flmap[196];
    const NNPDF::CommonData& cd = parent.GetCommonData();
    const std::string hcxfile = "data/FTDY/"+cd.GetSetName()+".hcx";
    const std::string path = setupDir(par.thID, "ftdy_temp_"+to_string(id)+'/');
    APFEL::ComputeFKTables(hcxfile, path.c_str(), par.Q0, flmap);

    // Read written table
    const std::string filename = path + "FK_"+parent.GetTargetName()+".dat";
    std::fstream f; f.open(filename.c_str(), std::ios::in);
    
    if (f.fail()) {
      std::cerr << "Error opening grid file: " << filename << std::endl;
      exit(-1);
    }

    std::string line;
    while(getline(f,line))
    {
      const std::vector<double> linesplit = dsplit(line);

      const int d = linesplit[0];
      const int i = linesplit[1];
      const int j = linesplit[2];

      int idx = 3;
      for (size_t k=0; k<14; k++) // loop over flavour 1
        for (size_t l=0; l<14; l++) // loop over flavour 2
        {
          FK->Fill( d-1, i, j, k, l, linesplit[idx] );
          idx=idx+1;
        }
    }

    f.close();
    remove(filename.c_str());

    return;
  }

}