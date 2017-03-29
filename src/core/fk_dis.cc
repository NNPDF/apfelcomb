
#include "apfelcomb/fk_dis.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>
#include "NNPDF/nnpdfdb.h"

using namespace std;

namespace DIS
{

  size_t SubGrid::GetNdat()  const 
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
    const double M2proton = 0.938; // GeV - should ideally be taken from APFEL

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

  void SubGrid::Compute(QCD::qcd_param const& par, vector<double>& xsec) const
  {
    const time_point t1 = std::chrono::system_clock::now();
    const NNPDF::CommonData& cd = parent.GetCommonData();
    for (int i=0; i<cd.GetNData(); i++)
    {
      const double x = cd.GetKinematics(i,0);
      const double Q = sqrt(cd.GetKinematics(i,1));
      const double y = cd.GetKinematics(i,2);

      if (min(par.xiF, par.xiR)*Q >= par.Q0)
        xsec[i] = nrmdat*QCD::disobs(process, x, Q, y);
      else
        xsec[i] = 0;
      StatusUpdate(t1, ((double)i+1)/(double)cd.GetNData(), cout );
    }
  }

  void SubGrid::Splash(ostream& o) const
  {
    FKSubGrid::Splash(o);
  }

  void SubGrid::Combine(QCD::qcd_param const& par, NNPDF::FKGenerator* FK) const
  {
    const time_point t1 = std::chrono::system_clock::now();
    const NNPDF::CommonData& cd = parent.GetCommonData();
    for (size_t d=0; d<cd.GetNData(); d++)
    {
      const double x = cd.GetKinematics(d,0);
      const double Q = sqrt(cd.GetKinematics(d,1));
      const double y = cd.GetKinematics(d,2);

      if (min(par.xiF, par.xiR)*Q < par.Q0)
        continue;

      for(size_t ix=0; ix<FK->GetNx(); ix++) 
        for(int ifl=0; ifl<14; ifl++) 
          FK->Fill(d, ix, ifl, nrmdat*QCD::diskernel(process, x, Q, y, ifl, ix) );
      StatusUpdate(t1, ((double)d+1)/(double)cd.GetNData(), cout );
    }
    return;
  }
}
