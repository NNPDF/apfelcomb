
#include "fk_ftdy.h"
#include "fk_qcd.h"
#include "nnpdfdb.h"

#include <sys/time.h>

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>
#include <NNPDF/utils.h>

#include "APFEL/APFEL.h"

using namespace std;

namespace FTDY
{

  // Populate FK table
  void computeGrid(QCD::qcd_param const& par, NNPDF::CommonData const& cd)
  {
    int flmap[196];
    const std::string hcxfile = "data/FTDY/"+cd.GetSetName()+".hcx";

    // Setup required directories
    stringstream theoryDir;
    theoryDir << resultsPath()<<"theory_" << par.thID<<"/apfelcomb/" << cd.GetSetName() <<"/";

    APFEL::ComputeFKTables(hcxfile, theoryDir.str().c_str(), par.Q0, flmap);
  }

  // Populate FK table
  void processFK(QCD::qcd_param const& par, NNPDF::CommonData const& cd, std::string const& fkname, NNPDF::FKGenerator* fk)
  {
    stringstream theoryDir;
    theoryDir << resultsPath()<<"theory_" << par.thID<<"/apfelcomb/" << cd.GetSetName() <<"/";
    const std::string filename = theoryDir.str() + "FK_"+fkname+".dat";

    // Read written table
    std::fstream f; f.open(filename.c_str(), std::ios::in);
    
    if (f.fail()) {
      std::cerr << "Error opening grid file: " << filename << std::endl;
      exit(-1);
    }

    std::string line;
    while(getline(f,line))
    {
      std::vector<NNPDF::real> linesplit;
      NNPDF::rsplit(linesplit,line);

      const int d = linesplit[0];
      const int i = linesplit[1];
      const int j = linesplit[2];

      int idx = 3;
      for (size_t k=0; k<14; k++) // loop over flavour 1
        for (size_t l=0; l<14; l++) // loop over flavour 2
        {
          fk->Fill( d-1, i, j, k, l, linesplit[idx] );
          idx=idx+1;
        }
    }

    f.close();
    remove(filename.c_str());
  }


  // Get minimum x value
  double getXmin(const NNPDF::CommonData& g)
  {
    double xmin = 1.0;
    
    for (int i=0; i<g.GetNData(); i++)
    {
      const double y     = g.GetKinematics(i, 0);
      const double m2    = g.GetKinematics(i, 1);
      const double sshad = g.GetKinematics(i, 2);
      
      const double STAUdat = sqrt(m2)/sshad;
    
      xmin = fmin(xmin, STAUdat * exp(-y));
    }
    
    return xmin;
  }

  // Get the maximum scale of an applgrid
  double getQ2max(const NNPDF::CommonData& g)
  {
    double q2max = 0;
    
    for (int i=0; i<g.GetNData(); i++)
     q2max = fmax(q2max, g.GetKinematics(i,1));
    
    return q2max;
  }

}