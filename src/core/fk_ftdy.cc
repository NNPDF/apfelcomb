
#include "apfelcomb/fk_ftdy.h"
#include "apfelcomb/fk_qcd.h"

#include <sys/time.h>

#include <NNPDF/common.h>
#include <NNPDF/commondata.h>
#include <NNPDF/utils.h>
#include "NNPDF/nnpdfdb.h"

#include "APFEL/APFEL.h"

using namespace std;
using NNPDF::FKHeader;

namespace FTDY
{

  // *********************** FKHeader population ************************

  void set_params(QCD::qcd_param const& par, std::string const& gridname, std::string const& setname, int const& ndata, NNPDF::FKHeader& FK)
  {
    std::stringstream desc;
    desc  <<  "-----------------------------------------------------------"<<std::endl
          <<  " FK_"<<gridname<<".dat"<<std::endl
          <<  "-----------------------------------------------------------";

    FK.AddTag(FKHeader::BLOB, "GridDesc", desc.str());
    FK.AddTag(FKHeader::GRIDINFO, "SETNAME", setname);
    FK.AddTag(FKHeader::GRIDINFO, "NDATA", ndata);
    FK.AddTag(FKHeader::GRIDINFO, "HADRONIC", true);

    // Full flavourmap
    stringstream fMapHeader;
    for (int i=0; i<14; i++)
    {
      for (int i=0; i<14; i++)
        fMapHeader << "1 ";
      fMapHeader<<std::endl;
    }
    FK.AddTag(FKHeader::BLOB, "FlavourMap", fMapHeader.str());

    // Set QCD parameters
    QCD::set_params(par, FK);
  }

  // Populate FK table
  void computeGrid(QCD::qcd_param const& par, NNPDF::CommonData const& cd)
  {
    int flmap[196];
    const std::string hcxfile = "data/FTDY/"+cd.GetSetName()+".hcx";

    // Setup required directories
    stringstream theoryDir;
    theoryDir << resultsPath()<<"theory_" << par.thID<<"/subgrids/";

    APFEL::ComputeFKTables(hcxfile, theoryDir.str().c_str(), par.Q0, flmap);
  }

  // Populate FK table
  void processFK(QCD::qcd_param const& par, NNPDF::CommonData const& cd, std::string const& fkname, NNPDF::FKGenerator* fk)
  {
    stringstream theoryDir;
    theoryDir << resultsPath()<<"theory_" << par.thID<<"/subgrids/";
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