// optgrid.cc
// Code to optimise x-grid nX
//
// nph  09/14

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <numeric>
#include <sys/time.h>

#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_qcd.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  if (argc!=3)
  {
    cout << "Usage: appl_comb <database id> <theory id>"<<endl;
    exit(1);
  }

  Splash();
  
  // APPLgrid and theory indices
  const int iDB = atoi(argv[1]);
  const int iTh = atoi(argv[2]);

  // Parse parameters
  APP::appl_param par;
  QCD::parse_input(iTh, par);
  APP::parse_input(iDB, par);
  par.datamap = APP::parse_map(par.common_subgrids, iDB);
  APP::grid sourceGrid(par);


  const std::string commonfile = dataPath() + "commondata/DATA_" + par.setname + ".dat"; 
  const std::string sysfile    = dataPath() + "commondata/systypes/SYSTYPE_" + par.setname + "_DEFAULT.dat";  
  const NNPDF::CommonData cd = NNPDF::CommonData::ReadFile(commonfile, sysfile);
  
  // Initialise QCD
  QCD::initQCD(par, APP::getQ2max(sourceGrid.g));
  QCD::initTruthGrid(par, APP::getXmin(sourceGrid.g,false)); 
  const std::string testPDF = "NNPDF30_nlo_as_0118.LHgrid";
  cout <<endl<< "--  Calculating truth values *************************************"<<endl;
  
  // Compute with applgrid interface (at NLO)
  int pto = par.pto-1;
  if (par.ptmin == 1)
    pto = -1;
  
  const int iCheck = 100;
  vector< vector<double > > truth;
  for (size_t n=0; n<iCheck; n++)
  {
    QCD::initPDF(testPDF, n);
    vector<double> itruth  = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto );
    truth.push_back(itruth);
  }
  
  cout <<endl<< "--  Calculating minimum grid size *************************************"<<endl;

  // // Targets and xmin
  const double xmin = APP::parse_xmin(par.common_subgrids);

  // // Output to file
  stringstream historyfile;
  historyfile << "./res/opt/ID_"<<iDB<<".hist";
  ofstream histout(historyfile.str().c_str());
  histout << par.setname <<"\t XMIN: "<<xmin<<endl;
  
  int n, maxn;
  double max, max0;
  for (n=15; n<100; n+=5)
  {
    QCD::initEvolgrid(n,xmin);

    max0 = 0.0;
    max = 0.0;
    maxn = 0;
    for (int imem=0; imem<truth.size(); imem++)
    {
      QCD::initPDF(testPDF, imem);
      const vector<double> xsec  = sourceGrid.g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto );
      vector<double> diff;
      
      size_t ibin=0;
      for (size_t o=0; o<sourceGrid.g->Nobs(); o++)
        if (par.mask[o])
        {
          diff.push_back(abs((truth[imem][o] - xsec[o])/truth[imem][o]));
          ibin++;
        }

      for (int i=0; i<par.ndata; i++)
      {
        const double truthpred = truth[imem][par.maskmap[i]];
        const double testpred = xsec[par.maskmap[i]];
        const double rel_err = abs((truthpred-testpred)/testpred);
        const double targetPrec = APP::computeTargetPrecision(par.datamap[i], cd);
        diff.push_back(rel_err/targetPrec);
      }
      
      if (imem == 0)
        max0 = *std::max_element(diff.begin(), diff.end());
      max = std::max(max,*std::max_element(diff.begin(), diff.end()));

      if ( max > 1.0) // Threshold already breached
      {
        maxn = imem;
        break;
      }
    }

    cout << "NX: "<<n<<"  MAX_0: "<<max0<<"  MAXERR: "<<max<< " ("<<maxn<<")"<<endl;
    histout<< n<<" \t "<<max0<<" \t "<<max<< " ("<<maxn<<")"<<endl;
    
    if (max < 1.0)
      break;
  }
  
  histout.close();
  
  // Output to file
  stringstream targetfile;
  targetfile << "./res/opt/ID_"<<iDB<<".dat";
  ofstream datout(targetfile.str().c_str());
  datout << par.setname <<"\t\t XMIN: "<<xmin<<"\t\t NX: "<< n << "\t\t MAX_0: "<<max0<<"\t\t MAXERR: "<<max<<endl;
  datout.close();
  
  cout <<endl<< "--  OptGrid Complete **********************************"<<endl;
  
  exit(0);
}

