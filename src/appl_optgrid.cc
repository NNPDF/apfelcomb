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

  appl::grid *g = NULL;
  fastnlo *fg   = NULL;
  
  if (par.fnlo)
  {
    fg = new fastnlo(par.gridfile);
    g = fg->grids()[par.fnlobin];
    
  } else{ g = new appl::grid(par.gridfile); }
  

  // Initialise QCD
  QCD::initQCD(par, APP::getQ2max(g));
  QCD::initTruthGrid(par, APP::getXmin(g,false)); 
  const std::string testPDF = "NNPDF30_nlo_as_0118.LHgrid";
  cout <<endl<< "--  Calculating truth values *************************************"<<endl;
  
  // Compute with applgrid interface (at NLO)
  int pto = par.pto-1;
  if (par.ptmin == 1)
    pto = -1;
  
  const int iCheck = 11;
  vector< vector<double > > truth;
  for (size_t n=0; n<iCheck; n++)
  {
    QCD::initPDF(testPDF, n);
    vector<double> itruth  = g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto );
    truth.push_back(itruth);
  }
  
  cout <<endl<< "--  Calculating minimum grid size *************************************"<<endl;

  // Targets and xmin
  const double target = par.tgtprec;
  const double xmin = par.xmin;
  const double xtest = APP::getXmin(g,true);

  if ( xtest < 0.99*xmin || xmin > 1)
  {
    cerr << "Error: minimum x value incorrectly set in database: should be "<<xtest<<endl;
    exit(-1);
  }

  // Output to file
  stringstream historyfile;
  historyfile << "./res/opt/ID_"<<iDB<<".hist";
  ofstream histout(historyfile.str().c_str());
  histout << par.setname <<"\t XMIN: "<<xmin<<endl;
  
  int n;
  double max, avg;
  cout << "Target precision: "<<target<<endl;
  for (n=15; n<100; n+=5)
  {
    QCD::initEvolgrid(n,xmin);

    avg = 0.0;
    max = 0.0;

    for (int imem=0; imem<truth.size(); imem++)
    {
      QCD::initPDF(testPDF, imem);
      vector<double> xsec  = g->vconvolute( QCD::evolpdf_applgrid, QCD::alphas, pto );
      vector<double> diff;
      
      size_t ibin=0;
      for (size_t o=0; o<par.nbins; o++)
        if (par.mask[o])
        {
          diff.push_back(abs((truth[imem][o] - xsec[o])/truth[imem][o]));
          ibin++;
        }
      
      avg += std::accumulate(diff.begin(),diff.end(),0.0);
      max = std::max(max,*std::max_element(diff.begin(), diff.end()));

      if ( max > target) // Threshold already breached
        break;
    }
    avg /= 100*par.nbins;

    cout << "NX: "<<n<<"  AVGERR: "<<avg<<"  MAXERR: "<<max<<endl;
    histout<< n<<" \t "<<avg<<" \t "<<max<<endl;
    
    if (max < target)
      break;
  }
  
  histout.close();
  
  // Output to file
  stringstream targetfile;
  targetfile << "./res/opt/ID_"<<iDB<<".dat";
  ofstream datout(targetfile.str().c_str());
  datout << par.setname <<"\t\t XMIN: "<<xmin<<"\t\t NX: "<< n << "\t\t AVGERR: "<<avg<<"\t\t MAXERR: "<<max<<endl;
  datout.close();
  
  if (par.fnlo)
    delete fg;
  else
    delete g;
  
  cout <<endl<< "--  OptGrid Complete **********************************"<<endl;
  
  exit(0);
}

