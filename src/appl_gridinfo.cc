// appl_gridinfo.cc
// Examine an applgrid - returns information
// required by convolution
//
// n.p.hartland@ed.ac.uk - 03/12

#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

#include "fk_appl.h"
#include "appl_grid/appl_grid.h"

#include "LHAPDF/LHAPDF.h"

using namespace std;

int main(int argc, char* argv[]) {
  
  if (argc!=2)
  {
    cout << "Invalid Parameters:"<<endl;
    cout <<" Syntax: ./appl_gridinfo <applgrid file> [CommonData file]"<<endl;
    exit(1);
  }

  Splash();
  
  string filename(argv[1]);
  bool fnlo=false;
  size_t fnlobin=4;
  
  // Init applgrid
  appl::grid *g = NULL;
  fastnlo *fg = NULL;
  
  if (fnlo)
  {
    fg = new fastnlo(filename);
    g = fg->grids()[fnlobin];
    cout << "FASTNLO GRID: "<<fnlobin<<endl;    
  } else{ g = new appl::grid(filename); }
  
  cout <<"Applgrid information:"<<endl;
  cout <<"Transform: "<<g->getTransform()<<endl;
  cout <<"GenPDF: "<<g->getGenpdf()<<endl;
  cout <<"SubProc: "<<g->subProcesses(0)<<endl;
  cout <<"Nbins: "<<g->Nobs()<<endl;
  cout << "Symmetric: " << (g->isSymmetric() ? "Y" : "N") <<endl;
  cout << "xMin: " << APP::getXmin(g, true) << endl;


  // Cleanup
  if (fnlo)
  {   delete fg; }
  else
  {   delete g;  }
  
  return 0;
}
