#include <sys/stat.h>

#include <string>
using std::string;

#include "nlogrid.h"
#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"
using appl::grid;
using appl::igrid;

#include <bits/hhc-process.h>

#include "pdf-genlha.h"

#include "TFile.h"
#include "TVectorT.h"
//
//     renormalisation and factorisation scales from file
//     for consistency
//
#include "scales.h"
//
//grid parameters
//
static const int NQ2bins = 25, iOrderQ2 = 5;
static const double Q2low = 99., Q2up = 1.225e7;

static const int  Nxbins = 40,  iOrderx = 5;
static const double xlow = 1e-6, xup = 1.0;

static const double apramval=5.;
static const bool pdfWeight = false;
// static const bool pdfWeight = true;

static const string pdf_function = "nlojet";
static const int lowest_order = 2, nloops = 1;



// clone a histogram, scale it and write it, then delete it
void ScaleWrite(TH1D* h, double d ) { 
  std::string name = std::string("N") + std::string(h->GetName());
  TH1D* _h = (TH1D*)h->Clone(name.c_str());
  _h->Scale(d);
  _h->Write("",TObject::kOverwrite);
  delete _h;
}

//
//
//
//  Constructor & destructor
//
//
//
nlogrid::nlogrid(const std::string& inputName, const int &binNumber) : isEventBad(false), debug(false)
{
  // number of events
  numOfEvents = 0;
  numOfBadEvents = 0;

  // file to store grid
  observableName = inputName;
  fullFileName = "./output/" +  observableName + ".root";

  refDirName = "addReference";

  Directory nlogridDir(observableName);
  nlogridDir.push();

  int _nbins = -1;
  const double * _bins = NULL;
  switch( binNumber)
    {
    case 0 :
      _nbins  = nPtBins0;
      _bins   = ptBins0;
      break;
    case 1 :
      _nbins  = nPtBins1;
      _bins   = ptBins1;
      break;
    case 2 :
      _nbins  = nPtBins2;
      _bins   = ptBins2;
      break;
    case 3 :
      _nbins  = nPtBins3;
      _bins   = ptBins3;
      break;
    case 4 :
      _nbins  = nPtBins4;
      _bins   = ptBins4;
      break;
    case 5 :
      _nbins  = nPtBins5;
      _bins   = ptBins5;
      break;
    default:
      Error("Rebinning","Wrong Rapidity bin ");
      exit(-1);
      break;
    }


  //  FILE * testfile; 
  //  testfile = fopen (fullFileName.c_str(),"r");
  //  if (testfile == NULL) 
  
  // test if file exists, don't need to try to open it, 
  struct stat stFileInfo;
  if ( stat(fullFileName.c_str(),&stFileInfo) )   
    {
      cout<<"Creating new grid... "<<endl;
      
      mode=0;
      

      appl::grid::transformvar(apramval);
      
      gridObject = new appl::grid(	       
		       _nbins, _bins,
		       NQ2bins, Q2low, Q2up, iOrderQ2,    
		       Nxbins, xlow, xup, iOrderx,       
		       pdf_function, lowest_order, nloops
		       );
      gridObject->reweight(pdfWeight);
      cout<<"DONE Creating new grid TEST : "<<_nbins<<endl;
    }
  else
    {
      cout<<"Creating optimized grid... "<<endl;
      
      mode=1;

      gridObject = new appl::grid(fullFileName);
      nlogridDir.push();
      if (gridObject->isOptimised())
	{
	  delete gridObject;
	  std::cout <<"Grid is aready optimised. Quitting ..."<<std::endl;
	  exit(0); // should really throw an exception
	}
      // reseting reference histgram
      TH1D* htemp = gridObject->getReference();
      for (int i = 0; i <= htemp->GetNbinsX() + 1; i++) htemp->SetBinContent(i,0);

      gridObject->optimise();

      cout<<"DONE Creating optimized grid"<<endl;
    }

  std::cout <<"------------------------------------"<< std::endl;
  std::cout <<"order = "<<gridObject->leadingOrder()<< std::endl;
  std::cout <<"nloops = "<<gridObject->nloops()<< std::endl;
  std::cout <<" doc = "<<gridObject->getDocumentation()<< std::endl;
  std::cout <<"------------------------------------"<< std::endl;
  
#ifdef REN_REFERENCE
  // booking reference histos
  bookReferenceHistograms(_nbins, _bins);  
#endif
  //Write out grid information to screen
  //cout<<*gridObject<<endl;
  nlogridDir.pop(); 
  applpdf = new appl_pdf();
}
//
// destructor
//
nlogrid::~nlogrid()
{
  //  cout<<"\t\t\t\t nlogrid destructor"<<endl;

#ifdef REN_REFERENCE
  deleteReferenceHistograms();
#endif

  delete gridObject;
  delete applpdf;
  currentEvent.clear();

  //  cout<<" \n\n\n\n\n\n \t\t\t Calculation finished!!!! \n\n\n\n\n\n"<<endl;
}
//
//    Grid persistency
//
void nlogrid::writeGrid(long long& nRuns)
{
  gridObject->run() = nRuns;

  Directory obs(observableName);
  obs.push();
  gridObject->Write(fullFileName);

  obs.pop();

#ifdef REN_REFERENCE
      writeReferenceHistograms();
#endif

  string goMode = (mode == 0 ? "Non-o": "O");
  goMode += "ptimised";

  std::cout<<"\tGridObject ( "<<goMode<<" ) saved after "<<nRuns
	   <<" number of NLOJET events. Accepted = "<<numOfEvents<<" ( "<<numOfBadEvents<<" bad )"
	   <<"\n Weightgrid File = "<<fullFileName<<" ."
	   <<std::endl;
}
//
//
//
void nlogrid::fill(const double &x1,
		   const double &x2,
		   const double &Q2,
		   const double &obs,
		   const nlo::amplitude_hhc& amp,
		   nlo::pdf_and_coupling_hhc* pdf)
{
  if (debug) std::cout << __PRETTY_FUNCTION__ << std::endl;
  //  if ( !gridObject->obsInRange( obs ) ) return;

  if ( isEventBad ) return; 

  double xs = conversion( amp( pdf, Q2, Q2, pb_fac) );
  if ( xs != xs )
    {
      isEventBad = true;
      return;
    }
  //
  //
  //
  appl_event _ce;
  _ce.x1 = x1;
  _ce.x2 = x2;
  _ce.Q2 = Q2;
  _ce.obs = obs;
  _ce.xSection = xs;
#ifdef REN_REFERENCE
  for(int ir = 0; ir < Nscales; ir++) _ce.xSectionScale[ir] = conversion ( amp(pdf, mur[ir]*mur[ir]*Q2, muf[ir]*muf[ir]*Q2, pb_fac ) );
#endif
  if (getMode() == 0 )
    for (int iSubProcess = 0; iSubProcess < nSubProcesses; iSubProcess++) _ce.weights[iSubProcess] =  1.0;
  else
    for (int iSubProcess = 0; iSubProcess < nSubProcesses; iSubProcess++) _ce.weights[iSubProcess] =   amp(applpdf, Q2, Q2, pb_fac )[iSubProcess];
  _ce.processOrder = amp.contrib()==0 ? 0 : 1 ;
  //
  //
  //
  currentEvent.push_back( _ce ); 

  return;
  
}
void nlogrid::endOfEvent()
{
  if (debug) std::cout << __PRETTY_FUNCTION__ << std::endl;
  
  if (currentEvent.size() == 0) return;
  
  numOfEvents++;
  if (isEventBad)
    {
      numOfBadEvents++;
      return;
    }
  
  for(std::vector<appl_event>::iterator i = currentEvent.begin(); i != currentEvent.end(); i++)
    {
      if ( getMode() == 0 )
	gridObject->fill_phasespace( (*i).x1, (*i).x2, (*i).Q2, (*i).obs, (*i).weights, (*i).processOrder);
      else
	gridObject->fill((*i).x1, (*i).x2, (*i).Q2, (*i).obs, (*i).weights, (*i).processOrder);
      
      fillReferenceHistograms( (*i).obs, (*i).xSection, (*i).xSectionScale );
    }
  
  
  currentEvent.resize(0);
  return;
}
//
//
//
#ifdef REN_REFERENCE

void nlogrid::bookReferenceHistograms(const int & _nbins, const double * _bins)
{
  Directory refDir(refDirName);
  refDir.push();
  
  char histname[3000]; char htit[3000];

  for(int ir = 0; ir <  Nscales; ir++)
    { // loop ren scal variations
      sprintf(histname,"referenceScale_%d", ir);
      sprintf(htit,"addRefHist. scales  #mu_{R} = %3.2f, #mu_{F} = %3.2f", mur[ir], muf[ir]);
      referenceScale[ir] = new TH1D(histname, htit, _nbins, _bins );
      referenceScale[ir]->Sumw2();
    }
  refDir.pop();

}                   // bookReferenceHistograms

void nlogrid::deleteReferenceHistograms()
{
  
  for(int ir = 0; ir < Nscales; ir++) 
    { 
      delete referenceScale[ir];
      referenceScale[ir] = NULL;
    }
}

void nlogrid::writeReferenceHistograms()
{
  Directory obs(observableName);
  obs.push();
  
  //cout << "nlogrid::writeReferenceHistograms() \t\t pwd=" << gDirectory->GetName() << endl;
  
  TFile myfile(fullFileName.c_str(),"UPDATE");
  
  // write the scales to the histogram;
  
  TVectorT<double> mur_v(Nscales);
  TVectorT<double> muf_v(Nscales);
  
  for ( int i=0 ; i<Nscales ; i++ ) 
    { 
      mur_v(i) = mur[i];
      muf_v(i) = muf[i];
    }
  
  //  mur_v.Write("mu_r");
  //  muf_v.Write("mu_f");
  
  mur_v.Write("mu_r");
  muf_v.Write("mu_f");

  Directory refDir(refDirName);
  refDir.push();
  //  cout << "nlogrid::writeReferenceHistograms() \t\t pwd=" << gDirectory->GetName() << endl;  
  
  double invNruns = 1;
  if ( gridObject->run()!=0 ) invNruns = 1/double(gridObject->run());

  
  for(int ir=0; ir< Nscales; ir++) ScaleWrite(referenceScale[ir], invNruns); 
  
  refDir.pop();
  myfile.Close();
  
  TFile refFile("./output/Reference.root","UPDATE");
  TDirectory* isDir = refFile.GetDirectory(observableName.c_str());
  if ( isDir == NULL ) refFile.mkdir(observableName.c_str());
  refFile.cd( observableName.c_str() );
  
  for(int ir=0; ir< Nscales; ir++) ScaleWrite(referenceScale[ir], invNruns);
  ScaleWrite(gridObject->getReference(), invNruns);
  refFile.Close();
  
  obs.pop();
}

#endif 

void nlogrid::fillReferenceHistograms( double &obs, 
				       double &xs,
				       double* xsScales
				       )
{
  
  bool debug=false;
  double binwidth = gridObject->deltaobs(gridObject->obsbin(obs));
  
  gridObject->getReference()->Fill(obs, xs/binwidth);
//  gridObject->getReferenceSquared()->Fill(obs, xs*xs/binwidth);
  
#ifdef REN_REFERENCE
  for(int ir = 0; ir < Nscales; ir++) referenceScale[ir]->Fill(obs, xsScales[ir]/binwidth);
#endif     

  return;
}




