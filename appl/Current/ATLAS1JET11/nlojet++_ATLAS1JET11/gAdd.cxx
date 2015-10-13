#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>

#include <sys/stat.h>

#include "appl_grid/appl_grid.h"
#include "appl_grid/Directory.h"


#include "appl_grid/appl_timer.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TFile.h>
#include <TPad.h>
#include "TSystem.h"

#include "LHAPDF.h"
#include "hoppet_v1.h"

#define DBG false


#include "scales.h"
#include "Numbers.hh"


using namespace appl;

Bool_t checkRun(const TH1D*& lim, const TH1D*& ref);

// extern double Escale;

int main(int argc, char** argv) { 
  //
  // run parameters
  // gAdd dirName gridName  flagPS nRuns subRange
  // gridName is supposed to be under ./dirName/runI/gridName.root
  //
  // ex. : gAdd output eta1 10 1 1
  // will produce output/int/gridName.root
  //
  //  castor modification
  //  ./gAdd PATH NAME nGRIDS FLAG
  //


  if ( argc < 6 ) return -1;

  TString gridDirectory(argv[1]);
  TString gridName(argv[2]);
  bool onlyPhaseSpace = (bool)atoi(argv[3]);
  int nRuns = atoi(argv[4]);
  int iSubSet = atoi(argv[5]);

  TString gridType = "grid";
  if (onlyPhaseSpace) gridType = "phasespace";

  grid* g[nRuns];
  TFile* gf[nRuns];
  FileStat_t filestatBuffer;


  std::cout <<"----------------------------------------------------------------------------------"<<std::endl;
  TFile* fileLimits = TFile::Open("obsLimits.root","read");
  if ( fileLimits == NULL)
   {
     std::cout << "File with limits is not found !!!!!!!!!!!!! " << std::endl;
     exit(-1);
   }

  fileLimits->ls();

  TString limitsHistName("");
  limitsHistName.Form("%s_Limits", gridName.Data());
  const TH1D* xsLimits = (TH1D*)fileLimits->Get( limitsHistName.Data() );
  if (xsLimits == NULL ) 
  {
   std::cout <<"Problems with getting limits histogram "<< limitsHistName << std::endl;
   exit(-1);
  }
  xsLimits->Print();

  std::cout <<"----------------------------------------------------------------------------------"<<std::endl;


  //
  //  opening files
  //
  
  for (int i = 0; i < nRuns; i++)
    {

      g[i]  = NULL;
      gf[i] = NULL;
      
      if ( iSubSet > 0 )
	{ 
	  if ( i <  ( iSubSet - 1 ) * 100 ) continue;
	  if ( i >= ( iSubSet ) * 100 )     continue;
	}

      
      TString fullGridName = "";
      //      fullGridName.Form("./%s/run%d/%s.root", 
      fullGridName.Form("%s/output.%s.%d/%s.root", 
			gridDirectory.Data(),
			gridType.Data(),
			i+1,
			gridName.Data());
      //      struct stat stfileinfo;


      TFile::SetOpenTimeout(10000);

      Int_t fileStatus = gSystem->GetPathInfo( fullGridName, filestatBuffer );
      
      if ( ( fileStatus == 0 ) && ( filestatBuffer.fSize != 0 ) )
	{	  //std::cout << "\n \t\t ****** Opening with grid " << i <<" "<<fullGridName<< " ****** \n"<< std::endl; 
	  TFileOpenHandle* fh = TFile::AsyncOpen(fullGridName, "READ NET TIMEOUT=15");
	  std::cout <<"opt = "<< fh->GetOpt() <<" m = "<<fh->Matches(fullGridName)<< std::endl;
	  gf[i] = TFile::Open(fh);
	  if ( gf[i] == 0 ) continue;
	  cout <<" --- file exists   "<<fullGridName<< endl;	

	  const TH1D* refHist = (TH1D*)gf[i]->Get("grid/reference");
	  
	  if ( !checkRun(xsLimits, refHist))
	    {
	      if (g[i]) delete g[i];
	      if (gf[i]) delete gf[i];
	      g[i] = NULL;
	      gf[i] = NULL;
	      cout <<" --- run is BAD !!! "<<fullGridName<< endl;
	      continue;
	    }
	  
	  
	  cout <<" --- run is good !!! "<<fullGridName<< endl;	

	  g[i] = new grid(fullGridName.Data());
	  bool gridPhaseSpace =  !g[i]->isOptimised();
	  if ( gridPhaseSpace == onlyPhaseSpace )
	    {
	      cout <<"MATCH "<< gridPhaseSpace <<" "<<onlyPhaseSpace<< endl;
	      g [i]->trim();
	      //gf[i] = TFile::Open(fullGridName, "read");
	    }
	  else
	    {
	      cout <<"NO MATCH "<< gridPhaseSpace <<" "<<onlyPhaseSpace<< endl;
	      if (g[i]) delete g[i];
	      if (gf[i]) delete gf[i];
	      g[i] = NULL;
	      gf[i] = NULL;
	    }
	  cout <<" ---------------------------------------  "<< endl;	
	}
      else
	{
	  std::cout <<"\t\t ZOMBI FILE = "<<fullGridName<< std::endl;
	}
    }

  //
  //  setting the first grid
  //

  int iGridStart =  0;
  if (iSubSet > 0) iGridStart = ( iSubSet - 1 ) * 100;

  if ( NULL == g[iGridStart] )
    for (int i = 0; i < nRuns; i++) 
      {
	
	if ( iSubSet > 0 )
	  { 
	    if ( i <  ( iSubSet - 1 ) * 100 ) continue;
	    if ( i >= ( iSubSet ) * 100 )     continue;
	  }
	if (i == iGridStart ) continue; 

	if (NULL != g[i])
	  {
	    g[iGridStart] = g[i]; gf[iGridStart] = gf[i];
	    g[i] = NULL; gf[i] = NULL;
	    break;
	  }

      }

  long long startRuns( g[iGridStart]->run() );
  long long nActualRuns(startRuns);

  TString histName = "";

  TH1D* reference ;
  TH1D* referenceScale[Nscales];

  reference = (TH1D*)gf[iGridStart]->Get("grid/reference");	
  reference->Scale( startRuns );
  
  for ( int j = 0 ; j < Nscales ; j++ ) 
    {
      histName.Form("addReference/NreferenceScale_%d", j); 
      referenceScale[j] = (TH1D*)gf[iGridStart]->Get( histName );
      referenceScale[j]->Scale( startRuns );
    }

  //
  //  addition
  //

  
  for (int i = 0; i < nRuns; i++) 
  {
    if ( iSubSet > 0 )
      { 
	if ( i <  ( iSubSet - 1 ) * 100 ) continue;
	if ( i >= ( iSubSet ) * 100 )     continue;
      }
    
    if ( i == iGridStart ) continue; 

    if (g[i]) 
      std::cout <<" -- adding -- "<<g[i]<<"  from run"<<i+1<<"  isOpt = "<<g[i]->isOptimised()<<" to "<<g[iGridStart];
    else
      std::cout <<" -- adding -- ZERO ---"<<std::endl;;


    if (  g[i] != NULL  ) 
      {
	*(g[iGridStart]) += *(g[i]);
	nActualRuns += g[i]->run();


	long long currentRuns = 0;
	currentRuns = g[i]->run();
	reference->Add((TH1D*)gf[i]->Get("grid/reference"), currentRuns);

	for ( int j = 0 ; j < Nscales ; j++ ) 
	  {
	    histName.Form("addReference/NreferenceScale_%d", j); 
	    referenceScale[j]->Add((TH1D*)gf[i]->Get(histName), currentRuns);
	  }
         std::cout <<" \t\t -- done -- "<< std::endl;
      }
    
  }

  TString outFileName = "";
  if (iSubSet > 0 ) 
    outFileName.Form( "output/interm/output.%s.%d/%s.root", gridType.Data(), iSubSet, gridName.Data() );
  else
    outFileName.Form( "output/%s.root", gridName.Data() );

  if ( g[iGridStart] )
    {
      g[iGridStart]->Write( outFileName.Data() );
      cout <<"writing grid done. "<< endl;

      TFile* fout = new TFile(outFileName, "update");
      gDirectory->cd("grid");
      reference->Scale(1./nActualRuns);
      reference->Write("reference",TObject::kOverwrite);
      fout->cd();

      Directory addref("addReference");
      addref.push();

      for ( int i = 0 ; i < Nscales ; i++  ) 
	{
	  referenceScale[i]->Scale(1./nActualRuns);
	  referenceScale[i]->Write("", TObject::kOverwrite);
	}

      addref.pop();
    
      fout->Write();
      fout->Close();
      
      cout <<"writing reference histograms done. "<< endl;
    }
  else
    {
      cout <<"writing is impossible due to zero file. "<< endl;
      exit( 0 );
    }
  

  
    for (int i = 0; i < nRuns; i++)
      {
	if ( iSubSet > 0 )
	  { 
	    if ( i <  ( iSubSet - 1 ) * 100 ) continue;
	    if ( i >= ( iSubSet ) * 100 )     continue;
	  }
	
	if ( g [i] ) delete g [i];
	if ( gf[i] ) delete gf[i];
      }
    
    return 0;
}








Bool_t checkRun(const TH1D*& limits, const TH1D*& hist)
{

  if ( isbad(hist->GetSumOfWeights()) ) return false;

  for (Int_t iBin = 1; iBin <= limits->GetNbinsX(); iBin++)
    {
      Double_t median = limits->GetBinContent(iBin);
      Double_t sigma  = limits->GetBinError(iBin);
      Double_t val    = hist  ->GetBinContent(iBin);
      if ( ( val < median - 10.*sigma ) || ( val > median + 10.*sigma ) ) return false;
    }

  return true;
}







