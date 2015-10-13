//
//
//
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//------ DON'T TOUCH THIS PART! ------//------ DON'T TOUCH THIS PART! ------
//
//
//
#include <bits/hhc-phasespace.h>
#include <bits/hhc-process.h>
#include <bits/hhc-jetfunc.h>
//#include <nlo++-module_add.h>
//----- used namespaces -----
using namespace nlo;
using namespace std;
//----- declaration of the user defined functons -----
void inputfunc(unsigned int&, unsigned int&, unsigned int&);
void psinput(phasespace_hhc *, double&);
user_base_hhc * userfunc();

//typedef unsigned long int (*module_add_type)(bool, const list<basic_string<char> >&, const basic_string<char>&);
//extern  module_add_type module_add;


//----- array of the symbols symbols -----
extern "C"{
struct { 
  const char *name;
  void *address;
} user_defined_functions[] = 
  {
    //   process index: hhc for hadron-hadron --> jets
    {"procindex", (void *) "hhc"},
    
    //   input function 
    {"inputfunc", (void *) inputfunc},
    
    //   phase space input function 
    {"psinput", (void *) psinput},
    
    //   user defined functions
    {"userfunc",  (void *) userfunc},
    
    //   module to generate the readable result
    //    {"main_module_add", (void *) module_add},
    
    //  end of the list
    {0, 0}
  };
};
//
//
//------ USER DEFINED PART STARTS HERE ------
//------ USER DEFINED PART STARTS HERE ------
//------ USER DEFINED PART STARTS HERE ------
//
//

#include <algorithm>
#include <kT_clus.h>

#include "pdf-genlha.h"

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <bits/phys-cone_seedless.h>

#include "nlogrid.h"
#include "fjClustering.h"

#include <stdlib.h>
//
//
//
const int      debug           = 0;
const long int saveAfterEvents = 100000000;
const long int maxNumberEvents = 500000000;

//
//   user class
//
class UserHHC : public user1h_hhc
{
public:
  UserHHC();
  ~UserHHC();
  void deleteObjects();
  void operations_at_the_end_of_event(); 

  long long eventNb;                    // number of event
  
  void initfunc(unsigned int);
  void userfunc(const event_hhc&, const amplitude_hhc&);

  bool isGoodEvent(){return goodEvent;};

  static const int numberOfJets = 1U;                       // MIN number of jets in event
  
  static const double sqrts     = 7000.;                    // 7 TeV initial energy
  

  
private:
  
  fjClustering* jetclus[nRadius];                           // jet clustering algorithm
  void clusterJets(const event_hhc&);
  
  pdf_and_coupling_hhc * pdf;                               // LHA
  
  double x1;                                                // initial parton1 momentum fraction
  double x2;                                                // initial parton2 momentum fraction
  double SCALE2[nRadius];                                   // highest PT^2 of jets in event
  bool goodEvent;
  //
  bool eventSelected[nRadius][nGrids];                      // event has passed the cuts
  nlogrid* mygrid[nRadius][nGrids];                         // grid for coeficients
  
  typedef lorentzvector<double> _Lv;
  bounded_vector<_Lv> cj[nRadius], jets[nRadius], pj[nRadius][nGrids];     // the jet structure in lab. frame
  bounded_vector<unsigned int> jet;

  struct pT_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.perp2() > p2.perp2();
    }
  };
  struct E_sort {
    bool operator()(const _Lv& p1, const _Lv& p2) const {
      return p1.T() > p2.T();
    }
  };
};
//
//  destructor
//
UserHHC::~UserHHC() 
{
  deleteObjects();
}
void UserHHC::deleteObjects()
{
  for (int iRadius = 0; iRadius < nRadius ; iRadius++) delete jetclus[iRadius];
  delete pdf;

  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
      if (mygrid[iRadius][iGrid]) delete  mygrid[iRadius][iGrid];

  return;
}
//
//  constructor
//
UserHHC::UserHHC() : pdf(0)  
{
  cout<<"NLOJET++ FillGrid started..."<<endl;
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    jetclus[iRadius] = new fjClustering(fastjet::antikt_algorithm, 
					jetSizes[iRadius], 
					fastjet::E_scheme, 
					fastjet::Best);
  
  // create grid structure
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
      {
	TString gridName = "";
	gridName.Form("ATLAS1JET11-R%02d-rap-%02d-%02d",
		      (int)(10*jetSizes[iRadius]+0.5),
		      (int)(10*etaBinsLow [iGrid]+0.5),
		      (int)(10*etaBinsHigh[iGrid]+0.5) );

	std::cout <<"creating new grid : "<<gridName<<".root"<< std::endl;
	mygrid[iRadius][iGrid] = new nlogrid( gridName.Data(), iGrid );
	mygrid[iRadius][iGrid]->getGridObject()->setCMSScale( sqrts );
      }
  //Zero the counter for number of events
  eventNb = 0;
  if (debug) 
    {
      cout<<"UserHHC::UserHHC() \t\t gDirectory = ";
      //      gDirectory->pwd();
    }
}
//
//
//
user_base_hhc* userfunc() {return new UserHHC;}
//
//
//
void inputfunc(unsigned int& nj, unsigned int& nu, unsigned int& nd)
{
  nj = UserHHC::numberOfJets;
  // number of UP quarks
  nu = 2U;
  // number of DOWN quarks
  nd = 3U;
}

void psinput(phasespace_hhc *ps, double& s)
{
  //  total c.m. energy square
  s = UserHHC::sqrts * UserHHC::sqrts;                  // unit:GeV
  
  //   You can use your own phase generator. 
  //   Here we use the default.
  ps = 0;
}



void UserHHC::initfunc(unsigned int)
{
  pdf = new pdf_genlha("PDFsets/cteq66.LHgrid",0);
//  pdf = new pdf_genlha("PDFsets/abm11_5n_nlo.LHgrid",0);
}

void UserHHC::clusterJets(const event_hhc& p)
{
// //   std::cout << "UserHHC::clusterJets(const event_hhc& p)" << std::endl;
// //   std::cout<< p <<std::endl;
// //   std::cout << "p[1] = " << p[1].perp() <<" p2 = "<<p[2].perp()<< std::endl;
  //
  // initialisation
  //
  x1 = x2 = -1.;
  goodEvent = false;

  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    {
      SCALE2[iRadius] = 0.;
      cj[iRadius].resize(1,0);
      jets[iRadius].resize(1,0);
      for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
	{
	  eventSelected[iRadius][iGrid] = false;
	  pj[iRadius][iGrid].resize(1,0);
	}
    }
  //
  //----- do the cluster analysis-----
  //
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    {
      cj[iRadius]=(*jetclus[iRadius])(p);                          
      
      if (debug)
	{
	  std::cout<<"\n"<<std::endl;
	  for(int jet = 1; jet <= cj[iRadius].upper(); jet++)
	    std::cout <<" radius = "<<jetSizes[iRadius]<<" jet "
		      <<jet<<" "<<cj[iRadius][jet]<<" pt = "<< cj[iRadius][jet].perp()<<std::endl;
	}

      for (int iJet = 1; iJet <=  cj[iRadius].upper(); iJet++)
	if ( ( std::abs(cj[iRadius][iJet].rapidity()) <= etaAccept ) && (std::abs(cj[iRadius][iJet].perp()) >= pT_min ) ) jets[iRadius].push_back(cj[iRadius][iJet]);
      
      std::sort( jets[iRadius].begin(), jets[iRadius].end(), pT_sort() );
    }

  //
  //  SELECT ONLY INTERESTING JETS
  //
  
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    for(int iJet = 1; iJet <= jets[iRadius].upper(); iJet++)
      {
	double jetRap = std::abs( jets[iRadius][iJet].rapidity() );
	for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
	  if ( ( jetRap  <  etaBinsHigh[iGrid] ) && ( jetRap >= etaBinsLow [iGrid] ) ) pj[iRadius][iGrid].push_back( jets[iRadius][iJet] );
      }
  
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
      if ( pj[iRadius][iGrid].upper() >= UserHHC::numberOfJets )  // there is at least one JET
	{
	  eventSelected[iRadius][iGrid] = true;
	  goodEvent = true;
	  // std::sort( pj[iRadius][iGrid].begin(), pj[iRadius][iGrid].end(), pT_sort() );
	}
 
  
  if ( this->isGoodEvent() )
    {
      x1 = p[-1].T()/(0.5*sqrts);
      x2 =  p[0].T()/(0.5*sqrts);
      for (int iRadius = 0; iRadius < nRadius ; iRadius++)
	if  ( jets[iRadius].upper() >= UserHHC::numberOfJets ) SCALE2[iRadius] = jets[iRadius][1].perp2();
      
      int scaleTest = 0;
      for (int iRadius = 0; iRadius < nRadius ; iRadius++) 
	if( SCALE2[iRadius] > 0. ) scaleTest++;
      
      if ( scaleTest <= 0 )
	{
	  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
	    {
	      std::cout <<"\n\n "<<  __PRETTY_FUNCTION__ 
			<< " SCALE is wrong " << SCALE2[iRadius] <<" RADIUS = "<<jetSizes[iRadius]<<" \n Event : \n"<< p
			<<std::endl;
	      
	      std::cout <<" Njets check : upper()"<<jets[iRadius].upper() <<" MinNumber = "<< UserHHC::numberOfJets<< std::endl;
	      std::cout <<"List of jets for radius "<<jetSizes[iRadius]<< std::endl;
	      for(int jet = 1; jet <= (jets[iRadius]).upper(); jet++)
		std::cout <<" jet "<<jet<<" "<<jets[iRadius][jet]<<" pt = "<<jets[iRadius][jet].perp()<<std::endl;
	      std::cout <<"**************************************"<<std::endl;
	      
	      goodEvent = false;
	      return;
	    }
	}

      if (debug)
	{
	  std::cout<<"\n"<<std::endl;
	  for (int iRadius = 0; iRadius < nRadius ; iRadius++) 
	    {
	      std::cout <<" radius = "<<jetSizes[iRadius]<<endl;
	      for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
		{
		  std::cout <<" grid = "<<iGrid<<endl;
		  for(int jet = 1; jet <= pj[iRadius][iGrid].upper(); jet++)
		    std::cout <<"jet # "<< jet<<" "<<pj[iRadius][iGrid][jet]<<" pt = "<< pj[iRadius][iGrid][jet].perp()<<std::endl;

		}
	    }
	} 
    } // isgood

  return;
}

//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//**********************************************************************
//
//   User analysis
//
void UserHHC::userfunc(const event_hhc& p, const amplitude_hhc& amp)
{ 
  if (debug) cout<<" evtNb = "<<eventNb<<"\t contrib = "<<amp.contrib()<<endl;

  if (debug) std::cout <<"*************************************************************"<< std::endl;
  if (debug) std::cout <<"*************************************************************"<< std::endl;
  if (debug) std::cout <<"*************************************************************"<< std::endl;
  
  clusterJets(p);
  
  if (debug) std::cout <<"*************************************************************"<< std::endl;
  if (debug) std::cout <<"*************************************************************"<< std::endl;
  if (debug) std::cout <<"*************************************************************"<< std::endl;
  
  if ( !( this->isGoodEvent() ) ) return;
  
  //
  //  grid works...
  //
  for (int iRadius = 0; iRadius < nRadius ; iRadius++) 
    for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
      {
	if ( eventSelected[iRadius][iGrid])
	  { 
	    for(int jet = 1; jet <= (pj[iRadius][iGrid]).upper(); jet++) 
	      { //loop over jets
		double ptJet  = ((pj[iRadius][iGrid])[jet]).perp();
		mygrid[iRadius][iGrid]->fill(x1, x2, ptJet * ptJet, ptJet, amp, pdf);
///////		mygrid[iRadius][iGrid]->fill(x1, x2, SCALE2[iRadius], ptJet, amp, pdf);
	      }
	  }             // if (eventSelected)
      }                 // for (;;iGrid++)
  return; 
}     

// end of UserHHC::userfunc
void UserHHC::operations_at_the_end_of_event()
{
  if (debug) std::cout << __PRETTY_FUNCTION__ << std::endl;
  
  for (int iRadius = 0; iRadius < nRadius ; iRadius++)
    for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
      mygrid[iRadius][iGrid]->endOfEvent();
  
  
  eventNb++;
  // persistency
  // do some things every saveAfterEvents
  if( ( eventNb % saveAfterEvents == 0 )  || ( eventNb == maxNumberEvents ) ) 
    {
      cout<<"saving grid after "<< eventNb <<" runs"<<endl;
      
      for (int iRadius = 0; iRadius < nRadius ; iRadius++)
	for (int iGrid = oGrids; iGrid < nGrids; iGrid++)
	  mygrid[iRadius][iGrid]->writeGrid( eventNb );
      
      
    }
  
  // remeber to run NLOjet's official end of event operations
  user1h_hhc::operations_at_the_end_of_event();

  if ( eventNb == maxNumberEvents ) 
    {
      this->deleteObjects();
      exit(0);
    }
  return;
}
