#ifndef NLOGRID__HH
#define NLOGRID_HH 1

#include <string>
#include "appl_grid/appl_grid.h"
#include "TH1D.h"
#include <bits/hhc-process.h>
#include <algorithm>

#include "pdf-genlha.h"

#include "scales.h"

class nlogrid
{
 public:
  nlogrid(const std::string& inputName = "test",  const int &binNumber = 0);
  ~nlogrid();

 public:

  inline int getMode(){return mode;};

  TH1D* getGridReference(){return gridObject->getReference();};

  appl::grid* getGridObject() { return gridObject; }


  void writeGrid(long long&);
  void endOfEvent();
  
  void fill(
	    const double &x1, 
	    const double &x2, 
	    const double &SCALE2, 
	    const double &obs,
	    const nlo::amplitude_hhc& amp,
	    nlo::pdf_and_coupling_hhc* pdf
	    );
  
#ifdef REN_REFERENCE

  void bookReferenceHistograms(const int&, const double *);
  void deleteReferenceHistograms();
  void writeReferenceHistograms();
#endif
  void fillReferenceHistograms( double &obs, 
				double &xs,
				double* xsScales
			       );
  
  

 private:

  appl::grid* gridObject;
  std::string fullFileName;      // name of the file to store grid
  std::string observableName;    // name of the file to store grid
  int mode;                      // mode of grid (0 = new ; 1 = optimisation)
  long int numOfEvents, numOfBadEvents;
  bool isEventBad, debug;
  
  std::string refDirName;

 private:
 struct appl_event
  {
    double x1;
    double x2;
    double Q2;
    double obs;
    double xSection;
    double xSectionScale[ nScales ];
    double weights[ nSubProcesses ];
    int processOrder;
  };

 std::vector<appl_event>  currentEvent;
 
 appl_pdf             * applpdf;                           // always returns 1.
 //  conversion form weight_hhc to double
 weight_conversion<weight_hhc, double> conversion;
 
#ifdef REN_REFERENCE
 
  // reference histograms to test renormalisation and factorisation scale dependence
  TH1D* referenceScale [Nscales];

#endif

};



#endif
