
#ifndef __SCALES_H
#define __SCALES_H

const int nSubProcesses = 7;

const int Nscales=7;
const int nScales=Nscales;
static const double mur[Nscales] = { 1.0, 0.5, 2.0, 1.0, 1.0, 0.5, 2.0 };
static const double muf[Nscales] = { 1.0, 1.0, 1.0, 0.5, 2.0, 0.5, 2.0 };


const double pb_fac             = 3.89379656e8 ;    // conversion GeV^2 -> pb  
const double pb_fac_Unity       = 1.0;              // unit factor

const int nRadius = 1;
const float jetSizes[nRadius] = { 0.6 };

//
//
//

const int oGrids = 0;
const int nGrids = 6;
const int netaBins = nGrids;

static const float etaBinsLow[ netaBins ] =
  {
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5  
  };

static const float etaBinsHigh[ netaBins ] =
  {
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0 
  };


static const int nPtBins = 20;                       //no of bins
static const double ptBins[nPtBins + 1] =            //bin-information, lower edges...
  { 
    10., 20., 30., 45., 60., 80., 110., 160., 210., 
    260., 310., 400., 500., 600., 800., 1000., 1200., 1500., 1800., 2500., 3500.
  };


static const int nPtBins0 = 31;
static const double ptBins0[ nPtBins0 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 346, 376, 408, 442, 478, 516, 556, 598, 642, 688, 736, 786, 838, 894, 952, 1012, 1076, 1162, 1310, 1530, 1992 };

static const int nPtBins1 = 29;
static const double ptBins1[ nPtBins1 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 346, 376, 408, 442, 478, 516, 556, 598, 642, 688, 736, 786, 838, 894, 952, 1012, 1162, 1310, 1992 };

static const int nPtBins2 = 26;
static const double ptBins2[ nPtBins2 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 346, 376, 408, 442, 478, 516, 556, 598, 642, 688, 736, 786, 838, 894, 1012, 1992 };

static const int nPtBins3 = 23;
static const double ptBins3[ nPtBins3 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 346, 376, 408, 442, 478, 516, 556, 598, 642, 688, 736, 894, 1992 };

static const int nPtBins4 = 19;
static const double ptBins4[ nPtBins4 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 346, 376, 408, 442, 478, 516, 556, 642, 894 };

static const int nPtBins5 = 12;
static const double ptBins5[ nPtBins5 + 1 ] = { 100, 116, 134, 152, 172, 194, 216, 240, 264, 290, 318, 376, 478 };

const double pT_min     = 0.95 * ptBins0[ 0 ] ;
const double etaAccept  = 3.0;


#endif // __SCALES_H

/*
const int nGrids = 9;
const int netaBins = nGrids;

static const float etaBinsLow[ netaBins ] =
  {
    0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0
  };

static const float etaBinsHigh[ netaBins ] =
  {
    0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.4
  };
*/
