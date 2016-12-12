/*
 *  appl_qcd.h
 *  QCD interfaces for applcomb
 * *  nph 09/14
 */

#pragma once

#include "APFEL/APFEL.h"
#include "APFEL/APFELdev.h"
 
#include "NNPDF/common.h"

namespace appl{
class grid;
}

namespace NNPDF{
class FKGenerator;
class FKHeader;
}

namespace QCD
{
  // ***************************************************************
  // ********************* BASIS ROTATION **************************

   // APFELCOMB basis is the same as the APFEL basis

     // Evolution basis
  //  γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35

     // Flavour basis
  //  γ, TB, BB, CB, UB, DB, G, D, U, C, B, T
  
  // Rotation matrix from LHA to EVLN
  static const double REVLN2LHA[14][14] = {
    { 1,  0,       0,        0,       0,        0,          0,         0,         0,       0,       0,         0,         0, 0       },
    { 0,  1/12.0,  0,  -1/12.0,       0,        0,          0,         0,    1/12.0,       0,       0,         0,         0, -1/12.0 },
    { 0,  1/12.0,  0,  -1/12.0,       0,        0,          0,  12/120.0,  -2/120.0,       0,       0,         0, -12/120.0, 2/120.0 },
    { 0,  1/12.0,  0,  -1/12.0,       0,        0,   15/120.0,  -3/120.0,  -2/120.0,       0,       0, -15/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,  -1/12.0,       0,   2/12.0,   -5/120.0,  -3/120.0,  -2/120.0,       0, -2/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,  -1/12.0, -3/12.0,  -1/12.0,   -5/120.0,  -3/120.0,  -2/120.0,  3/12.0,  1/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,  -1/12.0,  3/12.0,  -1/12.0,   -5/120.0,  -3/120.0,  -2/120.0, -3/12.0,  1/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  0,       1,        0,       0,        0,          0,         0,         0,       0,       0,         0,         0, 0       },
    { 0,  1/12.0,  0,   1/12.0, -3/12.0,   1/12.0,    5/120.0,   3/120.0,   2/120.0, -3/12.0,  1/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,   1/12.0,  3/12.0,   1/12.0,    5/120.0,   3/120.0,   2/120.0,  3/12.0,  1/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,   1/12.0,       0,  -2/12.0,    5/120.0,   3/120.0,   2/120.0,       0, -2/12.0,   5/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,   1/12.0,       0,        0,  -15/120.0,   3/120.0,   2/120.0,       0,       0, -15/120.0,   3/120.0, 2/120.0 },
    { 0,  1/12.0,  0,   1/12.0,       0,        0,          0, -12/120.0,   2/120.0,       0,       0,         0, -12/120.0, 2/120.0 },
    { 0,  1/12.0,  0,   1/12.0,       0,        0,          0,         0,   -1/12.0,       0,       0,         0,         0, -1/12.0 }
  };
  
  // Rotation matrix from EVLN to LHA
  static const double RLHA2EVLN[14][14] = {
    { 1, 0,   0,  0,  0,  0,  0, 0,  0, 0,  0,  0,  0, 0  },
    { 0, 1,   1,  1,  1,  1,  1, 0,  1, 1,  1,  1,  1, 1  },
    { 0, 0,   0,  0,  0,  0,  0, 1,  0, 0,  0,  0,  0, 0  },
    { 0, -1, -1, -1, -1, -1, -1, 0,  1, 1,  1,  1,  1, 1  },
    { 0, 0,   0,  0,  0, -1,  1, 0, -1, 1,  0,  0,  0, 0  },
    { 0, 0,   0,  0,  2, -1, -1, 0,  1, 1, -2,  0,  0, 0  },
    { 0, 0,   0,  3, -1, -1, -1, 0,  1, 1,  1, -3,  0, 0  },
    { 0, 0,   4, -1, -1, -1, -1, 0,  1, 1,  1,  1, -4, 0  },
    { 0, 5,  -1, -1, -1, -1, -1, 0,  1, 1,  1,  1,  1, -5 },
    { 0, 0,   0,  0,  0,  1, -1, 0, -1, 1,  0,  0,  0, 0  },
    { 0, 0,   0,  0, -2,  1,  1, 0,  1, 1, -2,  0,  0, 0  },
    { 0, 0,   0, -3,  1,  1,  1, 0,  1, 1,  1, -3,  0, 0  },
    { 0, 0,  -4,  1,  1,  1,  1, 0,  1, 1,  1,  1, -4, 0  },
    { 0, -5,  1,  1,  1,  1,  1, 0,  1, 1,  1,  1,  1, -5 }
  };

  
  // Basis transforms
  void LHA2EVLN ( const double* LHA, NNPDF::real* EVLN );
  void EVLN2LHA ( const double* EVLN, NNPDF::real* LHA );
  
  // ***************************************************************
  // ********************** Parameters *****************************

  // QCD data struct
  class qcd_param
  {  
  public:
    size_t thID;      //!< Theory ID
    std::map<std::string, std::string> thMap; //!< Map of theory parameters obtained from database

    size_t evol_pto;  //!< Perturbative order of the evolution
    size_t pto;       //!< Perturbative order of the produced grid

    double Q0;        //!< Initial Q^2
    double xiF;       //!< mu_F/Q_F
    double xiR;       //!< mu_R/Q_R
  };

  // Parse theory input
  void parse_input( int innum, qcd_param&  param);
  void set_params(qcd_param const& par, NNPDF::FKHeader& FK);

  // ***************************************************************
  
  // Initialise the APFEL interface
  void initQCD(qcd_param const& par, const double& Q2max);
  void initEvolgrid(const int& nx, const double& xmin);
  
  void initNNPDF30Grid(); // Initialise NNPDF30 style x-grid for test purposes
  void initTruthGrid(const double& xmin); // Initialise very large evolution grid for test purposes

  // APFEL PDF access functions
  void evolpdf(const double& x, const double& Q, double* pdf);
  void evolpdf_applgrid(const double& x, const double& Q, double* pdf); // Version for applgrid (no photon)
  void evolpdf_applgrid_pbar(const double& x, const double& Q, double* pdf); // Version for applgrid (no photon) with antiproton pdfs
  double alphas(const double& Q);
  double beta0();
  
  // APFEL FK functions
  void avals(const int& xi, const double& xo, const int& fi, const double& Q, double*);
  void avals_pbar(const int& xi, const double& xo, const int& fi, const double& Q, double*);
  void bvals(int const& xi, double const& xo, int const& fi, double* b);

  double diskernel(std::string const& obs, double const& x, double const& Q, double const& y, int const& i, int const& beta);
  double disobs(std::string const& obs, double const& x, double const& Q, double const& y);
  
  // mode setters
  void setDISmode(bool const& mode);
  void setFTDYmode(bool const& mode);
  void setSIAmode(bool const& mode);
}
