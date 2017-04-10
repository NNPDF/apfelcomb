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

  // Flavour basis (-7 to 6)
  //  γ, TB, BB, CB, SB, UB, DB, G, D, U, S, C, B, T

  // APPLgrid basis (0 to 12)
  //  TB, BB, CB, SB, UB, DB, G, D, U, S, C, B, T

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
    double Q0;        //!< Initial Q^2
    double xiF;       //!< mu_F/Q_F
    double xiR;       //!< mu_R/Q_R
  };

  // Parse theory input
  void parse_input( int innum, qcd_param&  param);
  void set_params(qcd_param const& par, NNPDF::FKHeader& FK);
  vector<size_t> active_flavours(qcd_param const& param);

  // ***************************************************************
  
  // Initialise the APFEL interface
  void initQCD(qcd_param& par, const bool& positivity, const double& Q2max);
  void initEvolgrid(const int& nx, const double& xmin);
  void initPDF(std::string const& setname, int const& i);

  // APFEL PDF access functions
  void evolpdf(const double& x, const double& Q, double* pdf);
  double alphas(const double& Q);
  double beta0();
  
  // PDF Operators
  void EvolutionOperator(const bool& ppbar, const int& xi, const double& xo, const int& fi, const double& Q, double*);
  void DerivativeOperator(const bool& ppbar, const int& beta, const double& alpha, const int& j, const double& Q, double* a);

  double diskernel(std::string const& obs, double const& x, double const& Q, double const& y, int const& i, int const& beta);
  double disobs(std::string const& obs, double const& x, double const& Q, double const& y);
  
  // mode setters
  void setFTDYmode(bool const& mode);
  void setSIAmode(bool const& mode);
}
