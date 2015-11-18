/*
 *  appl_qcd.h
 *  QCD interfaces for applcomb
 * *  nph 09/14
 */

#pragma once

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
    size_t evol_pto;  //!< Perturbative order of the evolution
    size_t pto;       //!< Perturbative order of the produced grid
    
    string FNS;       //!< Flavour number scheme to use in the evolution
    bool   damp;      //!< FONLL damping switch
    bool   IC;        //!< Intrinsic charm switch

    string MODEV;     //!< Mode of DGLAP solution
    double xiR;       //!< Xi_R renormalisation scale factor
    double xiF;       //!< Xi_F factorisation scale factor

    size_t nff;        //!< Number of flavours in FFN
    size_t nf_as;   //!< Number of flavours in alphas running
    size_t nf_pdf;  //!< Number of flavours in PDF running

    double Q0;        //!< Initial Q^2
    double alphas;    //!< Value of alpha_s
    double QREF;      //!< Reference QCD scale (Typically M_Z)

    bool QED;         //!< Flag to enable QED evolution
    double alpha_qed; //!< Value of alpha_qed
    double QEDREF;    //!< Reference QED scale

    bool SxRes;       //!< Small-x resummation switch
    string SxOrd;     //!< Small-x resummation order

    string HQMASS;    //!< Heavy quark mass (POLE/DGLAP)
    
    double mc;       //!< Mass of charm
    double mb;       //!< Mass of bottom
    double mt;       //!< Mass of top
    double Qmc;       //!< Charm mass reference scale
    double Qmb;       //!< Bottom mass reference scale
    double Qmt;       //!< Top mass reference scale

    double CKM[3][3]; //!< CKM matrix

    double mz;  //!< Z mass
    double mw;  //!< W mass

    double gf; //!< G_Fermi
    double sin2tw; //!< Sin^2 theta_w

    bool TMC;   //!< Target mass corrections  
    double Mt;  //!< Target mass

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
  double alphas(const double& Q);
  
  // APFEL FK functions
  void avals(const int& xi, const double& xo, const int& fi, const double& Q, double*);
  double diskernel(std::string const& obs, double const& x, double const& Q, double const& y, int const& i, int const& beta);
  double disobs(std::string const& obs, double const& x, double const& Q, double const& y);
  
  // Initialised x-grid values
  size_t  getNXGrid();
  double  getXVal(const int idx);

  // mode setters
  void setDISmode(bool const& mode);
  void setFTDYmode(bool const& mode);
  void setSIAmode(bool const& mode);

  // Proton mass fetcher (for xi)
  double getProtonMass();

}
