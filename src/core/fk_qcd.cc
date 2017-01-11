/*
 *  appl_qcd.cc
 *  QCD interfaces for applcomb
 * *  nph 09/14
 */

#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_xgrid.h"

#include "APFEL/APFEL.h"
#include "APFEL/APFELdev.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"

#include "NNPDF/common.h"
#include "NNPDF/fkgenerator.h"
#include "NNPDF/nnpdfdb.h"

#include <limits>

using NNPDF::FKHeader;

namespace QCD
{

  // Set static DIS mode
  static bool DIS_mode = false;  
  static bool FTDY_mode = false;
  static bool SIA_mode = false;

  static std::string FKObs; //!< Cached FK observable
  
  static double Q0 = 0; // Initial scale  Q0
  static double QM = 0; // Maximum scale  QM
  static double QC = 0; // Cached scale   QC

  // Set modes
  void setDISmode(bool const& mode) {DIS_mode=mode;};
  void setFTDYmode(bool const& mode) {FTDY_mode=mode;};
  void setSIAmode(bool const& mode) {SIA_mode=mode;};

  // *********************** BASIS ROTATION *****************************
  
  // Rotate evolution basis PDFs to flavour basis
  void EVLN2LHA(const double* EVLN, NNPDF::real* LHA)
  {
    for (int i=0; i<14; i++)
    {
      LHA[i] = 0;
      for (int j=0; j<14; j++)
        LHA[i] += REVLN2LHA[i][j]*EVLN[j];
    }
  }
  
  // Rotate flavour basis PDFs to evolution basis
  void LHA2EVLN(const double* LHA, NNPDF::real* EVLN)
  {
    for (int i=0; i<14; i++)
    {
      EVLN[i] = 0;
      for (int j=0; j<14; j++)
        EVLN[i] += RLHA2EVLN[i][j]*LHA[j];
    }
  }

  // *********************** FKGenerator Parameters *************************

  void parse_input(int innum, qcd_param& param)
  {
    // First of all, ensure that libnnpdf is compiled in safe mode
    #ifdef SSE_CONV
      cout << "Error: libnnpdf must be compliled in safe mode for use with apfelcomb"<<endl;
      exit(-1);
    #endif

    // Theory index
    param.thID = innum;

    // Load theory choices DB
    NNPDF::IndexDB theoryDB(resultsPath()+"theory.db", "theoryIndex");
    const int entries =theoryDB.GetNEntries();
    if (innum < 0 || innum > entries)
    {
      cerr << "Error: theory ID ("<<innum<<") must be between 0 and "<<entries-1<<endl;
      exit(-1);
    }

    // Fetch settings map
    theoryDB.ExtractMap(param.thID, APFEL::kValues, param.thMap);

    // Set perturbative order and initial scale
    param.Q0 = atof(param.thMap["Q0"].c_str());
    param.evol_pto = atoi(param.thMap["PTO"].c_str());
    if (param.evol_pto == 2 and DIS_mode == false and FTDY_mode == false)
    {
      cout << endl<<endl;
      cout << " ****************** WARNING ****************** "<<endl;
      cout << " *    APPLGRID CONVOLUTION WITH N2LO NOT     * "<<endl;
      cout << " *    ENABLED - CONVOLUTION FIXED TO N1LO    * "<<endl;
      cout << " ********************************************* "<<endl;
      cout << endl<<endl;
      
      // Set internal PTO to N1LO
      param.pto = 2;
    } else
    {
      param.pto = param.evol_pto + 1;
    }

    // Set scale variation parameters
    param.xiF = atof(param.thMap["XIF"].c_str());
    param.xiR = atof(param.thMap["XIR"].c_str());

    cout << "                FastKernel Grid Combination                 "<<endl;
    cout << "    - TheoryID: "<<param.thID << endl;
    cout << "    - Conv PTOrd: N"<<param.pto -1<<"LO"<<endl;
    cout << "    - xi_F: " << param.xiF<<endl;
    cout << "    - xi_R: " << param.xiR<<endl;

    return;
  }


  void set_params(qcd_param const& par, NNPDF::FKHeader& FK)
  {
    // Set FK table description
    FK.AddTag(FKHeader::VERSIONS, "APFEL", APFEL::GetVersion());
    FK.AddTag(FKHeader::VERSIONS, "libnnpdf", NNPDF::getVersion());

    FK.AddTag(FKHeader::GRIDINFO, "NX", APFEL::nIntervals());

    // x-grid header
    stringstream xGheader;
    for (int i=0; i<APFEL::nIntervals(); i++)
      xGheader << std::setprecision(16) << std::scientific << APFEL::xGrid(i) <<std::endl;

    FK.AddTag(FKHeader::BLOB, "xGrid", xGheader.str());

    // Theory configuration
    std::map<std::string, std::string>::const_iterator imap;
    for (imap = par.thMap.begin(); imap != par.thMap.end(); imap++)
      FK.AddTag(FKHeader::THEORYINFO, imap->first, imap->second);
  }
  
  // *********************** EVOLUTON FUNCTIONS *****************************

  // Initialise QCD according to parameters
  void initQCD(qcd_param const& par, const double& Q2max)
  {
    APFEL::SetParam(par.thMap);

    // Init Q0
    QCD::Q0 = par.Q0;

    // Truncated Epsilon
    APFEL::SetEpsilonTruncation(1E-1);

    // Set maximum scale
    QM = std::max(par.xiF, par.xiR)*std::sqrt(Q2max);
    if ( fabs(par.xiF - 1.0) > 1E-5)
      QM = std::max(15000., QM); // HOPPET max scale in APPLgrid
    APFEL::SetQLimits( Q0, QM );

    if (SIA_mode) 
    {
      APFEL::SetPDFSet("kretzer");
      APFEL::SetTimeLikeEvolution(true);
    }

    // Start APFEL
    if (DIS_mode)
      APFEL::InitializeAPFEL_DIS();
    else
      APFEL::InitializeAPFEL();

    return;
  }

 // Initialise APFEL for evolution factors
  void initTruthGrid(const double& xmin)
  {
    // Set truth grid
    APFEL::SetNumberOfGrids(2);
    APFEL::SetGridParameters(1,60,5,xmin);
    APFEL::SetGridParameters(2,60,5,1e-1);
    
    // Initialise
    APFEL::LockGrids(true);
    APFEL::EnableWelcomeMessage(false);

    // Start APFEL
    if (DIS_mode)
      APFEL::InitializeAPFEL_DIS();
    else
      APFEL::InitializeAPFEL();
  
    return;
  }
  
  // Initialise APFEL for evolution factors
  // Note here nx is the number of x-points to be output to the Fk table
  // Therefore it doesn't include x=1. This is added manually in this function
  void initEvolgrid(int const& nx, double const& xmin)
  {
    // Reset cache
    QC = 0;
    double* xg = new double[nx+1];

    // Special requirements for FTDY in APFEL
    if (FTDY_mode)
    {
      // Requires two x-grid points below x-min
      xg[0] = 0.9*xmin;
      xg[1] = 0.95*xmin;

      const double ymin = XGrid::appl_fy(xmin);
      const double ymax = XGrid::appl_fy(1.0);

      // Populate grid
      for (int i=2; i<=nx; i++)
        xg[i] = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx))*(i-2));

    }
    else
    { // Normal mode
      const double ymin = XGrid::appl_fy(0.99*xmin); 
      const double ymax = XGrid::appl_fy(1.0);
      
      // Populate grid
      for (int i=0; i<=nx; i++)
        xg[i] = XGrid::appl_fx(ymin + ((ymax-ymin)/((double) nx))*i);
    }
    
    // Just in case of numerical trouble
    xg[nx] = 1;
    
  
    // Set evolution operator parameters
    APFEL::SetFastEvolution(false);
    APFEL::EnableEvolutionOperator(true); // Enable the computation of the Evolution Operator
    APFEL::SetNumberOfGrids(1);           // The evolution will be done on a number of subgrids defined here
    APFEL::SetExternalGrid(1,nx,5,xg);    // Set the grid as external (np: number of intervals, 3: interpolation degree, xg: defined above)
    
    // Initialise
    APFEL::EnableWelcomeMessage(false);

    // Start APFEL
    if (DIS_mode)
      APFEL::InitializeAPFEL_DIS();
    else
      APFEL::InitializeAPFEL();

    // Needed to join the grids
    APFEL::EvolveAPFEL(Q0,Q0);
    
    delete[] xg;
    return;
  }

  // Recalculate evolution if not cached
  static void updateEvol(const double& Q){
    if (Q != QC)
    {
      if (fabs(Q-Q0) < 1E-7) // Fuzzy comp
        APFEL::EvolveAPFEL(Q0,Q0);
      else
        APFEL::EvolveAPFEL(Q0,Q);
      QC = Q;
    }
  }

  // APFEL PDF return for APPLgrid (no photon!)
  void evolpdf_applgrid(const double& x, const double& Q, double* pdf)
  {
    // A nice trick of APPLgrid is to request PDF x-values smaller than
    // are actually used
    if (x<APFEL::xGrid(0))
    {
      for (int i=-6; i<7; i++)
        pdf[i+6]=0;
      return;
    }
    updateEvol(Q);
    for (int i=-6; i<7; i++)
      pdf[i+6]=APFEL::xPDF(i,x);
  }

  // APFEL PDF return for APPLgrid (no photon!) - antiproton version
  void evolpdf_applgrid_pbar(const double& x, const double& Q, double* pdf)
  {
    // A nice trick of APPLgrid is to request PDF x-values smaller than
    // are actually used
    if (x<APFEL::xGrid(0))
    {
      for (int i=-6; i<7; i++)
        pdf[i+6]=0;
      return;
    }
    updateEvol(Q);
    for (int i=-6; i<7; i++)
      pdf[-i+6]=APFEL::xPDF(i,x);
  }

    // APFEL PDF return
  void evolpdf(const double& x, const double& Q, double* pdf)
  {
    // Recalculate if not cached
    if (Q != QC)
    {
      if (fabs(Q-Q0) < 1E-7) // Fuzzy comp
        APFEL::EvolveAPFEL(Q0,Q0);
      else
        APFEL::EvolveAPFEL(Q0,Q);
      QC = Q;
    }
    
    for (int i=-7; i<7; i++)
      pdf[i+7]= (i==-7 ? 0:APFEL::xPDF(i,x));
  }

  // APFEL strong coupling
  double alphas(const double& Q)
  {
    return APFEL::AlphaQCD(Q);
  }

  double beta0()
  {
    const double twopi = 2.*M_PI;
    const double nc = 3;
    const double nf = 5;  // This is what APPLgrid does
    const double beta0=(11.*nc-2.*nf)/(6.*twopi);
    return beta0;
  }
  
  // Evaluate A values
  void avals(const int& xi, const double& xo, const int& fi, const double& Q, double* a)
  {
    updateEvol(Q);
    for(int i=0; i<13; i++)
      a[i] = APFEL::ExternalEvolutionOperator(std::string("Ev2Ph"),i-6,fi,xo,xi);
    return;
  }
  // Evaluate A values -ppbar
  void avals_pbar(const int& xi, const double& xo, const int& fi, const double& Q, double* a)
  {
    updateEvol(Q);
    for(int i=0; i<13; i++)
      a[-1*(i-6) + 6] = APFEL::ExternalEvolutionOperator(std::string("Ev2Ph"),i-6,fi,xo,xi);
    return;
  }

  double bvals(int const& xo, int const& xi, int const& fo, int const& fi)
  {
    const int pt = 0;
    const int nf = 5;
    return APFEL::ExternalSplittingFunctions(std::string("Ph2Ph"),pt,nf,fo-6,fi-6,APFEL::xGrid(xo),xi);
  }

  double diskernel(std::string const& obs, double const& x, double const& Q, double const& y, int const& i, int const& beta)
  {
    // Out of range
    if (Q < Q0) return 0;

    // Recalculate if not cached
    if (Q != QC || FKObs.compare(obs) != 0)
    {
      APFEL::SetFKObservable(obs);
      APFEL::ComputeStructureFunctionsAPFEL(Q0,Q);

      // Set cache values
      QC = Q;
      FKObs = obs;
    }

    // Return the DIS kernel
    return APFEL::FKSimulator(x,Q,y,i,beta);
  }

  double disobs(std::string const& obs, double const& x, double const& Q, double const& y)
  {
    // Out of range
    if (Q < Q0) return 0;

    // Recalculate if not cached
    if (Q != QC || FKObs.compare(obs) != 0)
    {
      APFEL::SetFKObservable(obs);
      APFEL::ComputeStructureFunctionsAPFEL(Q0,Q);

      // Set cache values
      QC = Q;
      FKObs = obs;
    }

    return APFEL::FKObservables(x,Q,y);
  }


}

