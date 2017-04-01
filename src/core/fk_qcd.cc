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
#include <string>

using NNPDF::FKHeader;

namespace QCD
{

  // Set static DIS mode
  static bool FTDY_mode = false;
  static bool SIA_mode = false;

  static std::string FKObs; //!< Cached FK observable
  
  static double Q0 = 0; // Initial scale  Q0
  static double QM = 0; // Maximum scale  QM
  static double QC = 0; // Cached scale   QC

  // Set modes
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


  // Returns a list of active flavours in the APFEL evolution basis, for a given initial scale
  //  γ, Σ, g, V, V3, V8, V15, V24, V35, T3, T8, T15, T24, T35
  vector<size_t> active_flavours(qcd_param const& param)
  {
    vector<size_t> afl = {1,2,3,4,5,9,10};
    if (stoi((*param.thMap.find("QED")).second) == 1) afl.push_back(0); // γ

    if (param.Q0 > APFEL::HeavyQuarkThreshold(4)) afl.push_back(6); // V15
    if (param.Q0 > APFEL::HeavyQuarkThreshold(5)) afl.push_back(7); // V24
    if (param.Q0 > APFEL::HeavyQuarkThreshold(6)) afl.push_back(8); // V35

    if (param.Q0 > APFEL::HeavyQuarkThreshold(4)) afl.push_back(11); // V15
    if (param.Q0 > APFEL::HeavyQuarkThreshold(5)) afl.push_back(12); // V24
    if (param.Q0 > APFEL::HeavyQuarkThreshold(6)) afl.push_back(13); // V35   
    return afl;
  }

  // *********************** FKGenerator Parameters *************************

  void parse_input(int innum, qcd_param& param)
  {
    // Theory index
    param.thID = innum;

    // Load theory choices DB
    NNPDF::IndexDB theoryDB(dataPath()+"theory.db", "theoryIndex");
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

    // Set scale variation parameters
    param.xiF = atof(param.thMap["XIF"].c_str());
    param.xiR = atof(param.thMap["XIR"].c_str());

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

    const time_point st = std::chrono::system_clock::now();
    const std::time_t start_time = std::chrono::system_clock::to_time_t(st);
    FK.AddTag(FKHeader::VERSIONS, "GenTime", ctime(&start_time));
  }
  
  // *********************** EVOLUTON FUNCTIONS *****************************

  // Initialise QCD according to parameters
  void initQCD(qcd_param& par, const bool& positivity, const double& Q2max)
  {
    // Disable TMCs in positivity observables
    if (positivity)
    {
      par.thMap["TMC"] = '0'; 
      par.thMap["XIF"] = '1'; par.xiF = 1; 
      par.thMap["XIR"] = '1'; par.xiR = 1;
    }

    APFEL::SetParam(par.thMap);

    // Init Q0
    QCD::Q0 = par.Q0;
    QCD::QM = std::max(par.xiF, par.xiR)*std::sqrt(Q2max);
    APFEL::SetQLimits( Q0, QM );

    APFEL::SetEpsilonTruncation(1E-1);
    APFEL::SetFastEvolution(false);
    APFEL::LockGrids(false);
    APFEL::EnableEvolutionOperator(true); 

    if (SIA_mode) 
    {
      APFEL::SetPDFSet("kretzer");
      APFEL::SetTimeLikeEvolution(true);
    }

    // Start APFEL
    APFEL::InitializeAPFEL_DIS();

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
    std::cout << " Initialising  "<< nx << " points starting from x = " << xmin <<std::endl; 


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
    
    // Set scale limits
    APFEL::SetQLimits( Q0, QM );

    // Set evolution operator parameters
    APFEL::SetFastEvolution(false);
    APFEL::LockGrids(false);
    APFEL::EnableEvolutionOperator(true); 
    APFEL::SetNumberOfGrids(1);           
    APFEL::SetExternalGrid(1,nx,3,xg);  

    // Initialise
    APFEL::EnableWelcomeMessage(false);

    // Start APFEL
    APFEL::InitializeAPFEL_DIS();

    // Needed to join the grids
    APFEL::EvolveAPFEL(Q0,Q0);
    
    delete[] xg;
    return;
  }

  void initPDF(std::string const& PDFSet, int const& i) 
  {
    APFEL::SetPDFSet(PDFSet);
    APFEL::SetReplica(i);
  };


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
    
    // Fuzzy x
    const double fx = fabs(x-APFEL::xGrid(0))/x < 1E-6 ? APFEL::xGrid(0):x;    
    for (int i=-7; i<7; i++)
      pdf[i+7]= (i==-7 ? 0:APFEL::xPDF(i,fx));
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
  
  // PDF evolution+rotation operator
  // This functions computes the operators A such that
  // PDF(fo,xo,Q) = \sum_xi,fi A(fo,fi,xo,xi)*N(fi,xi,Q_0)
  // where N is the initial scale APFEL evolution basis PDF and PDF is in the APPLgrid physical basis.
  // inputs:
  // ppbar - bool to switch flavours for the case of incoming antiproton beam
  // xi    - index in APFEL x-grid corresponding to the input PDF
  // xo    - value of x for the evolved PDF
  // fi    - input flavour (in APFEL evolution basis)
  // Q     - scale of the final evolved PDF
  // a     - array in output flavour fo (APPLgrid flavour basis) of evolution operators
  void EvolutionOperator(const bool& ppbar, const int& xi, const double& xo, const int& fi, const double& Q, double* a)
  {
    updateEvol(Q);
    for(int i=0; i<13; i++)
      a[ppbar ? 12-i:i] = APFEL::ExternalEvolutionOperator(std::string("Ev2Ph"),i-6,fi,xo,xi);
    return;
  }

  // LO Evolved PDF derivative operator
  // This function computes the operators D such that
  // dPDF(fo,xo,Q)/dQ = \sum_xi,fi D(fo,fi,xo,xi,Q)*N(fi,xi,Q_0)
  // where N is the initial scale APFEL evolution basis PDF and the derivatives are evaluated in the APPLgrid physical basis.
  // inputs:
  // ppbar - bool to switch flavours for the case of incoming antiproton beam
  // xi    - index in APFEL x-grid corresponding to the input PDF
  // xo    - value of x for which the derivative is evaluated
  // fi    - input flavour (in APFEL evolution basis)
  // Q     - scale at which the derivative is evaluated
  // da    - array in output flavour fo (APPLgrid flavour basis) of derivative operators
  void DerivativeOperator(const bool& ppbar, const int& xi, const double& xo, const int& fi, const double& Q, double* da)
  {
    const int pt = 0;
    const int nf = 5;
    updateEvol(Q);
    
    for(int fo=0; fo<13; fo++)
      da[fo] = 0;

    for(int ixd=0; ixd<APFEL::nIntervals(); ixd++)
    {
      APFEL::ComputeExternalSplittingFunctions("Ev2Ph", pt, nf, xo, ixd);
      const double xd = APFEL::xGrid(ixd);
        for(int fd=0; fd<14; fd++)
        {
          const double A = APFEL::ExternalEvolutionOperator(std::string("Ev2Ev"), fd, fi, xd, xi);
          for(int fo=0; fo<13; fo++)
            da[ppbar ? 12-fo:fo] += 0.5*APFEL::ExternalSplittingFunctions(fo-6,fd)*A;  // 0.5 due to APFEL expansion parameter
        }
    }
  }


  // *********************************************** DIS functions ***********************************************

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

