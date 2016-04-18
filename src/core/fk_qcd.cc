/*
 *  appl_qcd.cc
 *  QCD interfaces for applcomb
 * *  nph 09/14
 */

#include "APFEL/APFEL.h"
#include "APFEL/APFELdev.h"

#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"

#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"
#include "apfelcomb/fk_xgrid.h"

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
  
  static int nx = 0; // Number of initial scale x-points
  static double* xg = 0; // Initial scale x-grid

  // Set modes
  void setDISmode(bool const& mode) {DIS_mode=mode;};
  void setFTDYmode(bool const& mode) {FTDY_mode=mode;};
  void setSIAmode(bool const& mode) {SIA_mode=mode;};

  // Return number of x-points in grid
  size_t getNXGrid() { return (size_t) nx; };
  
  // Return X value (APFEL coordinates)
  double getXVal( const int idx) { return xg[idx]; };

  // Return APFEL proton mass
  double getProtonMass()
  {
    return APFEL::GetProtonMass();
  }

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

    // Initial scale hardcoded

    cout << "                FastKernel Grid Combination                 "<<endl;
    cout << "    - TheoryID: "<<param.thID << endl;
    cout << "    - Conv PTOrd: N"<<param.pto -1<<"LO"<<endl;

    return;
  }


  void set_params(qcd_param const& par, NNPDF::FKHeader& FK)
  {
    // Set FK table description
    FK.AddTag(FKHeader::VERSIONS, "APFEL", APFEL::GetVersion());
    FK.AddTag(FKHeader::VERSIONS, "libnnpdf", NNPDF::getVersion());

    FK.AddTag(FKHeader::GRIDINFO, "NX", nx);

    // x-grid header
    stringstream xGheader;
    for (int i=0; i<nx; i++)
      xGheader << std::setprecision(16) << std::scientific << xg[i] <<std::endl;

    FK.AddTag(FKHeader::BLOB, "xGrid", xGheader.str());

    // Full flavourmap
    stringstream fMapHeader;
    for (int i=0; i<14; i++)
    {
      for (int i=0; i<14; i++)
        fMapHeader << "1 ";
      fMapHeader<<std::endl;
    }
    FK.AddTag(FKHeader::BLOB, "FlavourMap", fMapHeader.str());

    // Theory configuration
    std::map<std::string, std::string>::const_iterator imap;
    for (imap = par.thMap.begin(); imap != par.thMap.end(); imap++)
      FK.AddTag(FKHeader::THEORYINFO, imap->first, imap->second);
  }
  
  // *********************** EVOLUTON FUNCTIONS *****************************

  void initTimelike()
  {
    APFEL::SetPDFSet("MELA");
    APFEL::SetTimeLikeEvolution(true);
  }


  // Initialise QCD according to parameters
  void initQCD(qcd_param const& par, const double& Q2max)
  {
    APFEL::SetParam(par.thMap);

    // Init Q0
    QCD::Q0 = par.Q0;

    // Truncated Epsilon
    APFEL::SetEpsilonTruncation(1E-1);

    // Set maximum scale
    QM = sqrt(Q2max);
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
  void initEvolgrid(int const& _nx, double const& xmin)
  {
    // Reset cache
    QC = 0;
    
    // Re-init
    nx = _nx;
    if (xg != 0)
      delete[] xg;
    xg = new double[nx+1];

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
    
    return;
  }

  void initNNPDF30Grid()
  {
    // Reset cache
    QC = 0;
    nx = 50;

    if (xg != 0)
      delete[] xg;
    xg = new double[nx];

    xg[0] = 1E-07;
    xg[1] = 1.73780082874938E-07;
    xg[2] = 3.01995172040202E-07;
    xg[3] = 5.24807460249773E-07;
    xg[4] = 9.1201083935591E-07;
    xg[5] = 1.58489319246111E-06;
    xg[6] = 2.75422870333817E-06;
    xg[7] = 4.78630092322639E-06;
    xg[8] = 8.31763771102671E-06;
    xg[9] = 1.44543977074593E-05;
    xg[10] = 2.51188643150958E-05;
    xg[11] = 4.36515832240166E-05;
    xg[12] = 7.58577575029184E-05;
    xg[13] = 0.000131825673855641;
    xg[14] = 0.000229086765276777;
    xg[15] = 0.000398107170553497;
    xg[16] = 0.000691830970918937;
    xg[17] = 0.00120226443461741;
    xg[18] = 0.00208929613085404;
    xg[19] = 0.00363078054770101;
    xg[20] = 0.00630957344480194;
    xg[21] = 0.0109647819614318;
    xg[22] = 0.0190546071796325;
    xg[23] = 0.0331131121482591;
    xg[24] = 0.0575439937337157;
    xg[25] = 0.1;
    xg[26] = 0.1375;
    xg[27] = 0.175;
    xg[28] = 0.2125;
    xg[29] = 0.25;
    xg[30] = 0.2875;
    xg[31] = 0.325;
    xg[32] = 0.3625;
    xg[33] = 0.4;
    xg[34] = 0.4375;
    xg[35] = 0.475;
    xg[36] = 0.5125;
    xg[37] = 0.55;
    xg[38] = 0.5875;
    xg[39] = 0.625;
    xg[40] = 0.6625;
    xg[41] = 0.7;
    xg[42] = 0.7375;
    xg[43] = 0.775;
    xg[44] = 0.8125;
    xg[45] = 0.85;
    xg[46] = 0.8875;
    xg[47] = 0.925;
    xg[48] = 0.9625;
    xg[49] = 1;

    // Set evolution operator parameters
    APFEL::EnableWelcomeMessage(false);
    APFEL::SetFastEvolution(false);
    APFEL::EnableEvolutionOperator(true); // Enable the computation of the Evolution Operator
    APFEL::SetNumberOfGrids(1);           // The evolution will be done on a number of subgrids defined here
    APFEL::SetExternalGrid(1,nx-1,5,xg);    // Set the grid as external (np: number of intervals, 3: interpolation degree, xg: defined above)

    // Start APFEL
    if (DIS_mode)
      APFEL::InitializeAPFEL_DIS();
    else
      APFEL::InitializeAPFEL();
    
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
    if (xg != NULL && x<xg[0])
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
    if (xg != NULL && x<xg[0])
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

