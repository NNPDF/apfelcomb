// Classes for A matrices (combined evolution-rotation)
// and Sigma matrices (combined applgrid - A matrices)
// n.p.hartland@ed.ac.uk  - 03/12

#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include "NNPDF/common.h"
#include "NNPDF/nnpdfdb.h"
#include "NNPDF/commondata.h"

#include <chrono>
#include <ctime>

using namespace std;
using NNPDF::FKHeader;

namespace APP
{
  vector<std::string> splitpdf ( std::string const& str )
  {
    vector<std::string> outvec;
    std::stringstream ss(str);
    std::string s;

    while (getline(ss, s, ':')) 
      outvec.push_back(s);
    return outvec;
  }

  // ******************* APPLGrid parsing ******************************

  // Determine the minimum x-point required in a target grid
  double parse_xmin(std::vector<int> const& subgridIDs)
  {
    double xmin = 1;
    for ( const int i : subgridIDs )
    {
      APP::appl_param par; APP::parse_input(i, par, true);
      APP::grid* subgrid = new APP::grid(par);
      xmin = std::min(xmin, APP::getXmin(subgrid->g,true));
      delete subgrid;
    }

    return xmin;
  }

  // Determine the total number of datapoints in a target grid
  int parse_ndat(std::vector<int> const& subgridIDs)
  {
    int nDat = 0; 
    for (const int i: subgridIDs)
    {
      APP::appl_param par; APP::parse_input(i, par, true);
      nDat += par.incdat + par.muldat*par.ndata;
    }
    return nDat;
  }

  // Construct applgrid<->datafile map for target subgrid ID
  std::vector< std::vector<int> > parse_map(std::vector<int> const& subgridIDs, int const& tID)
  {
    int lastNdat = 0; 
    for (const int i: subgridIDs)
    {
      APP::appl_param par; APP::parse_input(i, par, true);
      lastNdat += par.incdat;
      if (i == tID)
      {
        std::vector< std::vector<int> > applDataMap(par.ndata, vector<int>(par.muldat, -1));
        for (size_t j=0; j<par.ndata; j++)
          for (size_t k=0; k<par.muldat; k++)
            applDataMap[j][k] = lastNdat +  j*par.muldat + k;
        return applDataMap;
      }
      lastNdat += par.muldat*par.ndata;
    }

    std::cerr << "Error: target ID " << tID << " not present in list of subgrids!"<<std::endl;
    exit(-1);
  }

  void parse_input(int innum, appl_param& param, bool silent)
  {
    // Setup db connection
    NNPDF::IndexDB grid_db(databasePath()+"applgrid.db", "grids");
    NNPDF::IndexDB subgrid_db(databasePath()+"applgrid.db", "subgrids");

    // Fetch number of entries
    const int entries =subgrid_db.GetNEntries();
    if (innum < 0 || innum > entries)
    {
      cerr << "Error: applgrid ID ("<<innum<<") must be between 1 and "<<entries<<endl;
      exit(-1);
    }

    // Read grid information
    const std::string fktarget = NNPDF::dbquery<string>(subgrid_db,innum,"fktarget");
    const int target = NNPDF::dbmatch(grid_db, "name", fktarget)[0];
    param.nx      =  NNPDF::dbquery<int>(grid_db,target,"nx");
    param.desc    =  NNPDF::dbquery<string>(grid_db,target,"description");
    param.setname =  NNPDF::dbquery<string>(grid_db,target,"setname");
    
    // Read subgrid information
    param.applgrid  = applPath() + param.setname + "/" + NNPDF::dbquery<string>(subgrid_db,innum,"applgrid");
    param.fnlobin =  NNPDF::dbquery<int>(subgrid_db,innum,"fnlobin");
    param.ptmin   =  NNPDF::dbquery<int>(subgrid_db,innum,"ptmin");
    param.pdfwgt    = NNPDF::dbquery<bool>(subgrid_db,innum,"pdfwgt");
    param.ppbar     = NNPDF::dbquery<bool>(subgrid_db,innum,"ppbar");
    param.common_subgrids = NNPDF::dbmatch(subgrid_db, "fktarget", fktarget);
    param.gridname  = fktarget + "_" + to_string(innum) + ".subgrid";


    // Fetch datapoint mask
    string mask = NNPDF::dbquery<string>(subgrid_db,innum,"mask");
    vector<string> masksplit = ssplit(mask);
    for (size_t i=0; i<masksplit.size(); i++)
      param.mask.push_back((bool) atoi(masksplit[i].c_str()));
    param.ndata   =   std::count(param.mask.begin(), param.mask.end(), 1);

    // Fill map
    param.maskmap.clear();
    for (size_t i=0; i<param.mask.size(); i++)
      if (param.mask[i]==true)
        param.maskmap.push_back(i);

    // Read operators
    param.incdat = 0; param.muldat = 1; param.nrmdat = 1.0;
    const std::string operators = NNPDF::dbquery<string>(subgrid_db,innum,"operators");
    const bool activeOperator = operators != "";
    if (activeOperator)
    {
      const vector<string> tokens = gsplit(operators, ",");
      for (auto s:tokens) 
      {
        const vector<string> subtoken = gsplit(s, ":");
        if (subtoken.size() != 2)
        {
          std::cerr << "Cannot parse operator: " <<s <<std::endl;
          exit(-1);
        }

        if (subtoken[0] == "*") param.muldat = atoi(subtoken[1].c_str());
        if (subtoken[0] == "+") param.incdat = atoi(subtoken[1].c_str());
        if (subtoken[0] == "N") param.nrmdat = atof(subtoken[1].c_str());
      }
    }

    // See if there is a valid README
    string readmefilename  = applPath() + param.setname + "/README_appl_" + param.setname;
    std::ifstream readmefile; readmefile.open(readmefilename.c_str());
    if (!readmefile.good())
    {
      std::cerr << "APPL::parse_input Error: Cannot read README from: "<<std::endl<<readmefilename<<std::endl;
      exit(-1);
    }

    string readme((std::istreambuf_iterator<char>(readmefile)),
                std::istreambuf_iterator<char>());
    param.readme = readme;

    // Get common grids for inventory
    vector<int> commonGrids = NNPDF::dbmatch(grid_db, "setname", param.setname);
    for ( auto i : commonGrids)
    {
      const std::string commonTarget = NNPDF::dbquery<string>(grid_db,i,"name");
      vector<int> commonSubgrids = NNPDF::dbmatch(subgrid_db, "fktarget", commonTarget);
      for ( auto j : commonSubgrids) param.inventory.push_back(fktarget+"_"+to_string(j)+".subgrid");
    }

    /*
     *        ***    VERIFICATION    ***
     */
    
    if (param.ndata<1)
    {
      cerr <<"Error: invalid ndata: "<<param.ndata<<endl;
      exit(-1);
    }
    
    if (param.muldat != 1 && param.ndata != 1)
    {
      cerr << "Error: cannot use multiplicative operator unless there is only one datapoint" <<endl;
      cerr << param.muldat <<"  "<<param.ndata <<endl;
      exit(-1);
    }
    
    if (param.ptmin >= param.pto && !silent )
      cout << "Warning: minimum perturbative order is greater than the maximum perturbative order!"<<endl;

    // Set number of active ptords
    param.pto -= param.ptmin;
    param.pto = std::max(param.pto, (size_t) 1);

    if (!silent)
    {
      cout <<endl;
      cout << "    - OutputFile: "<<param.gridname<<endl;
      cout << "    - SetName: "<<param.setname<<endl;
      cout <<endl;
      cout << "    - PTMin: "<<param.ptmin<<endl;
      cout << "    - NData: "<<param.ndata<<endl;
      cout << "    - fnlobin: " <<param.fnlobin <<endl;
      cout <<endl;
      cout << "    - PDFWeight: "<<param.pdfwgt<<endl;
      cout << "    - ppbar transform: "<<param.ppbar<<endl;
      cout <<endl;
      if (activeOperator)
      {
        cout << "    - data increment: "<<param.incdat<<endl;
        cout << "    - data multiplier: "<<param.muldat<<endl;
        cout << "    - data normalisation: "<<param.nrmdat<<endl;
        cout <<endl; 
      }
      cout << "    - Common grids: "<<param.inventory.size()<<endl;
      cout <<endl;
      DisplayHR();
    }
    
    return;
  }

// *********************** FKHeader population ************************

  void set_params(appl_param const& par, NNPDF::FKHeader& FK)
  {
      FK.AddTag(FKHeader::BLOB, "GridDesc", par.desc);
      FK.AddTag(FKHeader::BLOB, "Readme", par.readme);
      FK.AddTag(FKHeader::GRIDINFO, "SETNAME", par.setname);
      FK.AddTag(FKHeader::GRIDINFO, "NDATA", parse_ndat(par.common_subgrids));
      FK.AddTag(FKHeader::GRIDINFO, "HADRONIC", true);
      FK.AddTag(FKHeader::VERSIONS, "APPLrepo", applCommit() );

      // Full flavourmap
      stringstream fMapHeader;
      for (int i=0; i<14; i++)
      {
        for (int i=0; i<14; i++)
          fMapHeader << "1 ";
        fMapHeader<<std::endl;
      }
      FK.AddTag(FKHeader::BLOB, "FlavourMap", fMapHeader.str());

      // Set QCD parameters
      QCD::set_params(par, FK);
  }

  // *********************** APPLgrid helpers ****************************

  // Get the maximum scale of an applgrid
  double getQ2max(const appl::grid* g)
  {
    // Find maximum required scale
    double Q2max = g->weightgrid(0, 0)->getQ2max();
    for(int i=0; i<2; i++)  // pto
      for (int j=0; j<g->Nobs(); j++) // bin
      {
        appl::igrid const *igrid = g->weightgrid(i, j);
        Q2max = max(Q2max, igrid->getQ2max());
      }
    return Q2max;
  }


  // Translates 'loop' order to appl::grid index
  // This is specifically in order to translate aMC@NLO-like four-part grids
  // into LO and NLO components.
  // aMC@NLO-like convolution uses Born = grid 3
  //                               NLO  = grid 0
  int get_grid_idx( appl::grid *g, int const& pto )
  {
    if (g->calculation() == appl::grid::AMCATNLO)
      return (pto==0) ? 3:0;
    return pto;
  }

  // Returns the APPLgrid PDF object associated with the ith subgrid of g
  appl::appl_pdf* get_appl_pdf( appl::grid *g, int const& i )
  {
    const std::string pdfnames = g->getGenpdf();
    std::vector<std::string> pdfvec = splitpdf( pdfnames );
    const size_t isubproc = pdfvec.size() == 1 ? 0:i;
    return appl::appl_pdf::getpdf( pdfvec[isubproc] );
  }

  // Returns the minimum and maximum x-grid points for a specified subgrid slice.
  // igrid is the requested subgrid, nsubproc the number of subprocesses held within igrid.
  // tau specified the bin in scale to be investigated, and alpha specifies the bin in x1.
  // return.first and return.second return the minumum and maxiumum bins in x2 respectively.
  std::pair<int,int> get_slice_limits(appl::igrid const* igrid, int const& nsubproc, int const& tau, int const& alpha)
  {
    std::pair<int,int> limits(igrid->Ny2(), 0);
    appl::igrid* igrid_nc = const_cast<appl::igrid*>(igrid);
    for (int tsp=0; tsp<nsubproc; tsp++)
      if ((*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))[tau] != NULL)
        for (int ix2=0; ix2<igrid->Ny2(); ix2++)
          if ( (*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))(tau,alpha,ix2) != 0 )
          {                
            limits.first  = std::min(ix2, limits.first);
            limits.second = std::max(ix2, limits.second);
          }

    return limits;
  }

  // Returns the minimum and maximum used values of x1 across all grid slices
  std::pair<int,int> get_igrid_limits_x1(appl::igrid const* igrid, int const& nsubproc, int const& tau)
  {
    std::pair<int,int> limits(igrid->Ny1(), 0);
    for (int alpha = 0; alpha < igrid->Ny1(); alpha++)
    { 
      const std::pair<int,int> sl = get_slice_limits(igrid, nsubproc, tau, alpha);  
      if (sl.first <= sl.second) // This alpha is not trimmed  
      {
        limits.first  = std::min(alpha, limits.first);
        limits.second = std::max(alpha, limits.second);
      }
    }
    return limits;
  }

  // Returns the minimum and maximum used values of x2 across all grid slices
  std::pair<int,int> get_igrid_limits_x2(appl::igrid const* igrid, int const& nsubproc, int const& tau)
  {
    std::pair<int,int> limits(igrid->Ny2(), 0);
    for (int alpha = 0; alpha < igrid->Ny1(); alpha++)
    {
      const std::pair<int,int> sl = get_slice_limits(igrid, nsubproc, tau, alpha);  
      limits.first  = std::min(sl.first, limits.first);
      limits.second = std::max(sl.second, limits.second);
    }   
    return limits;
  }

  // If nonzero is true, get the smallest x-value regardless of weight
  // if false, get the smallest x-value associated with a non-zero weight
  double getXmin(const appl::grid* g, const bool& nonzero)
  {
    // Run over all grids to verify kinematic limits
    double xmin = 1.0;
    for(int i=0; i<=g->nloops(); i++) 
      for (int j=0; j<g->Nobs(); j++)
      {
        appl::igrid const *igrid = g->weightgrid(i, j);
        const size_t nsubproc = g->subProcesses(i);

        if (!nonzero)
        {
          xmin = min(xmin, igrid->fx(igrid->gety1(0)));
          xmin = min(xmin, igrid->fx(igrid->gety1(igrid->Ny1()-1)));
          xmin = min(xmin, igrid->fx(igrid->gety2(0)));
          xmin = min(xmin, igrid->fx(igrid->gety2(igrid->Ny2()-1)));
        }
        else
        {
          for (int t=0; t<igrid->Ntau(); t++) // Loop over Q^2 integral
          {
            const std::pair<int,int> l1 = get_igrid_limits_x1(igrid, nsubproc, t);  
            if (l1.first <= l1.second)
            {
              xmin = min(xmin, igrid->fx(igrid->gety1(l1.first)));
              xmin = min(xmin, igrid->fx(igrid->gety1(l1.second)));
            }
            const std::pair<int,int> l2 = get_igrid_limits_x2(igrid, nsubproc, t);  
            if (l2.first <= l2.second)
            {
              xmin = min(xmin, igrid->fx(igrid->gety2(l2.first)));
              xmin = min(xmin, igrid->fx(igrid->gety2(l2.second)));
            }
          }
        }
      }
    
    return xmin;
  }

  // Computes the normalisation required for translating APPLgrid weights to FK ones
  // e.g including factors of alpha_S, bin width etc.
  // g is the APPLgrid being combined with evolution factors.
  // pto specified the perturbative order being combined, as the value of alpha_S in the current bin,
  // and x1/x2 specify the numerical values of the PDF x-values for the first and second PDF respectively.
  double compute_wgt_norm(appl::grid *g, appl_param const& par, int const& d, double const& pto, double const& as, double const& x1, double const& x2)
  {
    // PDF x and bin width normalisation
    double norm = 1.0/(x1*x2*g->deltaobs(d));

    // Normalisation by number of runs
    appl::grid* g_nc = const_cast<appl::grid*>(g);
    if ( !g->getNormalised() && g_nc->run() )
      norm*=1.0/double(g_nc->run());

    // Factor of alpha_S
    const double LO = g->leadingOrder();
    if (g->calculation() == appl::grid::AMCATNLO)
    { norm*=pow(as*(4.0*M_PI), LO+pto ); }
    else
    { norm*=pow(as/(2.0*M_PI), LO+pto ); }
    
    return norm;
  }

   // Progress update ****************************************************

  int countElements(appl_param const& par, appl::grid* g)
  {
    // Counter
    int nElm = 0;
    for (auto bin : par.maskmap )
      for (size_t pto=0; pto<par.pto; pto++) 
      {
        int gidx = get_grid_idx(g, pto); // APPLgrid grid index
        appl::igrid const *igrid = g->weightgrid(gidx, bin);
        const size_t nsubproc = g->subProcesses(gidx);     // Number of subprocesses
        for (int t=0; t<igrid->Ntau(); t++)     // Loop over scale bins
          for (int a=0; a<igrid->Ny1(); a++  )  // Loop over x1 bins
          {
            const std::pair<int,int> sl = get_slice_limits(igrid, nsubproc, t, a);  
            nElm += std::max(0, sl.second - sl.first + 1);
          }        
      }
    return nElm;
  } 

  void statusUpdate( time_point const& t1, int const& totElements, int& compElements)
  {
    // Increment computed elements
    compElements++;
    const int interval = compElements< 100 ? 5:10;
    if (compElements % interval == 0)
    {
      // Elapsed time update
      const time_point t2 = std::chrono::system_clock::now();
      const double percomp= 100.0*((double)compElements/(double)totElements);
      const time_point t3 = t2 + std::chrono::duration_cast<time_span>( (t2-t1) * ( 100.0 / percomp - 1.0 ));
      const std::time_t end_time = std::chrono::system_clock::to_time_t(t3);

      char eta[80]; strftime (eta,80,"ETA: %R %x.", localtime(&end_time));
      cout << "-- "<< setw(6) << setprecision(4)  << percomp << "\% complete. "
           << eta <<"\r";
      cout.flush();
    }
  }


// ********************* Evolution factors ******************************

  class EvolutionFactors
  {
  public:
    EvolutionFactors(const int nxin, const int nxout):
    b1(13),
    b2(13*14),
    b3(13*14*nxin),
    data(new double[nxout*b3]) 
    {
      for (int i=0; i<nxout*b3; i++)
        data[i] = 0;
    };
    ~EvolutionFactors() {delete[] data;};

    double* operator()(int const& ox, int const& ix, int const& fi ) {return data+b3*ox+b2*ix+b1*fi;};
    const double* operator()(int const& ox, int const& ix, int const& fi ) const {return data+b3*ox+b2*ix+b1*fi;};
  private:
    const int b1;
    const int b2;
    const int b3;
    double* data;
  };

  // ******************* FK Table computation ****************************

  void computeFK(appl_param const& par, appl::grid* g, NNPDF::FKGenerator* fk)
  {    
    // Progress monitoring
    int completedElements = 0;
    const int nXelements = countElements(par, g);
    const time_point t1 = std::chrono::system_clock::now();

    const vector<size_t> afl = QCD::active_flavours(par);
    std::cout << "APFELcomb: "<<afl.size() <<" active flavours in evolution." <<std::endl;
   
    for (size_t d=0; d<par.ndata; d++)
    {    
      // Fetch associated applgrid info
      const size_t bin = par.maskmap[d];
      for (size_t pto=0; pto<((size_t) par.pto); pto++) // Loop over perturbative order
      {
        // Determine grid index, allocate subprocess arrays
        const int gidx = get_grid_idx(g, pto+par.ptmin);
        appl::appl_pdf *genpdf = get_appl_pdf( g, gidx );
        const size_t nsubproc = g->subProcesses(gidx);
        double *W = new double[nsubproc];   // Weight array
        double *H = new double[nsubproc];   // Evolution factor array
        double *H1 = new double[nsubproc];  // Split evolution factor array 1
        double *H2 = new double[nsubproc];  // Split evolution factor array 2

        // Fetch grid pointer and loop over Q
        appl::igrid const *igrid = g->weightgrid(gidx, bin);
        for (int t=0; t<igrid->Ntau(); t++) 
        {
          // Scales and strong coupling
          const double Q2  = igrid->fQ2( igrid->gettau(t));
          const double QF  = sqrt(Q2)*par.xiF;
          const double QR  = sqrt(Q2)*par.xiR;
          const double as  = QCD::alphas(QR);

          // Renormalisation and factorisation scale variation terms
          const bool vary_ren = pto == 0 && par.evol_pto == 1 && par.xiR != 1.0;
          const bool vary_fac = pto == 0 && par.evol_pto == 1 && par.xiF != 1.0;
          const double renscale =  (as/(2.0*M_PI))*2.0*M_PI*QCD::beta0()*g->leadingOrder()*log(par.xiR*par.xiR);
          const double facscale = -(as/(2.0*M_PI))*log(par.xiF*par.xiF);

          // define evolution factor arrays
          const int nxin = fk->GetNx();
          EvolutionFactors fA1(nxin, igrid->Ny1());   // PDF 1 Evolution factors
          EvolutionFactors fA2(nxin, igrid->Ny2());   // PDF 2 Evolution factors
          EvolutionFactors fdA1(nxin, igrid->Ny1());  // PDF 1 Splitting function x Evolution factors
          EvolutionFactors fdA2(nxin, igrid->Ny2());  // PDF 2 Splitting function x Evolution factors

          // Compute nonzero evolution factors
          const std::pair<int,int> l1 = get_igrid_limits_x1(igrid, nsubproc, t);  
          const std::pair<int,int> l2 = get_igrid_limits_x2(igrid, nsubproc, t);  

          for (size_t ix = 0; ix < nxin; ix++)
          for (size_t fl : afl)
          {
            for (int ox=l1.first; ox<=l1.second; ox++) // Loop over applgrid x1
            {
              QCD::EvolutionOperator(false, ix,igrid->fx(igrid->gety1(ox)),fl,QF,fA1(ox,ix,fl));
              if (vary_fac) QCD::DerivativeOperator(false, ix,igrid->fx(igrid->gety1(ox)),fl,QF,fdA1(ox,ix,fl));
            }
            for (int ox=l2.first; ox<=l2.second; ox++) // Loop over applgrid x2
            {
              QCD::EvolutionOperator(par.ppbar, ix,igrid->fx(igrid->gety2(ox)),fl,QF,fA2(ox,ix,fl));
              if (vary_fac) QCD::DerivativeOperator(par.ppbar, ix,igrid->fx(igrid->gety2(ox)),fl,QF,fdA2(ox,ix,fl));
            }
          }

          for (int a=l1.first; a<=l1.second; a++ )
          {
            const double x1 = igrid->fx(igrid->gety1(a));           
            const std::pair<int,int> limits = get_slice_limits(igrid, nsubproc, t, a);  
            for (int b=limits.first; b<=limits.second; b++) // Loop over applgrid x2
            {
              // fetch weight values
              for (size_t ip=0; ip<nsubproc; ip++)
                W[ip] = (*(const SparseMatrix3d*) const_cast<appl::igrid*>(igrid)->weightgrid(ip))(t,a,b);
              
              // Calculate normalisation factors
              const double x2 = igrid->fx(igrid->gety2(b));
              const double pdfnrm = par.pdfwgt ? igrid->weightfun(x1)*igrid->weightfun(x2) : 1.0;
              const double norm = pdfnrm*compute_wgt_norm(g, par, bin, pto+par.ptmin, as, x1, x2)*par.nrmdat;

              for (size_t i=0; i<nxin; i++)    // Loop over input pdf x1
                for (size_t j=0; j<nxin; j++)  // Loop over input pdf x2
                  for (size_t k : afl)         // loop over flavour 1
                    for (size_t l : afl)       // loop over flavour 2
                    {
                      // Rotate to subprocess basis
                      genpdf->evaluate(fA1(a,i,k),fA2(b,j,l),H);

                      if (vary_fac)
                      {
                        genpdf->evaluate(fdA1(a,i,k),fA2(b,j,l),H1);
                        genpdf->evaluate(fA1(a,i,k),fdA2(b,j,l),H2);
                      }

                      for (size_t ip=0; ip<nsubproc; ip++)
                        if (W[ip] != 0)
                        {
                          double fill = H[ip];                              // Basic fill
                          if (vary_ren) fill += renscale*H[ip];             // Ren. scale variation
                          if (vary_fac) fill += facscale*(H1[ip] + H2[ip]); // Fac. scale variation
                          if ( fill != 0.0 )
                            for (int const& td : par.datamap[d])
                              fk->Fill( td, i, j, k, l, norm*fill*W[ip]);
                        }
                    }

              // Update progress
              statusUpdate(t1, nXelements, completedElements);
            }
          }
        }

        // Free subprocess arrays
        delete[] W;
        delete[] H;
        delete[] H1;
        delete[] H2;

      } // /pto
    } // /data

      
    cout << "FastKernel table computed."<<endl;
    
    return;

  }
}