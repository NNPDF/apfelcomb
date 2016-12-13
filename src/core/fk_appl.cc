// Classes for A matrices (combined evolution-rotation)
// and Sigma matrices (combined applgrid - A matrices)
// n.p.hartland@ed.ac.uk  - 03/12

#include "apfelcomb/fk_appl.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include <sys/time.h>

#include <NNPDF/common.h>
#include "NNPDF/nnpdfdb.h"

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

  void parse_input(int innum, appl_param& param)
  {
    // Setup db connection
    NNPDF::IndexDB db(databasePath()+"applgrid.db", "sets");

    // Fetch number of entries
    const int entries =db.GetNEntries();
    if (innum < 0 || innum > entries)
    {
      cerr << "Error: applgrid ID ("<<innum<<") must be between 1 and "<<entries<<endl;
      exit(-1);
    }

    // Set/grid names
    param.gridname  = NNPDF::dbquery<string>(db,innum,"gridname");
    param.setname   = NNPDF::dbquery<string>(db,innum,"setname");
    param.gridfile  = applPath() + param.setname + "/" + NNPDF::dbquery<string>(db,innum,"gridfile");

    // Set up values
    param.nx      =  NNPDF::dbquery<int>(db,innum,"nxpt");
    param.ndata   =  NNPDF::dbquery<int>(db,innum,"ndata");
    param.nbins   =  NNPDF::dbquery<int>(db,innum,"nbins");
    param.fnlobin =  NNPDF::dbquery<int>(db,innum,"fnlobin");
    param.ptmin   =  NNPDF::dbquery<int>(db,innum,"ptmin");

    param.xmin    =  NNPDF::dbquery<double>(db,innum,"xmin");
    param.maxprec =  NNPDF::dbquery<double>(db,innum,"maxprec");

    param.setname   = NNPDF::dbquery<string>(db,innum,"setname");
    param.desc      = NNPDF::dbquery<string>(db,innum,"description");

    param.fnlo      = NNPDF::dbquery<bool>(db,innum,"fnlo");
    param.pdfwgt    = NNPDF::dbquery<bool>(db,innum,"pdfwgt");
    param.ppbar    = NNPDF::dbquery<bool>(db,innum,"ppbar");

    // Fetch datapoint mask
    string mask = NNPDF::dbquery<string>(db,innum,"mask");
    vector<string> masksplit = ssplit(mask);
    for (size_t i=0; i<masksplit.size(); i++)
      param.mask.push_back((bool) atoi(masksplit[i].c_str()));
    
    // Fill map
    param.map.clear();
    for (size_t i=0; i<param.mask.size(); i++)
      if (param.mask[i]==true)
        param.map.push_back(i);

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

    // Get common grids
    vector<int> commonGrids = NNPDF::dbmatch(db, "setname", param.setname);
    for ( auto i : commonGrids)
      param.inventory.push_back(NNPDF::dbquery<string>(db,i,"gridname"));

    /*
     *        ***    VERIFICATION    ***
     */
    
    if (param.ndata<1)
    {
      cerr <<"Error: invalid ndata: "<<param.ndata<<endl;
      exit(-1);
    }
    
    if (param.nbins<1)
    {
      cerr <<"Error: invalid nbins: "<<param.nbins<<endl;
      exit(-1);
    }
    
    if (param.ndata>param.nbins)
    {
      cerr <<"Error: ndata: "<<param.ndata<<" > nbins: "<<param.nbins<<endl;
      exit(-1);
    }
    
    if (param.mask.size()!=param.nbins)
    {
      cerr << "Error: mask size: "<<param.mask.size()<<" does not match nbins: "<<param.nbins<<endl;
      exit(-1);
    }
    
    size_t count=0;
    for (size_t i=0; i<param.mask.size(); i++)
      if (param.mask[i])
        count++;
    
    if (count!=param.ndata)
    {
      cerr << "Error: mask does not match ndata: "<<param.ndata<<endl;
      exit(-1);
    }
    
    if (param.ptmin >= param.pto)
      cout << "Warning: minimum perturbative order is greater than the maximum perturbative order!"<<endl;

    // Set number of active ptords
    param.pto -= param.ptmin;
    param.pto = std::max(param.pto, (size_t) 1);
    cout <<endl;
    cout << "    - Gridname: "<<param.gridname<<endl;
    cout << "    - SetName: "<<param.setname<<endl;
    cout <<endl;
    cout << "    - PTMin: "<<param.ptmin<<endl;
    cout << "    - NBins: "<<param.nbins<<endl;
    cout << "    - NData: "<<param.ndata<<endl;
    cout <<endl;
    cout << "    - PDFWeight: "<<param.pdfwgt<<endl;
    cout << "    - ppbar transform: "<<param.ppbar<<endl;
    cout <<endl;
    cout << "    - Common grids: "<<param.inventory.size()<<endl;
    cout <<endl;
    DisplayHR();
    
    return;
  }

// *********************** FKHeader population ************************

  void set_params(appl_param const& par, NNPDF::FKHeader& FK)
  {
      FK.AddTag(FKHeader::BLOB, "GridDesc", par.desc);
      FK.AddTag(FKHeader::BLOB, "Readme", par.readme);
      FK.AddTag(FKHeader::GRIDINFO, "SETNAME", par.setname);
      FK.AddTag(FKHeader::GRIDINFO, "NDATA", par.ndata);
      FK.AddTag(FKHeader::GRIDINFO, "HADRONIC", true);

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

  // If nonzero is true, get the smallest x-value regardless of weight
  // if false, get the smallest x-value associated with a non-zero weight
  double getXmin(const appl::grid* g, const bool& nonzero)
  {
    // Run over all grids to verify kinematic limits
    double xmin = 1.0;
    
    for(int i=0; i<2; i++)  // pto
      for (int j=0; j<g->Nobs(); j++)
      {
        appl::igrid const *igrid = g->weightgrid(i, j);
        const size_t nsubproc = g->subProcesses(0);

        for (int ix1=0; ix1<igrid->Ny1(); ix1++)
          for (int ix2=0; ix2<igrid->Ny2(); ix2++)
            for (int t=0; t<igrid->Ntau(); t++) // Loop over Q^2 integral
              for (size_t ip=0; ip<nsubproc; ip++)
                {
                  // Associated weight
                  const bool zero_weight = (*(const SparseMatrix3d*) const_cast<appl::igrid*>(igrid)->weightgrid(ip))(t,ix1,ix2) == 0;

                  if (!zero_weight || !nonzero) 
                  {
                    xmin = min(xmin, igrid->fx(igrid->gety1(ix1)));
                    xmin = min(xmin, igrid->fx(igrid->gety2(ix2)));
                  }
                }
      }
    
    return xmin;
  }

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

  // Returns the minimum and maxiumum x-grid points for a specified subgrid slice.
  // igrid is the requested subgrid, nsubproc the number of subprocesses held within igrid.
  // tau specified the bin in scale to be investigated, and alpha specifies the bin in x1.
  // nxlow and nxhigh return the minumum and maxiumum bins in x2 respectively.
  void get_igrid_limits(appl::igrid const* igrid, int const& nsubproc, int const& tau, int const& alpha, int& nxlow, int& nxhigh)
  {
    // Need to remove const due to APPLgrid non-const methods
    appl::igrid* igrid_nc = const_cast<appl::igrid*>(igrid);
    nxlow=igrid->Ny2(); nxhigh=0;
    tsparse1d<double> *ts1;
    for (int tsp=0; tsp<nsubproc; tsp++)
      if ((*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))[tau] != NULL)
        if (( ts1 = (*(const SparseMatrix3d*) igrid_nc->weightgrid(tsp))[tau]->trimmed(alpha) ) != NULL)
        {                
          nxlow  = std::min(ts1->lo(), nxlow);
          nxhigh = std::max(ts1->hi(), nxhigh);
        }
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
    else if (pto == 0 && par.evol_pto == 1) // NLO scale variation
    {
      const double pas = as/(2.0*M_PI);
      norm*=pow(pas, LO)*(1+pas*2.0*M_PI*QCD::beta0()*LO*log(par.xiR*par.xiR));
    }
    else
    { norm*=pow(as/(2.0*M_PI), LO+pto ); }
    
    return norm;
  }

  // Progress update ****************************************************

  int countElements(appl_param const& par, appl::grid* g)
  {
    // Counter
    int nElm = 0;
    for (auto bin : par.map )
      for (size_t pto=0; pto<par.pto; pto++) 
      {
        int gidx = get_grid_idx(g, pto); // APPLgrid grid index
        appl::igrid const *igrid = g->weightgrid(gidx, bin);
        const size_t nsubproc = g->subProcesses(gidx);     // Number of subprocesses
        for (int t=0; t<igrid->Ntau(); t++)     // Loop over scale bins
          for (int a=0; a<igrid->Ny1(); a++  )  // Loop over x1 bins
          {
            int nxlow, nxhigh; get_igrid_limits(igrid, nsubproc, t, a, nxlow, nxhigh);
            nElm += std::max(0,nxhigh - nxlow + 1);
          }        
      }
    return nElm;
  } 

  void statusUpdate(timeval const& t1, int const& totElements, int& compElements)
  {
    // Increment computed elements
    compElements++;

    const int interval = compElements< 100 ? 5:10;
    if (compElements % interval == 0)
    {
      // Elapsed time update
      timeval t2; gettimeofday(&t2, NULL);
      double elapsedTime = (t2.tv_sec - t1.tv_sec); 
      elapsedTime += (t2.tv_usec - t1.tv_usec) / 1E6f; 

      // Percentage complete, ETA
      const double percomp= 100.0*((double)compElements/(double)totElements);
      const double eta = ( elapsedTime / percomp ) * ( 100.0 - percomp );

      cout << "-- "<< setw(6) << setprecision(4)  << percomp << "\% complete."
           << " T Elapsed: "  << setw(6)<<setprecision(4) << elapsedTime/60.0 
           << " min. ETA: "   << setw(6)<<setprecision(4) << eta/60.0<<" min.\r";
      cout.flush();
    }
  }


// ********************* Evolution factors ******************************

  class EvolutionFactors
  {
  public:
    EvolutionFactors(const int _nxin):
    nxin(_nxin),
    data(new double[nxin*13*14]) {};
    ~EvolutionFactors() {delete[] data;};

    double* operator()(int const& ix, int const& fi ) {return data+13*14*ix+13*fi;};
  private:
    const int nxin;
    double* data;
  };

  // ******************* FK Table computation ****************************

  void computeFK(appl_param const& par, appl::grid* g, NNPDF::FKGenerator* fk)
  {    
    // Progress monitoring
    int completedElements = 0;
    const int nXelements = countElements(par, g);

     // define evolution factor arrays
    const int nxin = fk->GetNx();
    EvolutionFactors fA1(nxin);
    EvolutionFactors fA2(nxin);
    
    // Begin progress timer
    timeval t1;
    gettimeofday(&t1, NULL);
   
    for (size_t d=0; d<fk->GetNData(); d++)
    {    
      // Fetch associated applgrid info
      const size_t bin      = par.map[d];
      for (size_t pto=0; pto<((size_t) par.pto); pto++) // Loop over perturbative order
      {
        // Determine grid index, allocate subprocess arrays
        const int gidx = get_grid_idx(g, pto+par.ptmin);
        appl::appl_pdf *genpdf = get_appl_pdf( g, gidx );
        const size_t nsubproc = g->subProcesses(gidx);
        double *W = new double[nsubproc];
        double *H = new double[nsubproc];
        
        // Fetch grid pointer and loop over Q
        appl::igrid const *igrid = g->weightgrid(gidx, bin);
        for (int t=0; t<igrid->Ntau(); t++) 
        {
          const double Q2  = igrid->fQ2( igrid->gettau(t));
          const double QF  = sqrt(Q2)*par.xiF;
          const double QR  = sqrt(Q2)*par.xiR;
          const double as  = QCD::alphas(QR);
          
          for (int a=0; a<igrid->Ny1(); a++  )
          {
            const double x1 = igrid->fx(igrid->gety1(a));           
            int nxlow, nxhigh; get_igrid_limits(igrid, nsubproc, t, a, nxlow, nxhigh);

            // Compute nonzero evolution factors
            if (nxlow <= nxhigh) 
              for (size_t ix = 0; ix < nxin; ix++)
                for (size_t fl = 0; fl < 14; fl++)
                  QCD::avals(ix,x1,fl,QF,fA1(ix,fl));
            
            for (int b=nxlow; b<=nxhigh; b++) // Loop over applgrid x2
            {
              // Second x values
              const double x2 = igrid->fx(igrid->gety2(b));
              
              // fetch weight values
              bool nonzero=false;
              for (size_t ip=0; ip<nsubproc; ip++)
                if (( W[ip] = (*(const SparseMatrix3d*) const_cast<appl::igrid*>(igrid)->weightgrid(ip))(t,a,b) )!=0)
                  nonzero=true;
              
              // If nonzero, perform combination
              if (nonzero)
              {
                // Calculate normalisation factor
                const double pdfnrm =  par.pdfwgt ? igrid->weightfun(x1)*igrid->weightfun(x2) : 1.0;
                const double norm = pdfnrm*compute_wgt_norm(g, par, bin, pto+par.ptmin, as, x1, x2);

                // Compute evolution factors for second PDF
                for (size_t ix = 0; ix < nxin; ix++)
                  for (size_t fl = 0; fl < 14; fl++)
                  {
                    if (par.ppbar == true)
                      QCD::avals_pbar(ix,x2,fl,QF,fA2(ix,fl));
                    else
                      QCD::avals(ix,x2,fl,QF,fA2(ix,fl));
                  }
                
                for (size_t i=0; i<nxin; i++)    // Loop over input pdf x1
                  for (size_t j=0; j<nxin; j++)  // Loop over input pdf x2
                    for (size_t k=0; k<14; k++)         // loop over flavour 1
                      for (size_t l=0; l<14; l++)       // loop over flavour 2
                      {
                        // Rotate to subprocess basis
                        genpdf->evaluate(fA1(i,k),fA2(j,l),H);

                        for (size_t ip=0; ip<nsubproc; ip++)
                          if (W[ip] != 0 and H[ip] != 0)
                            fk->Fill( d, i, j, k, l, norm*W[ip]*H[ip] );
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

      } // /pto
    } // /data

      
    cout << "FastKernel table computed."<<endl;
    
    return;

  }
}