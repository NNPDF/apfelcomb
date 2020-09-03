// Classes for A matrices (combined evolution-rotation)
// and Sigma matrices (combined applgrid - A matrices)
// n.p.hartland@ed.ac.uk  - 03/12

#include "apfelcomb/fk_pine.h"
#include "apfelcomb/fk_qcd.h"
#include "apfelcomb/fk_utils.h"

#include "NNPDF/common.h"
#include "NNPDF/nnpdfdb.h"
#include "NNPDF/commondata.h"

#include <chrono>
#include <ctime>
#include <limits>
#include <stdexcept>

using namespace std;
using NNPDF::FKHeader;

namespace PINE
{

// // ********************* Evolution factors ******************************


  class EvolutionFactors
  {
  public:
    EvolutionFactors(const int nxin, const int nxout):
    b1(14),
    b2(14*14),
    b3(14*14*nxin),
    data(nxout*b3)
    {}

    double* operator()(size_t ox, size_t ix, size_t fi ) {return &data.at(b3*ox+b2*ix+b1*fi);}
    const double* operator()(size_t ox, size_t ix, size_t fi) const {return &data.at(b3*ox+b2*ix+b1*fi);}
  private:
    const size_t b1;
    const size_t b2;
    const size_t b3;
    std::vector<double> data;
  };

// ********************************* Basis rotation helpers *************************************

  // Rotates APFEL flavour basis into APPLgrid flavour basis (photon moves from 0 to 13)
  double evolpdf_pineapplgrid(int32_t pdgid, double x, double q2, void*)
  {
    if (x < APFEL::xGrid(0))  // A nice trick of APPLgrid is to request PDF x-values smaller than are actually used
      return 0;
    return QCD::flvpdf(pdgid, x, q2);
  }

  double alphas(double q2, void*)
  {
    return QCD::alphas(sqrt(q2));
  }

  // convert pineappl pdg to applgrid
  int32_t pdgid_to_apfel_array_id(int32_t pdgid)
  {
    switch (pdgid)
    {
    case 22:
      return 13;

    case 21:
      return 6;

    case -6:
    case -5:
    case -4:
    case -3:
    case -2:
    case -1:
    case 6:
    case 5:
    case 4:
    case 3:
    case 2:
    case 1:
      return pdgid + 6;

    default:
      throw std::runtime_error("pdgid not available.");
    }
  }

  // ********************************* Kinematics Helpers *************************************
/*
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


  // ********************************* Combination Helpers *************************************

  // Translates 'loop' order to appl::grid index
  // This is specifically in order to translate aMC@NLO-like four-part grids
  // into LO and NLO components.
  // aMC@NLO-like convolution uses Born = grid 3
  //                               NLO  = grid 0
  int get_grid_idx( const appl::grid *g, int const& pto )
  {
    if (g->calculation() == appl::grid::AMCATNLO)
      return (pto==0) ? 3:0;
    return pto;
  }

  // Returns the APPLgrid PDF object associated with the ith subgrid of g
  appl::appl_pdf* get_appl_pdf( const appl::grid *g, int const& i )
  {
    // Split PDF string
    const std::string pdfnames = g->getGenpdf();
    std::vector<std::string> pdfvec;
    std::stringstream ss(pdfnames); std::string s;
    while (getline(ss, s, ':')) pdfvec.push_back(s);

    const size_t isubproc = pdfvec.size() == 1 ? 0:i;
    return appl::appl_pdf::getpdf( pdfvec[isubproc] );
  }

  // Computes the normalisation required for translating APPLgrid weights to FK ones
  // e.g including factors of alpha_S, bin width etc.
  // g is the APPLgrid being combined with evolution factors.
  // pto specified the perturbative order being combined, as the value of alpha_S in the current bin,
  // and x1/x2 specify the numerical values of the PDF x-values for the first and second PDF respectively.
  double compute_wgt_norm( const appl::grid *g, int const& d, double const& pto, double const& as, double const& x1, double const& x2)
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

//    // Progress update ****************************************************

  int countElements(vector<int> const& maskmap, int const& min_pto, int const& max_pto, const appl::grid* g)
  {
    // Counter
    int nElm = 0;
    for (auto bin : maskmap )
      for (int pto=min_pto; pto<=max_pto; pto++)
      {
        const int gidx = get_grid_idx(g, pto); // APPLgrid grid index
        const appl::igrid *igrid = g->weightgrid(gidx, bin);
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
*/

 // *********************************** APPLGrid convolutions **************************************

  SubGrid::SubGrid(FKTarget const& parent, NNPDF::IndexDB const& db, int const& iDB):
    FKSubGrid(parent, iDB, NNPDF::dbquery<string>(db, iDB, "operators")),
    pineapplfile(NNPDF::dbquery<string>(db,iDB,"pineapplgrid")),
    readme(),
    maskmap(parse_maskmap(NNPDF::dbquery<string>(db,iDB,"mask"))),
    ndata(maskmap.size()),
    grid(pineappl_grid_read((pineapplPath() + parent.GetSetName() + "/" + pineapplfile).c_str()))
  {
    const size_t bin_count = pineappl_grid_bin_count(grid);
    if (ndata > bin_count)
    {
      cerr <<"Error: number of datapoints in settings: "<<ndata<<" > appl grid Nobs: "<<bin_count<<endl;
      cerr <<"Please check the provided mask" <<endl;
      exit(-1);
    }
  }

  void SubGrid::Compute(qcd_param const& par, vector<double>& results) const
  {
    if (par.evol_pto == 2)
      std::cout << "WARNING: APPLgrid does not currently support NNLO convolutions, fixing convolution to NLO" <<std::endl;

    // TODO: fix pto
    vector<double> xsec(pineappl_grid_bin_count(grid));
    pineappl_grid_convolute(grid, evolpdf_pineapplgrid, evolpdf_pineapplgrid, alphas,
                            nullptr, nullptr, nullptr, par.xiR, par.xiF, xsec.data());

    // apply database normalization rule
    for (double& obs : xsec)
      obs *= nrmdat;

    // Relate back to results
    for (size_t i=0; i<maskmap.size(); i++)
      for (int const& j : datamap[i])
        results[j] += xsec[maskmap[i]];
  }


//  *************************************** FK Table computation ********************************************

  void SubGrid::Combine(QCD::qcd_param const& par, NNPDF::FKGenerator* fk) const
  {
    if (par.evol_pto == 2)
      std::cout << "WARNING: APPLgrid does not currently support NNLO convolutions, fixing convolution to NLO" <<std::endl;

    // TODO: fit pto
    // Progress monitoring
    const vector<size_t> afl = QCD::active_flavours(par);

    vector<uint32_t> orders(pineappl_grid_order_count(grid));
    pineappl_grid_order_params(grid, orders.data());

    auto *grid_lumi = pineappl_grid_lumi(grid);
    const size_t lumi_count = pineappl_lumi_count(grid_lumi);

    vector<double> q2values(pineappl_subgrid_q2_grid_count(grid));
    pineappl_subgrid_q2_grid(grid, q2values.data());

    vector<double> xvalues(pineappl_subgrid_x_grid_count(grid));
    pineappl_subgrid_x_grid(grid, xvalues.data());
    const size_t nxpi = xvalues.size();

    vector<double> bin_sizes(pineappl_grid_bin_count(grid));
    pineappl_grid_bin_sizes(grid, bin_sizes.data());

    vector<double> weight_matrix(nxpi * nxpi);

    int completedElements = 0;
    int nXelements = 0;
    for (size_t d=0; d<maskmap.size(); d++)
    {
      const size_t bin = maskmap[d];
      for (size_t pto = 0; pto < orders.size()/4; pto++) // Loop over perturbative order
      {
        // TODO: unSkip renormalization and factorization scale variation.
        if (orders.at(4 * pto + 2) > 0 || orders.at(4 * pto + 3) > 0)
          continue;

        for (size_t lumi = 0; lumi < lumi_count; lumi++)
        {
          size_t slice_indices[2];
          pineappl_subgrid_filled_q2_slices(grid, pto, bin, lumi, slice_indices);
          nXelements += (slice_indices[1] - slice_indices[0]) * nxpi * nxpi;
        }
      }
    }

    const time_point t1 = std::chrono::system_clock::now();

    for (size_t d=0; d<maskmap.size(); d++)
    {
      // Fetch associated applgrid info
      const size_t bin = maskmap[d];

      for (size_t pto = 0; pto < orders.size()/4; pto++) // Loop over perturbative order
      {
        // TODO: unSkip renormalization and factorization scale variation.
        if (orders.at(4 * pto + 2) > 0 || orders.at(4 * pto + 3) > 0)
          continue;

        for (size_t lumi = 0; lumi < lumi_count; lumi++)
        {
          // prepare luminosity combinations
          const size_t combinations = pineappl_lumi_combinations(grid_lumi, lumi);
          vector<int32_t> pdgids(2 * combinations);
          vector<double> factors(combinations);
          pineappl_lumi_entry(grid_lumi, lumi, pdgids.data(), factors.data());

          // Fetch grid pointer and loop over Q
          size_t slice_indices[2];
          pineappl_subgrid_filled_q2_slices(grid, pto, bin, lumi, slice_indices);
          for (size_t t = slice_indices[0]; t < slice_indices[1]; t++)
          {
            // Scales and strong coupling
            const double Q2  = q2values.at(t);
            const double QF  = sqrt(Q2)*par.xiF;
            const double QR  = sqrt(Q2)*par.xiR;
            const double as  = QCD::alphas(QR);

            // Renormalisation and factorisation scale variation terms
            //const bool vary_ren = pto == 0 && par.evol_pto == 1 && par.xiR != 1.0;
            //const bool vary_fac = pto == 0 && par.evol_pto == 1 && par.xiF != 1.0;
            //const double renscale =  (as/(2.0*M_PI))*2.0*M_PI*QCD::beta0()*g->leadingOrder()*log(par.xiR*par.xiR);
            //const double facscale = -(as/(2.0*M_PI))*log(par.xiF*par.xiF);

            // define evolution factor arrays
            const int nxin = fk->GetNx();
            EvolutionFactors fA(nxin, nxpi);   // PDF 1 and 2 Evolution factors
            //EvolutionFactors fdA1(nxin, igrid->Ny1());  // PDF 1 Splitting function x Evolution factors
            //EvolutionFactors fdA2(nxin, igrid->Ny2());  // PDF 2 Splitting function x Evolution factors

            // Compute nonzero evolution factors
            for (int ix = 0; ix < nxin; ix++)
              for (size_t fl : afl)
              {
                for (size_t ox = 0; ox < nxpi; ox++) // Loop over applgrid x1
                {
                  QCD::EvolutionOperator(false, true, ix, xvalues.at(ox), fl, QF, fA(ox, ix, fl));
                  //if (vary_fac) QCD::DerivativeOperator(false,applgrid_nfl==14,ix,igrid->fx(igrid->gety1(ox)),fl,QF,fdA1(ox,ix,fl));
                }
              }

            // take the grid from pineappl
            pineappl_subgrid_q2_slice(grid, pto, bin, lumi, t, weight_matrix.data());

            for (size_t a = 0; a < nxpi; a++)
              for (size_t b = 0; b < nxpi; b++) // Loop over applgrid x2
              {
                // fetch weight values
                const double W = weight_matrix[a * nxpi + b];

                // Calculate normalisation factors
                const double norm = pow(as, orders.at(4 * pto + 0))/bin_sizes.at(bin)*nrmdat;

                for (int i = 0; i < nxin; i++)    // Loop over input pdf x1
                  for (int j = 0; j < nxin; j++)  // Loop over input pdf x2
                    for (size_t k : afl)         // loop over flavour 1
                      for (size_t l : afl)       // loop over flavour 2
                      {
                        double H = 0.0;
                        for (std::size_t c = 0; c != combinations; ++c)
                        {
                            auto const ai = pdgid_to_apfel_array_id(pdgids.at(2 * c + 0));
                            auto const bi = pdgid_to_apfel_array_id(pdgids.at(2 * c + 1));
                            H += factors.at(c) * (fA(a, i, k)[ai]) * (fA(b, j, l)[bi]);
                        }

                        // if (vary_fac)
                        // {
                        //   genpdf->evaluate(fdA1(a,i,k),fA2(b,j,l),H1);
                        //   genpdf->evaluate(fA1(a,i,k),fdA2(b,j,l),H2);
                        // }

                        double fill = H;                              // Basic fill
                        //if (vary_ren) fill += renscale*H[ip];             // Ren. scale variation
                        //if (vary_fac) fill += facscale*(H1[ip] + H2[ip]); // Fac. scale variation
                        if ( fill != 0.0 )
                          for (int const& td : datamap[d])
                            fk->Fill(td, i, j, k, l, norm*fill*W);
                    }

                // Update progress
                completedElements++;
                StatusUpdate(t1, (double)completedElements/(double)nXelements, cout);
              }
          }
        }
      } // /pto
    } // /data

    pineappl_lumi_delete(grid_lumi);

    return;
  }

  //  *************************************** Metadata methods ********************************************


  void SubGrid::Splash(ostream& o) const
  {
    FKSubGrid::Splash(o);
    o << "- PineAPPL: " << pineapplfile << endl;
  }

    // Get the maximum scale of an applgrid
  double SubGrid::GetQ2max() const
  {
    vector<double> q2values(pineappl_subgrid_q2_grid_count(grid));
    pineappl_subgrid_q2_grid(grid, q2values.data());
    return std::max(q2values.front(), q2values.back());
  }

  // Return the minimum x used in the subgrid
  double SubGrid::GetXmin() const
  {
    // TODO: check this implementation.
    return GetComputeXmin();
  };

  double SubGrid::GetComputeXmin() const
  {
    vector<double> xvalues(pineappl_subgrid_x_grid_count(grid));
    pineappl_subgrid_x_grid(grid, xvalues.data());
    return std::min(xvalues.front(), xvalues.back());
  };

  vector<int> SubGrid::parse_maskmap(string const& mask)
  {
    const vector<string> masksplit = ssplit(mask);
    vector<int> _maskmap;
    for (size_t i=0; i<masksplit.size(); i++)
      if ((bool) atoi(masksplit[i].c_str()))
        _maskmap.push_back(i);
    return _maskmap;
  }

}
