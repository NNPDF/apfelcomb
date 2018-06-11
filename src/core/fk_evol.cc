#include "apfelcomb/fk_evol.h"
#include <array>


namespace EVL
{
    EvolutionData::EvolutionData( dataInfoRaw const& di,
            vector<double> const& xgrid,
            vector<double> const& q2grid ):
        CommonData(di)
    {
        const int nx = xgrid.size();
        const int nq = q2grid.size();
        const int nf = 14;

        if (nx*nq*nf != fNData)
            throw std::runtime_error("Mismatch in nData for Evolution grid kinematics.");

        // Fill data kinematics
        int iDat = 0;
        for (auto x : xgrid)
            for (auto q2 : q2grid)
                for (int f=0; f < nf; f++)
                {
                    fData[iDat] = 0;
                    fStat[iDat] = 0;
                    fKin2[iDat] = x;
                    fKin2[iDat] = q2;
                    fKin3[iDat] = f;
                    iDat++;
                }
    }

    // allocate grid in x
    vector<double> generate_xgrid( const int nx,
            const double xmin,
            const double xmed,
            const double xmax)
    {
        vector<double> xg;
        const int nxm = nx/2;
        for (int ix = 1; ix <= nx; ix++)
        {
            if (ix <= nxm)
                xg.push_back(xmin*pow(xmed/xmin,2.0*(ix-1.0)/(nx-1.0)));
            else
                xg.push_back(xmed+(xmax-xmed)*((ix-nxm-1.0)/(nx-nxm-1.0)));
        }
        return xg;
    };

    // Generates a list of nq2 Q2 points between q2min and q2max.
    // Distributed as linear in tau = log( log( Q2 / lambda2 ) )
    vector<double> generate_q2grid( const int nq2,
            const double q2min,
            const double q2max,
            const double lambda2
            )
    {
        vector<double> q2grid;
        const double tau_min   = log( log( q2min / lambda2 ));
        const double tau_max   = log( log( q2max / lambda2 ));
        const double delta_tau = (tau_max - tau_min) /( (double) nq2 - 1);

        for (int iq2 = 0; iq2 < nq2; iq2++)
            q2grid.push_back(lambda2 * exp(exp( tau_min + iq2*delta_tau)));
        return q2grid;
    }

    // Compute Q2 grid
    // This is suitable for PDFs not AlphaS (see maxNF used below)
    vector< vector<double> > generate_q2subgrids( QCD::qcd_param par, const int nq2,
            const double q2min, const double q2max)
    {
        const double eps     = 1E-4;
        const double lambda2 = 0.0625e0;
        const int nfpdf = atoi(par.thMap["MaxNfPdf"].c_str());

        const std::array<double, 3> hqth = {
            pow(atof(par.thMap["Qmc"].c_str())*atof(par.thMap["kcThr"].c_str()),2),
            pow(atof(par.thMap["Qmb"].c_str())*atof(par.thMap["kbThr"].c_str()),2),
            pow(atof(par.thMap["Qmt"].c_str())*atof(par.thMap["ktThr"].c_str()),2),
        };

        // Determine subgrid edges
        vector<double> subgrid_edges={q2min};
        for (int i=0; i<(nfpdf-3); i++)
            if (hqth[i] > q2min) subgrid_edges.push_back(hqth[i]);
        subgrid_edges.push_back(q2max);

        // Determine point distribution across subgrids and generate subgrids
        const std::vector< double > point_dist = generate_q2grid(nq2, q2min, q2max, lambda2);
        std::vector< vector< double > > q2grid_subgrids;
        for (int i=0; i<subgrid_edges.size()-1; i++)
        {
            const double min_edge = subgrid_edges[i];
            const double max_edge = subgrid_edges[i+1];
            auto in_sub = [=](double q2pt){return (q2pt - min_edge) > -eps && (q2pt - max_edge) < eps;};
            const int npoints_sub = std::count_if(point_dist.begin(), point_dist.end(), in_sub);

            if (npoints_sub < 2)
                throw std::runtime_error("Too few points to sufficiently populate subgrids. More Q points required");

            const vector< double > subgrid = generate_q2grid(npoints_sub, min_edge, max_edge, lambda2);
            q2grid_subgrids.push_back(subgrid);
        }

        return q2grid_subgrids;
    }

}

