#include "apfelcomb/fk_evol.h"


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

}

