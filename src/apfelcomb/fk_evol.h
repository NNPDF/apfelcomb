#pragma once

#include <stdlib.h>
#include <fstream>
#include <vector>

// NNPDF
#include "NNPDF/fkgenerator.h"
#include "NNPDF/commondata.h"
#include "NNPDF/nnpdfdb.h"

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_grids.h"

using NNPDF::dataInfoRaw;
using NNPDF::CommonData;

namespace EVL
{
    // Class storing the required kinematic information for PDF evolution
    class EvolutionData: public CommonData
    {
        public:
        EvolutionData( dataInfoRaw const& di,
                       vector<double> const& xgrid,
                       vector<double> const& q2grid );

    };

    // Generate output (evolved scale) x-grid
    vector<double> generate_xgrid( const int nx,
            const double xmin,
            const double xmed,
            const double xmax);

    // Generate output q2 subgrids
    // One subgrid is generated for every HQ threshold crossed
    vector< vector<double> > generate_q2subgrids( QCD::qcd_param par, const int nq2,
            const double q2min, const double q2max);

}

