#pragma once
/*
 *  fk_ftdy.h
 *  APFEL FTDY conversion to FK
 * *  nph 04/15
 */

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

namespace FTDY
{
  // A FTDY SubGrid handles the fixed-target drell-yan chaos.
  class SubGrid: public FKSubGrid
  {
  public:
  void Compute(QCD::qcd_param const&, vector<double>&) const;     //!< Compute APPLgrid results mapped to Commondata
  void Combine(QCD::qcd_param const&, NNPDF::FKGenerator*) const;   //!< Perform the FK combination on a subgrid
  private:
    friend class ::FKTarget;
    SubGrid(FKTarget const& parent, NNPDF::IndexDB const& db, int const& iDB):
    FKSubGrid(parent, iDB, NNPDF::dbquery<string>(db, iDB, "operators"))
    {
      if (incdat != 0 || muldat != 1 || nrmdat != 1.0)
      {
        cerr << "Error: DYP grids do not support operators" << endl;
        exit(-1);
      }
    };

    size_t 	GetNdat() 			const; 	//!< Number of datapoints 
    void 	Splash(ostream&)   	const;  //!< Print metadata to stream
    double 	GetQ2max()     		const;  //!< Return maximum scale used in a subgrid
    double 	GetXmin()    		const;  //!< Return minimum x-value used in this sub grid

    // **********************************************************

  };
}

