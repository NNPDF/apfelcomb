#pragma once
/*
 *  fk_sia.h
 *  APFEL SIA conversion to FK
 * *  sc
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

namespace SIA
{
	// A SIA SubGrid essentially just aliases the parent FKSubGrid's CommonData attribute
	// I put this into a 'SubGrid' such that we can use the same testing/optimisation procedures as in APPLgrids
	class SubGrid: public FKSubGrid
	{
	public:
	  void Compute(QCD::qcd_param const&, vector<double>&) const;	       	//!< Compute APPLgrid results mapped to Commondata
	  void Combine(QCD::qcd_param const&, NNPDF::FKGenerator*) const;	//!< Perform the FK combination on a subgrid
	private:
	  friend class ::FKTarget;
	SubGrid(FKTarget const& parent, NNPDF::IndexDB const& db, int const& iDB):
	  FKSubGrid(parent, iDB, NNPDF::dbquery<string>(db, iDB, "operators")),
	    process(NNPDF::dbquery<string>(db,iDB,"process"))
	      {
		if (incdat != 0 || muldat != 1)
		  {
		    cerr << "Error: SIA grids do not support operators other than normalisation" << endl;
		    exit(-1);
		  }
	      };
	  
	  void Splash(ostream&) 	const;	//!< Print metadata to stream
	  size_t GetNdat() 		const;	//!< Return number of datapoints in subgrid
	  double GetQ2max() 		const;	//!< Return maximum scale used in a subgrid
	  double GetXmin() 		const;	//!< Return minimum x-value used in this sub grid
	  
	  // **********************************************************
	  
	  const string process;			//!< Process string of the observable
	};
}

/*
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sys/stat.h>

// NNPDF
#include "NNPDF/fkgenerator.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/commondata.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_dis.h"

namespace SIA
{
  // Parse SIA data
  void parse_input(int innum, DIS::dis_param& param);
}
*/
