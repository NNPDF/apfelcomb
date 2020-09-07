#pragma once

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sys/stat.h>

// NNPDF
#include "NNPDF/fkgenerator.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/nnpdfdb.h"

// PineAPPL includes
#include "pineappl_capi.h"

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_grids.h"

using QCD::qcd_param;

namespace PINE
{
	class SubGrid: public FKSubGrid
	{
	public:
	void Compute(qcd_param const&, vector<double>&) 			const;	//!< Compute APPLgrid results mapped to Commondata
	void Combine(QCD::qcd_param const&, NNPDF::FKGenerator*) 	const;	//!< Perform the FK combination on a subgrid
	private:
		friend class ::FKTarget;
		SubGrid(FKTarget const& parent, NNPDF::IndexDB const& db, int const& iDB);

		void Splash(ostream&) 	const;					//!< Print metadata to stream
		size_t GetNdat() 		const {return ndata;};	//!< Return number of datapoints in subgrid
		double GetQ2max() 		const;					//!< Return maximum scale used in a subgrid
		double GetXmin() 		const;					//!< Return minimum x-value used in this sub grid
		double GetComputeXmin() const;					//!< Return minimal x-value used in computation of this observable

		// ***********************************************************

		const string 		pineapplfile; //!< Path for the pineapplg file
		const string 		readme;		  //!< APPLgrid README

		const vector<int> 	maskmap;	//!< Map of masked applgrid points to datapoints
		const size_t 		ndata;     	//!< Number of selected datapoints

		const pineappl_grid* grid;	//!< PineAPPL class

		static vector<int>  parse_maskmap(string const&);
	};

}

