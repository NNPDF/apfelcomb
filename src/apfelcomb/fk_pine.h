#pragma once

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <sys/stat.h>

// NNPDF
#include "NNPDF/fkgenerator.h"
#include "NNPDF/fastkernel.h"
#include "NNPDF/nnpdfdb.h"

// APPLGrid main includes
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/appl_pdf.h"
#include "appl_grid/fastnlo.h"

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_grids.h"

using QCD::qcd_param;

namespace PINE
{

	// // Wrapper for appl::grids (handles fastNLO based grids)
	class grid
	{
	public:
		grid(std::string const& filename, int const& fnlobin);
		~grid();

		fastnlo* 	fg;
		appl::grid* g;
	};


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

		const string 		applfile;	//!< Path for the applgrid file
		const string 		readme;		//!< APPLgrid README

		const vector<int> 	maskmap;	//!< Map of masked applgrid points to datapoints
		const size_t 		ndata;     	//!< Number of selected datapoints
		const size_t 		ptmin;     	//!< Minimum perturbative order to contribute

		const int    		fnlobin;   	//!< Which bin in the fastNLO grid? (yes, this needs to be signed)
		const bool 			pdfwgt;     //!< Using a PDF weight? Typically used for fastNLO grids
		const bool 			ppbar;		//!< Does the grid need a pp -> ppbar transform?

		const grid 			applgrid;	//!< APPLgrid class

		static vector<int>  parse_maskmap(string const&);
	};

}

