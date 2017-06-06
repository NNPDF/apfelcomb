#pragma once
/*
 *  fk_appl.h
 *  APPLgrid conversion to FK
 * *  nph 09/14
 */

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

namespace APP
{

	// // Wrapper for appl::grids (handles fastNLO based grids)
	class grid
	{
	public:
		grid(std::string const& filename, int const& fnlobin):
		fg(fnlobin >= 0 ? new fastnlo(filename):0),
		g(fg ? fg->grids()[fnlobin]:new appl::grid(filename)){};
		~grid() {fg ? delete fg : delete g;};

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
		SubGrid(FKTarget const& parent, NNPDF::IndexDB const& db, int const& iDB):
		FKSubGrid(parent, iDB, NNPDF::dbquery<string>(db, iDB, "operators")),
		applfile(NNPDF::dbquery<string>(db,iDB,"applgrid")),
		readme(),
		maskmap(parse_maskmap(NNPDF::dbquery<string>(db,iDB,"mask"))),
		ndata(	maskmap.size() ),
		ptmin(	NNPDF::dbquery<size_t>(db,iDB,"ptmin")),
		fnlobin(NNPDF::dbquery<int>(db,iDB,"fnlobin")),
		pdfwgt(	NNPDF::dbquery<bool>(db,iDB,"pdfwgt")),
		ppbar(	NNPDF::dbquery<bool>(db,iDB,"ppbar")),
		applgrid(applPath() + parent.GetSetName() + "/" + applfile , fnlobin)
		{
		  if (ndata>static_cast<size_t>(applgrid.g->Nobs()))
		  {
		    cerr <<"Error: number of datapoints in settings: "<<ndata<<" > appl grid Nobs: "<<applgrid.g->Nobs()<<endl;
		    cerr <<"Please check the provided mask" <<endl;
		    exit(-1);
		  }
		};

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

