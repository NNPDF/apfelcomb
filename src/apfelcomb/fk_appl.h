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
#include "NNPDF/commondata.h"
#include "NNPDF/nnpdfdb.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// APPLGrid main includes
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/appl_pdf.h"
#include "appl_grid/fastnlo.h"

#include "fk_utils.h"
#include "fk_qcd.h"
#include "fk_grids.h"

namespace APP
{

	class SubGrid: public FKSubGrid
	{
	public:
		SubGrid(NNPDF::IndexDB const& db, int const& iDB):
		FKSubGrid( iDB, NNPDF::dbquery<string>(db, iDB, "operators")),
		applgrid(NNPDF::dbquery<string>(db,iDB,"applgrid")),
		readme(),
		maskmap(parse_maskmap(NNPDF::dbquery<string>(db,iDB,"mask"))),
		ndata(	maskmap.size() ),
		ptmin(	NNPDF::dbquery<size_t>(db,iDB,"ptmin")),
		fnlobin(NNPDF::dbquery<int>(db,iDB,"fnlobin")),
		pdfwgt(	NNPDF::dbquery<bool>(db,iDB,"pdfwgt")),
		ppbar(	NNPDF::dbquery<bool>(db,iDB,"ppbar"))
		{

		};

		size_t GetNdat() {return ndata;};	//!< Return number of datapoints in subgrid
		double GetXmin();					//!< Return minimum x-value used in this sub grid
		double GetComputeXmin();			//!< Return minimal x-value used in computation of this observable

		void Splash(ostream&);

	private:
	  const string 			applgrid;	//!< Path for the applgrid file
	  const string 			readme;		//!< APPLgrid README
	  
	  const vector<int> 	maskmap;	//!< Map of masked applgrid points to datapoints
	  const size_t 			ndata;     	//!< Number of selected datapoints
	  const size_t 			ptmin;     	//!< Minimum perturbative order to contribute
	  
	  const int    			fnlobin;   	//!< Which bin in the fastNLO grid? (yes, this needs to be signed)
	  const bool 			pdfwgt;     //!< Using a PDF weight? Typically used for fastNLO grids
	  const bool 			ppbar;		//!< Does the grid need a pp -> ppbar transform?

	  static vector<int>  parse_maskmap(string const&);

	};


	// Grid data struct
	class appl_param: public QCD::qcd_param
	{
	public:
	  string setname;   //!< Name of the set to which the grid belongs
	  string gridname;  //!< Name of the grid to be generated
	  string applgrid;	//!< Path for the applgrid file

	  string desc;      //!< FKTable description
	  string readme;	//!< APPLgrid README
	  
	  size_t ndata;     //!< Number of selected datapoints
	  size_t ptmin;     //!< Minimum perturbative order to contribute
	  
	  int    fnlobin;   //!< Which bin in the fastNLO grid? (yes, this needs to be signed)
	  bool pdfwgt;      //!< Using a PDF weight? Typically used for fastNLO grids
	  bool ppbar;		//!< Does the grid need a pp -> ppbar transform?

	  vector<bool> mask;				//!< Mask for which bins to pick out as datapoints
	  vector<int>  maskmap; 			//!< Map to list which bins enter as datapoints
	  vector< vector<int> > datamap;	//!< Map of applgrid point to data point

	  size_t nx;           //!< Number of interpolation grid x-points to be used
	  double xmin;      //!< Minimum x-value to be used in interpolation

	  size_t incdat;    //!< Datapoint increment operator
	  size_t muldat;	//!< Datapoint multiplicative operator
	  double nrmdat;	//!< Datapoint normalisation operator
	  
	  vector<int> 	 common_subgrids;	//!< IDs of common subgrids required for grid
	};

	// Wrapper for appl::grids (handles fastNLO based grids)
	class grid
	{
	public:
		grid(appl_param const& par):
		fg(par.fnlobin >= 0 ? new fastnlo(par.applgrid):0),
		g(fg ? fg->grids()[par.fnlobin]:new appl::grid(par.applgrid)){};
		~grid() {fg ? delete fg : delete g;};

		fastnlo* 	fg;
		appl::grid* g;
	};

	// Parse APPLgrid data into FKHeader form
	void parse_input(int innum, appl_param& param, bool silent = false);
	double parse_xmin(std::vector<int> const&);
  	std::vector< std::vector<int> > parse_map(std::vector<int> const&, int const&);

  	void set_params(appl_param const& par, NNPDF::FKHeader& FK);
  	double computeTargetPrecision(std::vector<int> const& targetPoints, NNPDF::CommonData const&);

	// Populate FK table
	void computeFK(appl_param const&, appl::grid*, NNPDF::FKGenerator*);

	// Kinematic info fetchers
	double getXmin (const appl::grid* g, const bool& nonzero); // Gets smallest x-grid value - nonzero flag specifies whether to take absolute smallest value (for applgrid conv) or smallest nonzero (for FK conv)
	double getQ2max(const appl::grid* g);
}

