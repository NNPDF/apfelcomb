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

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// APPLGrid main includes
#include "appl_grid/appl_grid.h"
#include "appl_grid/appl_igrid.h"
#include "appl_grid/appl_pdf.h"
#include "appl_grid/fastnlo.h"

#include "fk_utils.h"
#include "fk_qcd.h"

namespace APP
{
	// Grid data struct
	class appl_param: public QCD::qcd_param
	{
	public:
	  string setname;   //!< Name of the set to which the grid belongs
	  string gridname;  //!< Name of the grid to be generated
	  string gridfile;	//!< Path for the applgrid file

	  string desc;      //!< FKTable description
	  string readme;	//!< APPLgrid README
	  
	  size_t nbins;     //!< Number of bins in the grid
	  size_t ndata;     //!< Number of selected datapoints
	  size_t ptmin;     //!< Minimum perturbative order to contribute
	  
	  bool   fnlo;      //!< Using fastNLO grids
	  size_t fnlobin;   //!< Which bin in the fastNLO grid?
	  bool pdfwgt;      //!< Using a PDF weight? Typically used for fastNLO grids
	  bool ppbar;		//!< Does the grid need a pp -> ppbar transform?
	  vector<bool> mask;//!< Mask for which bins to pick out as datapoints
	  vector<int>  map; //!< Map to list which bins enter as datapoints
	  
	  int nx;           //!< Number of interpolation grid x-points to be used
	  double xmin;      //!< Minimum x-value to be used in interpolation

	  vector<string> inventory; //!< List of common grids required for set
	};

	// Parse APPLgrid data into FKHeader form
	void parse_input(int innum, appl_param& param);
  	void set_params(appl_param const& par, NNPDF::FKHeader& FK);

	// Populate FK table
	void computeFK(appl_param const&, appl::grid*, NNPDF::FKGenerator*);

	// Kinematic info fetchers
	double getXmin (const appl::grid* g, const bool& nonzero); // Gets smallest x-grid value - nonzero flag specifies whether to take absolute smallest value (for applgrid conv) or smallest nonzero (for FK conv)
	double getQ2max(const appl::grid* g);
}

