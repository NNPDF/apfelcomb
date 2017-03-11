#pragma once
/*
 *  fk_dis.h
 *  APFEL DIS conversion to FK
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

// LHAPDF
#include "LHAPDF/LHAPDF.h"

#include "fk_utils.h"
#include "fk_qcd.h"

namespace DIS
{
	// Grid data struct
	class dis_param: public QCD::qcd_param
	{
	public:
	  std::string setname;   	//!< Name of the set to which the grid belongs
	  std::string gridname;  	//!< Name of the grid to be generated
  	  std::string process;  	//!< Process of the grid

	  std::string commonfile;	//!< Path for the commondata file
	  std::string sysfile;		//!< Path for the SYSTYPE file
	  	  
	  std::string desc;      //!< FKTable description
	  
	  size_t ndata;     //!< Number of selected datapoints

	  bool positivity;	//!< Is a positivity observable
	  	  
	  int nx;           //!< Number of interpolation grid x-points to be used
	  double xmin;      //!< Minimum x-value to be used in interpolation
	  double maxprec;   //!< Highest experimental precision reached in dataset
	};

	// Parse DIS data
	void parse_input(int innum, dis_param& param, std::string dbfile);
  	void set_params(dis_param const& par, NNPDF::FKHeader& hd);

	// Populate FK table
	void computeFK(dis_param const&, NNPDF::CommonData const&, NNPDF::FKGenerator*);

	// Kinematic info fetchers
	double getXmin (const NNPDF::CommonData& cdata);
	double getQ2max(const NNPDF::CommonData& cdata);
	double getQ2min(const NNPDF::CommonData& cdata);
}

