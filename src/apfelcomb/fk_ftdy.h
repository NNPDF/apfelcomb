#pragma once
/*
 *  fk_ftdy.h
 *  APFEL FTDY conversion to FK
 * *  nph 04/15
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

namespace FTDY
{
	// Populate FK header
	void set_params(QCD::qcd_param const& par, std::string const& gridname, std::string const& setname, int const& ndata, NNPDF::FKHeader& FK);

	// Populate FK table
	void computeGrid(QCD::qcd_param const&, NNPDF::CommonData const&);
  	void processFK(QCD::qcd_param const&, NNPDF::CommonData const&, std::string const&, NNPDF::FKGenerator*);

	// Kinematic info fetchers
	double getXmin (const NNPDF::CommonData& cdata);
	double getQ2max(const NNPDF::CommonData& cdata);
}

