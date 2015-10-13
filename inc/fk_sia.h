#pragma once
/*
 *  fk_sia.h
 *  APFEL SIA conversion to FK
 * *  sc
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
#include "fk_dis.h"

namespace SIA
{
	// Parse SIA data
  void parse_input(int innum, DIS::dis_param& param);
}

