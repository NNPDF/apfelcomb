#pragma once
#include "NNPDF/pdfset.h"
#include "fk_qcd.h"

/*
 *  fk_pdf.h
 *  APFEL interface to NNPDF::PDFSet
 * *  nph 09/14
 */

  class APFELPDFSet : public NNPDF::PDFSet
  {
  public:
    //!< Constructor
    APFELPDFSet():
    PDFSet("ApfelSet",1, ER_MCT0)
    {

    };

    virtual ~APFELPDFSet() {};   //!< Destructor

    virtual void GetPDF(NNPDF::real const& x, NNPDF::real const& Q2, int const& n, NNPDF::real* pdf) const
    {
      double* LHA = new double[14];
    	QCD::evolpdf(x,sqrt(Q2),LHA);
      QCD::LHA2EVLN(LHA,pdf);
      delete[] LHA;
    	return;
    };

  private:
  };