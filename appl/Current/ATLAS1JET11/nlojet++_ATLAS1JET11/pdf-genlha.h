#ifndef __pdf_genlha_h__
#define __pdf_genlha_h__ 1

#include <iostream>
#include <bits/hhc-process.h>
#include "LHAPDF.h"


//----- used namespaces -----
using namespace std;
using namespace nlo;



class pdf_genlha
  : public pdf_and_coupling_hhc
{
public:
  //   constructor
  inline explicit pdf_genlha(const string & name, unsigned int mem = 0)
    { 
      std::cout << "pdfset >>" << name << "<<" << std::endl; 
      initPDFset(name);
      initPDF(mem);
    }
  
  //   strong coupling
  double alpha_qcd(unsigned int nf, double mr2) {
    //    return alphasPDF(sqrt(mr2))/6.28318530717958647692;
    return alphasPDF(sqrt(mr2))/(2*M_PI);
  }
  
  //   the parton distribution function
  void hadronA(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13]; 
    evolvePDF(x, sqrt(Q2), __f);
    for(int i=-6; i <= 6; i++) f[i] = __f[6+i]/x;
  }
  
  void hadronB(double x, double Q2, unsigned int, unsigned int, double *f) {
    double __f[13];
    evolvePDF(x, sqrt(Q2), __f);
    for(int i=-6; i <= 6; i++) {
      f[i] = __f[6+i]/x; //p
      //f[i] = __f[6-i]/x; //pbar
    }
  }
};


class appl_pdf : public pdf_and_coupling<weight_hhc,2U,0U>
{
 public:  
  //   strong coupling
  inline double alpha_qcd(unsigned int, double) 
    {
      return 1.0;
    }
  
  
  inline weight_hhc pdf(double, double, double, unsigned int, unsigned int) 
    {
      weight_hhc retval;
      retval[0] = retval[1] = retval[2] = retval[3] = retval[4] = retval[5] = retval[6] = 1.0;
      return retval;
    }
  
};


#endif // __pdf_genlha_h__
