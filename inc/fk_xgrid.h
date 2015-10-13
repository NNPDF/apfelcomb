/*
 *  fk_xgrid.h
 *  PDF x grid distributions
 * *  nph 09/14
 */

#ifndef __nnpdfcpp__appl_xgrid__
#define __nnpdfcpp__appl_xgrid__

#include <stdio.h>
#include <cmath>

namespace XGrid
{
  static double m_transvar = 6.0;

  static double appl_fy(double x) { return -std::log(x)+m_transvar*(1-x); }
  static double appl_fx(double y) {
    // use Newton-Raphson: y = ln(1/x)
    // solve   y - yp - a(1 - exp(-yp)) = 0
    // deriv:  - 1 -a exp(-yp)
    
    if ( m_transvar==0 )  return std::exp(-y);
    
    const double eps  = 1e-12;  // our accuracy goal
    const int    imax = 100;    // for safety (avoid infinite loops)
    
    double yp = y;
    double x, delta, deriv;
    for ( int iter=imax ; iter-- ; ) {
      x = std::exp(-yp);
      delta = y - yp - m_transvar*(1-x);
      if ( std::fabs(delta)<eps ) return x; // we have found good solution
      deriv = -1 - m_transvar*x;
      yp  -= delta/deriv;
    }
    // exceeded maximum iterations
    std::cerr << "_fx2() iteration limit reached y=" << y << std::endl;
    std::cout << "_fx2() iteration limit reached y=" << y << std::endl;
    return std::exp(-yp);
}

}

#endif /* defined(__nnpdfcpp__appl_xgrid__) */
