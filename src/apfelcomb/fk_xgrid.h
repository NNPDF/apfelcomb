/*
 *  fk_xgrid.h
 *  PDF x grid distributions
 *  Code adapted from applgrid 1.4.70
 * *  nph 09/14
 */

#ifndef __nnpdfcpp__appl_xgrid__
#define __nnpdfcpp__appl_xgrid__

#include <stdio.h>
#include <cmath>

namespace XGrid
{
  // Generates an x-grid according to the APPLgrid prescription.
  // Points are distributed linearly in y(x) = ln(1/x) + a(1-x).
  // Parameter `a` controls point density at high-x.
  class Generator
  {
    public:
      Generator(const double a = 6):
      m_transvar(a){}

      double appl_fy(double x) const { return -std::log(x)+m_transvar*(1-x); }
      double appl_fx(double y) const {
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
        return std::exp(-yp);
    }

    private:
      const double m_transvar;

  };


}

#endif /* defined(__nnpdfcpp__appl_xgrid__) */
