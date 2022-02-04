/*
    Copyright (C) 2020
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must
    be included in all copies or substantial portions of APEX.
*/

#include "Rmath.h"
#include <limits>

const long double LD_PI = acosl(-1.0L);

namespace rmath {
  double pf(double x, double df1, double df2, int lower_tail, int log_p) {
    return ::pf(x, df1, df2, lower_tail, log_p);
  }

  long double qcauchyl(long double p, long double location, long double scale, int lower_tail, int log_p) {
    if (isnan(p) || isnan(location) || isnan(scale)) {
      return p + location + scale;
    }

    if ((log_p && p > 0) || (!log_p && (p < 0 || p > 1))) {
      return NAN;
    }

    if (scale <= 0 || !isfinite(scale)) {
      if (scale == 0) return location;
      return NAN;
    }

    #define my_inf (location + (lower_tail ? scale : -scale) * std::numeric_limits<long double>::infinity())

    if (log_p) {
      if (p > -1) {
        if (p == 0.)  { /* needed, since 1/tan(-0) = -Inf  for some arch. */
          return my_inf;
        }

        lower_tail = !lower_tail;
        p = -expm1(p);
      }
      else {
        p = exp(p);
      }
    }
    else if (p == 1.) {
      return my_inf;
    }

    return location + (lower_tail ? -scale : scale) / tanl(LD_PI * p);
    /*	-1/tan(pi * p) = -cot(pi * p) = tan(pi * (p - 1/2))  */
  }
}