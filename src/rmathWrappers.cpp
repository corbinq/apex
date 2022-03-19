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

#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)
#define R_D__0	(log_p ? ML_NEGINF : 0.)
#define R_D__1	(log_p ? 0. : 1.)
#define R_D_Clog(p)	(log_p	? log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */
#define R_D_val(x)	(log_p	? log(x) : (x))		/*  x  in pF(x,..) */

namespace rmath {
  double pf(double x, double df1, double df2, int lower_tail, int log_p) {
    return ::pf(x, df1, df2, lower_tail, log_p);
  }

  double qcauchy(double p, double location, double scale, int lower_tail, int log_p) {
    return ::qcauchy(p, location, scale, lower_tail, log_p);
  }

  long double qcauchyl(long double p, long double location = 0, long double scale = 1, int lower_tail = 0, int log_p = 0) {
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

  long double pcauchyl(long double x, double location, double scale, int lower_tail, int log_p) {
    static long double ML_NEGINF = -std::numeric_limits<long double>::infinity();

    if (isnan(x) || isnan(location) || isnan(scale)) {
      return x + location + scale;
    }

    if (scale <= 0) {
      return NAN;
    }

    x = (x - location) / scale;

    if (isnan(x)) {
      return NAN;
    }

    if (!isfinite(x)) {
      if (x < 0) {
        return R_DT_0;
      }
      else {
        return R_DT_1;
      }
    }

    if (!lower_tail) {
      x = -x;
    }

    if (fabs(x) > 1) {
      long double y = atanl(1/x) / LD_PI;
      return (x > 0) ? R_D_Clog(y) : R_D_val(-y);
    }
    else {
      return R_D_val(0.5 + atanl(x) / LD_PI);
    }
  }
}