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

#include "rmathWrappers.hpp"

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

  long double bounded_stdqcauchy(long double p) {
    long double q = qcauchyl(p, 0, 1, 1, 0);
    if (q == DBL_NEGINF) {
      return -9.4675e+4930L;
    }
    return q;
  }

  long double bounded_stdpcauchy(long double q) {
    long double value = pcauchyl(q, 0, 1, 1, 0);
    if (value == 0.0L) {
      return std::numeric_limits<long double>::min();
    }
    return value;
  }

  double bounded_log_pf(double x, double df1, double df2) {
    double log_pf = ::pf(x, df1, df2, 0, 1);
    if (log_pf == DBL_NEGINF) {
      // If the value of x grows too large, in log space, pf() will return -inf.
      // Here we catch that and return log(minimum long double value) to represent the smallest
      // p-value possible.
      return LOG_MIN_DBL;
    }
    return log_pf;
  }

  long double bounded_expl(long double x) {
    long double value = expl(x);
    if (value == 0.0L) {
      return std::numeric_limits<long double>::min();
    }
    return value;
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