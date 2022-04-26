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

#ifndef APEX_RMATHWRAPPERS_HPP
#define APEX_RMATHWRAPPERS_HPP

#include <limits>
#include <cmath>
#include "Rmath.h"

const long double LDBL_NEGINF = -std::numeric_limits<long double>::infinity();
const long double LDBL_POSINF =  std::numeric_limits<long double>::infinity();
const double DBL_NEGINF       = -std::numeric_limits<double>::infinity();
const double DBL_POSINF       =  std::numeric_limits<double>::infinity();
const long double LD_PI = acosl(-1.0L);
const long double LOG_MIN_LDBL = log(std::numeric_limits<long double>::min());
const double LOG_MIN_DBL = log(std::numeric_limits<double>::min());

namespace rmath {
  double pf(double x, double df1, double df2, int lower_tail, int log_p);
  double qcauchy(double p, double location, double scale, int lower_tail, int log_p);
  long double qcauchyl(long double p, long double location = 0, long double scale = 1, int lower_tail = 1, int log_p = 0);
  long double pcauchyl(long double x, double location = 0, double scale = 1, int lower_tail = 1, int log_p = 0);

  /**
   * Returns the standard cauchy(0,1) quantile function, but catches underflow if p is too small (or 0) and replaces
   * the return value with the smallest return value possible.
   * in long double precision.
   */
  long double bounded_stdqcauchy(long double p);

  /**
   * Returns standard cauchy(0, 1) CDF, also catching underflow to 0 and replacing with minimum long double value.
   * @param q
   * @return
   */
  long double bounded_stdpcauchy(long double q);

    /**
     * Returns log of F-distribution survival function (1 - CDF). If x grows too large,
     * this function will catch the underflow in the return value, and instead replace it
     * with the log of the minimum long double possible.
     */
  double bounded_log_pf(double x, double df1, double df2);

  /**
   * Returns exp(x), except if the value would cause exp to underflow, then
   * exp(log(minimum long double value)) is returned instead.
   */
  long double bounded_expl(long double x);
}
#endif //APEX_RMATHWRAPPERS_HPP
