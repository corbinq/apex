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

namespace rmath {
  double pf(double x, double df1, double df2, int lower_tail, int log_p);
  double qcauchy(double p, double location, double scale, int lower_tail, int log_p);
  long double qcauchyl(long double p, long double location = 0, long double scale = 1, int lower_tail = 1, int log_p = 0);
  long double pcauchyl(long double x, double location = 0, double scale = 1, int lower_tail = 1, int log_p = 0);
}

#endif //APEX_RMATHWRAPPERS_HPP
