/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NONINFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and this permission notice shall 
    be included in all copies or substantial portions of YAX.
*/


#ifndef MATHSTATS_HPP
#define MATHSTATS_HPP

#include <vector>
#include <numeric>

#include "setOptions.hpp"

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/cauchy.hpp>

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <Spectra/SymEigsSolver.h>

// ------------------------------------
//  Eigen matrix printing formats
// ------------------------------------
const static Eigen::IOFormat EigenCSV(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
const static Eigen::IOFormat EigenTSV(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");

// ------------------------------------
//  R-like pdf and cdf (TODO: switch to libRmath)
// ------------------------------------
double qnorm(double, bool lower = false);
double pnorm(double, bool lower = false);
double qt(double, double, bool lower = false);
double pt(double, double, bool lower = false);
double qf(double, double, double, bool lower = false);
double pf(double, double, double, bool lower = false);
double qchisq(double, double, bool lower = false);
double pchisq(double, double, bool lower = false);
double qcauchy(double, bool lower = false);
double pcauchy(double, bool lower = false);

double ACAT(const std::vector<double>&);
double ACAT(const std::vector<double>&,const std::vector<double>&);

std::vector<double> filter_lt( const std::vector<double>&, double);
double ACAT_non_missing( const std::vector<double>&);

std::vector<int> rank_vector(const std::vector<double>&);

#endif

