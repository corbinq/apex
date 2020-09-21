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

