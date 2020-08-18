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

using namespace std;

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

double ACAT(const vector<double>&);
double ACAT(const vector<double>&,const vector<double>&);

vector<double> filter_lt( const vector<double>&, double);
double ACAT_non_missing( const vector<double>&);

vector<int> rank_vector(const vector<double>&);

#endif

