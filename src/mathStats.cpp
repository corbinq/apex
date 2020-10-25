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


#include "mathStats.hpp"


double p_bd = 1e-300;
double q_bd = 3e+299;

double qnorm(double p, bool lower){
	boost::math::normal N01(0.0, 1.0);
	if( lower ) return boost::math::quantile(boost::math::complement(N01, p));
	return boost::math::quantile(N01, p);
}

double pnorm(double x, bool lower){
	boost::math::normal N01(0.0, 1.0);
	if( lower ) return boost::math::cdf(boost::math::complement(N01, x));
	return boost::math::cdf(N01, x);
}

double qcauchy(double p, bool lower){
	p = p > p_bd ? p : p_bd;
	p = p < 1 - p_bd ? p : 1 - p_bd;
	
	boost::math::cauchy C01(0.0, 1.0);
	if( lower ) return boost::math::quantile(boost::math::complement(C01, p));
	return boost::math::quantile(C01, p);
}

double pcauchy(double x, bool lower){
	x = x < q_bd ? x : q_bd;
	x = x > -q_bd ? x : -q_bd;
	
	boost::math::cauchy C01(0.0, 1.0);
	if( lower ) return boost::math::cdf(boost::math::complement(C01, x));
	return boost::math::cdf(C01, x);
}

double qt(double p, double df, bool lower){
	boost::math::students_t TDIST(df);
	if( lower ) return boost::math::quantile(boost::math::complement(TDIST, p));
	return boost::math::quantile(TDIST, p);
}

double pt(double x, double df, bool lower){
	boost::math::students_t TDIST(df);
	if( lower ) return boost::math::cdf(boost::math::complement(TDIST, x));
	return boost::math::cdf(TDIST, x);
}

double qf(double p, double df1, double df2, bool lower){
	boost::math::fisher_f FDIST(df1, df2);
	if( lower ) return boost::math::quantile(boost::math::complement(FDIST, p));
	return boost::math::cdf(FDIST, p);
}

double pf(double x, double df1, double df2, bool lower){
	boost::math::fisher_f FDIST(df1, df2);
	if( lower ) return boost::math::cdf(boost::math::complement(FDIST, x));
	return boost::math::cdf(FDIST, x);
}

double qchisq(double p, double df, bool lower){
	boost::math::chi_squared CHISQ(df);
	if( lower ) return boost::math::quantile(boost::math::complement(CHISQ, p));
	return boost::math::quantile(CHISQ, p);
}

double pchisq(double x, double df, bool lower){
	boost::math::chi_squared CHISQ(df);
	if( lower ) return boost::math::cdf(boost::math::complement(CHISQ, x));
	return boost::math::cdf(CHISQ, x);
}

double ACAT(const std::vector<double>& pvals){
	long double sum_c = 0.0;
	double n = pvals.size();
	for( const double& p: pvals ){
		if( p >= 1 ){
			sum_c += (qcauchy(1 - 1/n, true)/n);
		}else if( p <= 0 ){
			std::cerr << "ACAT failed; input pval <= 0. \n";
			exit(1);
		}else{
			sum_c += (qcauchy(p, true)/n);
		}
	}
	return pcauchy(sum_c, true);
}

double ACAT(const std::vector<double>& pvals,const std::vector<double>& weights){
	long double sum_c = 0.0;
	long double denom = 0.0;
	double n = pvals.size();
	int i = 0;
	for( const double& w: weights ){
		denom += weights[i];
	}
	if( denom <= 0){
		return -99;
	}
	for( const double& p: pvals ){
		if( p >= 1 ){
			sum_c += (weights[i] * qcauchy(1 - 1/n, true) / denom);
		}else if( p <= 0 ){
			std::cerr << "ACAT failed; input pval <= 0. \n";
			exit(1);
		}else{
			sum_c += (weights[i] * qcauchy(p, true) / denom);
		}
		i++;
	}
	return pcauchy(sum_c, true);
}

std::vector<double> filter_lt( const std::vector<double>& p, double thresh){
	std::vector<double> out;
	for( const double& x : p ){
		if( x > thresh ) out.push_back(x);
	}
	return out;
}

double ACAT_non_missing( const std::vector<double>& pvals){
	long double sum_c = 0.0;
	double n = 0;
	double n_p1 = 0;
	for( const double& p: pvals ){
		if( !std::isnan(p) ){
			if( p >= 1 ){
				n_p1 += 1;
				n += 1;
			}else if( p > 0 ){
				sum_c += qcauchy(p, true);
				n += 1;
			}
		}
	}
	if( n_p1 > 0 ){
		if( n < 4 ){
			n = 4;
		}
		sum_c += n_p1 * qcauchy(1 - 1/n, true);
	}else if(n < 1){
		n = 1;
	}
	return pcauchy(sum_c/n, true);
}


std::vector<int> rank_vector(const std::vector<double>& v)
{
	// This code is adapted from stackoverflow:
	// https://stackoverflow.com/questions/30822729/create-ranking-for-vector-of-double
	
    std::vector<size_t> w(v.size());
    iota(begin(w), end(w), 0);
    sort(begin(w), end(w), 
        [&v](size_t i, size_t j) { return v[i] < v[j]; });

    std::vector<int> r(w.size());
    for (size_t n, i = 0; i < w.size(); i += n)
    {
        n = 1;
        while (i + n < w.size() && v[w[i]] == v[w[i+n]]) ++n;
        for (size_t k = 0; k < n; ++k)
        {
            r[w[i+k]] = i + (n + 1) / 2.0; // average rank of n tied values
            // r[w[i+k]] = i + 1;          // min 
            // r[w[i+k]] = i + n;          // max
            // r[w[i+k]] = i + k + 1;      // random order
        }
    }
    return r;
}
