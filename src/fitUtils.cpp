/*  fitUtils:

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

#include "fitUtils.hpp"


Eigen::MatrixXd resid( const Eigen::MatrixXd& Y, const Eigen::MatrixXd& X ){
	return Y - X * (X.transpose() * X).ldlt().solve(X.transpose() * Y);
}

void rank_normalize(Eigen::MatrixXd& Y){
	double m = Y.cols();
	double n = Y.rows();
	
	std::vector<double> z((int) n);
	std::vector<double> rk((int) n);
	
	double mu = 0;
	double sd = 0;
	for(int i = 0; i < n; ++i){
		z[i] = qnorm( ((double)i+1.0)/((double)n+1.0) );
		mu += z[i];
		sd += z[i]*z[i];
	}
	sd = std::sqrt( sd/(n - 1) - mu*mu/( n*(n - 1.0) ) );
	mu = mu/n;
	for(int i = 0; i < n; ++i){
		z[i] = (z[i] - mu)/sd;
	}
	
	for(int j = 0; j < m; ++j){

		std::vector<double> v(n);
		for(int i = 0; i< n; ++i){
			v[i] = Y(i,j);
		}
		
		std::vector<int> ranks = rank_vector(v);
		
		for(int i = 0; i< n; ++i){
			 Y(i,j) = z[ ranks[i]-1 ];
		}
	}
}

void scale_and_center(Eigen::MatrixXd& Y, std::vector<double>& sd_vec){
	double m = Y.cols();
	double n = Y.rows();
	
	if( sd_vec.size() != m ) sd_vec = std::vector<double>(m, 0.0);
	
	for(int j = 0; j < m; ++j){
		double mu = 0;
		double sd = 0;
		for(int i = 0; i < n; ++i){
			mu += Y(i,j);
			sd += Y(i,j)*Y(i,j);
		}
		sd = std::sqrt( sd/(n - 1) - mu*mu/( n*(n - 1.0) ) );
		sd_vec[j] = sd;
		mu = mu/n;
		for(int i = 0; i < n; ++i){
			Y(i,j) = (Y(i,j) - mu)/sd;
		}
	}
}

void scale_and_center(Eigen::MatrixXd& Y){
	double m = Y.cols();
	double n = Y.rows();
	
	for(int j = 0; j < m; j++){
		double mu = 0;
		double sd = 0;
		for(int i = 0; i < n; i++){
			mu += Y(i,j);
			sd += Y(i,j)*Y(i,j);
		}
		sd = std::sqrt( sd/(n - 1) - mu*mu/( n*(n - 1.0) ) );
		mu = mu/n;
		for(int i = 0; i < n; i++){
			Y(i,j) = (Y(i,j) - mu)/sd;
		}
	}
}

Eigen::MatrixXd get_half_hat_matrix(const Eigen::MatrixXd& X){
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> XtX_es(X.transpose() * X);
	Eigen::VectorXd lambda = XtX_es.eigenvalues();
	for( auto& a : lambda ){
		a = 1/std::sqrt(a);
	}
	Eigen::MatrixXd U = X * XtX_es.eigenvectors() * lambda.asDiagonal();
	return U;
}

void make_half_hat_matrix(Eigen::MatrixXd& X){
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> XtX_es(X.transpose() * X);
	Eigen::VectorXd lambda = XtX_es.eigenvalues();
	for( auto& a : lambda ){
		a = 1/std::sqrt(a);
	}
	X *= XtX_es.eigenvectors();
	X *= lambda.asDiagonal();
	return;
}

Eigen::MatrixXd resid_from_half_hat( const Eigen::MatrixXd& Y, const Eigen::MatrixXd& C ){
	Eigen::MatrixXd CtY = C.transpose() * Y;
	return (Y - C * CtY).eval();
}

void make_resid_from_half_hat( Eigen::MatrixXd& Y, const Eigen::MatrixXd& C ){
	Eigen::MatrixXd CtY = C.transpose() * Y;
	Y.noalias() -= (C *CtY);
	return;
}

Eigen::VectorXd resid_vec_from_half_hat( const Eigen::VectorXd& Y, const Eigen::MatrixXd& C ){
	Eigen::VectorXd CtY = C.transpose() * Y;
	Eigen::VectorXd Yhat = C * CtY;
	return (Y - Yhat).eval();
}

void appendInterceptColumn( Eigen::MatrixXd &X ){
	X.conservativeResize(Eigen::NoChange, X.cols()+1);
	X.col(X.cols()-1) = Eigen::VectorXd::Ones(X.rows());
	return;
}
