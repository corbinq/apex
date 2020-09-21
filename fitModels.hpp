/*
	Corbin's notes on fitModels:
	
	Model fitting routines ...

*/

#ifndef FITMODELS_HPP
#define FITMODELS_HPP

#include <math.h>
#include <vector>
#include <functional>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <brent_fmin.hpp>

#include "setOptions.hpp"
#include "metaAnalysis.hpp"
#include "miscUtils.hpp"
#include "mathStats.hpp"


double get_neg_logLik_REML(const double& delta, Eigen::MatrixXd& X_tilde, Eigen::VectorXd& y_tilde, Eigen::VectorXd& lambda );

class LMM_fitter{
	
	private:
		const Eigen::Ref<Eigen::MatrixXd> X_tilde;
		const Eigen::Ref<Eigen::VectorXd> y_tilde;
		const Eigen::Ref<Eigen::VectorXd> lambda;
		double df_resid;
		
	public:
	
		Eigen::DiagonalMatrix<double,Eigen::Dynamic> Vi;
		double sigma2;
		double phi;
	
		LMM_fitter(const Eigen::Ref<Eigen::MatrixXd> X_, const Eigen::Ref<Eigen::VectorXd> y_, const Eigen::Ref<Eigen::VectorXd> lam) : 
			X_tilde(X_),
			y_tilde(y_),
			lambda(lam)
		{df_resid = X_tilde.rows() - X_tilde.cols();};

		double neg_logLik_REML(double phi_)
		{
			Eigen::VectorXd vals = 1.00 + phi_*lambda.array();
			
			Eigen::DiagonalMatrix<double,Eigen::Dynamic> Di = vals.asDiagonal().inverse();
			
			Eigen::MatrixXd XtDX = X_tilde.transpose() * Di * X_tilde;
			Eigen::VectorXd XtDy = X_tilde.transpose() * Di * y_tilde;
			Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
			
			double sigma2_ = (y_tilde.dot(Di * y_tilde) - XtDy.dot(b))/df_resid;
			
			double ll = 0.5*(
				df_resid*std::log(sigma2_) + 
				vals.array().log().sum() + std::log(XtDX.determinant()) + 1.0
			);
			
			return ll;
		};
	
		void fit_REML(){
			
			std::function<double(double)> f = [this](double x){return neg_logLik_REML(x);};
			phi = Brent_fmin(0.00, 5000, f, 2e-5);
			Eigen::VectorXd vals = 1.00 + phi*lambda.array();
			
			Vi = vals.asDiagonal().inverse();
			
			Eigen::MatrixXd XtDX = X_tilde.transpose() * Vi * X_tilde;
			Eigen::VectorXd XtDy = X_tilde.transpose() * Vi * y_tilde;
			Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
			
			sigma2 = (y_tilde.dot(Vi * y_tilde) - XtDy.dot(b))/df_resid;
			
			return;
		}
		
};
	
class vcov_getter;

static const Eigen::VectorXd vec0 = Eigen::VectorXd(0);
static const Eigen::MatrixXd mat0 = Eigen::MatrixXd(0,0);

class ss_lm_single
{
	public:
		double U;
		double V;
		double SSR;

		double beta;
		double se;
		double df;
		double pval;
		
		ss_lm_single () {};
		
		ss_lm_single(const double& U_, const double& V_, const double& df_, const double& SSR_, const double& scale_ = 1.0){
			U = U_;
			V = V_;
			SSR = SSR_;
			df = df_;
			beta = U_/V_;
			se = std::sqrt((SSR_/V_ - beta*beta)/( df - 1 ));
			beta *= scale_;
			se *= scale_;
			pval = -99;
			if( se > 0 ){
				double tstat = beta/se;
				if( !std::isnan(tstat) && tstat*tstat > 0 && df > 2 ){
					pval = pf( tstat*tstat, 1, df - 1, true );
				}
			}
		};
		
		ss_lm_single(const std::vector<ss_lm_single>& ss_per_study){
			beta = 0; se = 0; df = 0; 
			U = 0; V = 0; SSR = 0;
			pval = 1.0;
			double denom = 0.0, denom_beta = 0.0;
			for( const ss_lm_single& i : ss_per_study ){
				if( i.se > 0 && !std::isnan(i.beta) ){
					double w_beta = 1/(i.se*i.se);
					beta += w_beta*i.beta;
					df += i.df;
					denom_beta += w_beta;
					
					double w = w_beta/i.V;
					U += w*i.U;
					V += w*i.V;
					SSR += w*i.SSR;
					denom += denom;
				}
			}
			beta /= denom_beta;
			se = std::sqrt(1/denom_beta);
			U /= denom;
			V /= denom;
			SSR /= denom;
			pval = -99;
			if( se > 0 ){
				double tstat = beta/se;
				if( !std::isnan(tstat) && tstat*tstat > 0 ){
					pval = pf( tstat*tstat, 1, df, true );
				}
			}
		}
};


class svar_sumstat
{
	public:
	Eigen::VectorXd U;
	Eigen::VectorXd V;
	
	double DF;
	double SSR;
	
	double SCALE_Y_SD;
	
	std::vector<int> c;
	
	Eigen::MatrixXd V_c;
	Eigen::MatrixXd C_c;
	
	
	Eigen::VectorXd U_0;
	Eigen::VectorXd V_0;
	
	double DF_0;
	double SSR_0;
	
	svar_sumstat(){};
	svar_sumstat(const Eigen::VectorXd& U_, const Eigen::VectorXd& V_, const double& DF_, const double& SSR_, const double& SCALE_Y_SD_ = 1.0 ){
		U = U_; V = V_; 
		U_0 = U; V_0 = V;
		DF = DF_; SSR = SSR_;
		DF_0 = DF; SSR_0 = SSR;
		SCALE_Y_SD = SCALE_Y_SD_;
		// std::cout << DF << "\t" << SSR << "\n";
		return;
	};
	
	void condition_on( const std::vector<int>& k, const Eigen::MatrixXd& C_k ){
		
		int nk = k.size();
		int nc = c.size();
		if( C_k.rows() != U.size() || C_k.cols() != nk ){
			std::cerr << "\n\nFatal covar dim: " << C_k.rows() << "," << C_k.cols() << " != " << U.size() << "," << C_k.cols() << "\n\n";
			abort();
		}
		
		Eigen::MatrixXd V_k = Eigen::MatrixXd::Zero(nk, nk);
		Eigen::MatrixXd U_k = Eigen::VectorXd::Zero(nk);
		Eigen::MatrixXd V_ck = Eigen::MatrixXd::Zero(nc, nk);
		for(int i = 0; i < nk; i++){
			V_k(i, i) = V(k[i]);
			for( int j = 0; j < c.size(); j++){
				// omit this later
				if( std::abs(C_k(c[j], i) - C_c(k[i], j)) > 0.0001 ){
					std::cerr << "\n\nFatal covar mismatch: " << C_k(c[j], i) << ", " << C_c(k[i], j) << "\n\n";
					abort();
				}
				V_ck(j, i) = C_k(c[j], i);
			}
			for( int j = i; j < nk; j++){
				// omit this later
				if( std::abs(C_k(k[j], i) - C_k(k[i], j)) > 0.0001 ){
					std::cerr << "\n\nFatal covar mismatch: " << C_k(k[j], i) << ", " << C_k(k[i], j) << "\n\n";
					abort();
				}
				V_k(i, j) = C_k(k[j], i);
				V_k(j, i) = V_k(i, j);
			}
		}
		
		if( c.size() == 0 ){
			V_c = V_k;
			C_c = C_k;
		}else{
			V_c.conservativeResize(nc + nk, nc + nk);
			V_c.block(nc,nc,nk,nk) = V_k;
			V_c.block(nc,0,nk,nc) = V_ck.transpose();
			V_c.block(0,nc,nc,nk) = V_ck;
			C_c.conservativeResize(Eigen::NoChange, nc + nk);
			C_c.rightCols(nk) = C_k;
		}
		for( const double& kk : k ){
			c.push_back(kk);
		}
		Eigen::VectorXd U_c = U_0(c);
		Eigen::MatrixXd V_c_i = V_c.inverse();
		
		DF = DF_0 - nc - nk;
		SSR = SSR_0 - (double) U_c.dot(V_c_i * U_c);
			
		Eigen::MatrixXd CVc_i = C_c * V_c_i;
			
		U = U_0 - CVc_i * U_c;
		
		// std::cout << U(k[0]) << "\n";
		
		V = V_0 - C_c.cwiseProduct(CVc_i).rowwise().sum();
		
		return;
	};
	
	void condition_on_subset( const std::vector<int>& k ){
		
		int nk = k.size();
		
		Eigen::MatrixXd V_k = V_c(k, k);
		Eigen::VectorXd U_k = (U_0(c))(k);
		Eigen::MatrixXd C_k = C_c(Eigen::all, k);
		
		Eigen::MatrixXd V_k_i = V_k.inverse();
		DF = DF_0 - nk;
		SSR = SSR_0 - (double) U_k.dot(V_k_i * U_k);
			
		Eigen::MatrixXd CVk_i = C_k * V_k_i;
			
		U = U_0 - CVk_i * U_k;
		V = V_0 - C_k.cwiseProduct(CVk_i).rowwise().sum();
		
		return;
	};
	
	void drop_snp( int i ){
		
		std::vector<int> k = seq_int(c.size());
		k.erase(k.begin() + i);
		
		condition_on_subset(k);
		
		V_c = (V_c(k, k)).eval();
		C_c = (C_c(Eigen::all, k)).eval();
		c.erase(c.begin() + i);
	}
	
	ss_lm_single final_model(const int& c_i){
		
		if( c_i > c.size() ){
			std::cerr << "Attempted SNP #" << c_i << " in a " << c.size() << " SNP model.\n";
			abort();
		}
		int ii = c[c_i];
		if( c.size() == 1 ){	
			return ss_lm_single( U_0(ii), V_0(ii), DF_0, SSR_0 );
		}
		
		std::vector<int> cc = c;
		std::vector<int> cc_i = seq_int(c.size());
		
		cc.erase (cc.begin()+c_i);
		cc_i.erase (cc_i.begin()+c_i);
		
		// std::cout << V_c << "\n\n";
		
		Eigen::MatrixXd V_cc_inv = V_c(cc_i,cc_i).inverse();
		
		// std::cout << V_cc_inv << "\n\n";
		
		Eigen::VectorXd C_c_i = (C_c.col(c_i))(cc);
		
		// std::cout << C_c_i << "\n\n";
		
		return ss_lm_single(
			U_0(ii) - C_c_i.dot(V_cc_inv* U_0(cc)),
			V_0(ii) - C_c_i.dot(V_cc_inv*C_c_i), 
			DF_0 - cc.size(), 
			SSR_0 - U_0(cc).dot(V_cc_inv*U_0(cc))
		);
	}
	
	ss_lm_single marginal_model(const int& c_i){
		
		if( c_i > c.size() ){
			std::cerr << "Attempted SNP #" << c_i << " in a " << c.size() << " SNP model.\n";
			abort();
		}
		int ii = c[c_i];
		return ss_lm_single( U_0(ii), V_0(ii), DF_0, SSR_0 );
	}
	
	
	void min_pval(int& k, double& min_p){
		// Eigen::VectorXd P_VAR = U.cwiseProduct(U).cwiseQuotient(V);
		double MAX_P_VAR = -1;
		k = -1;
		min_p = -99;
		for(int i = 0; i < U.size(); i++){
			if( V(i) > 0 && V_0(i) > 0 ){
				if( V(i)/V_0(i) > 1.0 - global_opts::RSQ_PRUNE ){
					if( U(i)*U(i)/V(i) > MAX_P_VAR ){
						MAX_P_VAR = U(i)*U(i)/V(i);
						k = i;
					}
				}
			}
		}
		double STAT = MAX_P_VAR*DF/( SSR - MAX_P_VAR );
		if( STAT > 0 ){
			min_p = pf(STAT, 1.0, DF - 1, true );
		}
	}
	
	double single_snp_pval(const int& i){
		double P_VAR = U(i)*U(i)/V(i);
		double STAT = P_VAR*DF/( SSR - P_VAR );
		if( SSR > P_VAR && !std::isnan(STAT) ){
			return pf(STAT, 1.0, DF - 1, true );
		}else{
			return (double) -99.0;
		}
	}
	
	void acat_min_pval(int& k, double& min_p, double& acat_p){
		Eigen::VectorXd P_VAR = U.cwiseProduct(U).cwiseQuotient(V);
		
		double MAX_P_VAR = -1.0, NUMER = 0.0, DENOM = 0.0, N_PVAL1 = 0.0;
		
		k = -1;
		min_p = -99;
		for(int i = 0; i < U.size(); i++){
			if( V(i) > 0 && V_0(i) > 0 ){
				if( V(i)/V_0(i) > 1.0 - global_opts::RSQ_PRUNE ){
					
					DENOM += 1;
					double PVAL = 1.0;
					
					double STAT = 0;
					
					if( P_VAR(i) > 0 && !std::isnan( P_VAR(i) ) ){
						STAT = P_VAR(i)*DF/( SSR - P_VAR(i) );
						
						if( STAT > 0 && !std::isnan(STAT) ){
							PVAL = pf( STAT, 1.0, DF - 1, true );
						}
					}
					/*else{
						std::cerr << PVAL << ", " << P_VAR(i) << ", " << SSR << ", " << DF << ", " << P_VAR(i)*DF/( SSR - P_VAR(i) ) << "\n\n";
					}*/
					
					if( PVAL >= 1){
						N_PVAL1 += 1;
					}else if( PVAL < 0 || std::isnan(PVAL) ){
						
						std::cerr << "PVALUE OUT OF BOUNDS:\n";
						std::cerr << PVAL << ", " << P_VAR(i) << ", " << SSR << ", " << DF << ", " << P_VAR(i)*DF/( SSR - P_VAR(i) ) << "\n\n";
						abort();
						
					}else if( PVAL <= 1e-300 ){
						PVAL = 1e-300;
						NUMER += qcauchy(PVAL);
					}else{
						NUMER += qcauchy(PVAL);
					}
					
					if( P_VAR(i) > MAX_P_VAR ){
						MAX_P_VAR = P_VAR(i);
						
						k = i;
					}
				}
			}
		}
		// std::cout << U(k)/V(k) << "\t" << DF << "\t" << SSR << "\t" << P_VAR(k) << "\n";
		if( DENOM <= 0 ) DENOM++;
		double NUDGED_PVAL1 = DENOM >= 4 ? DENOM/(DENOM + 1) : 0.80;
		
		if( NUDGED_PVAL1 >= 1 - 1e-8 ){
			NUDGED_PVAL1 = 1 - 1e-8;
		}
		NUMER += N_PVAL1 * qcauchy( NUDGED_PVAL1 );
		if ( DENOM <= 0 || std::isnan(NUMER) || NUMER == 0.0 ){
			acat_p = -99;
			min_p = -99;
		}else{
			acat_p = pcauchy(NUMER/DENOM);
			double STAT = MAX_P_VAR*DF/( SSR - MAX_P_VAR );
			if( STAT > 0 && !std::isnan(STAT) && DF > 2 ){
				min_p = pf(STAT, 1.0, DF - 1, true );
			}
		}
	}
	
};


class meta_svar_sumstat
{
	private:
		std::vector<svar_sumstat> ss;
		
		vcov_getter& vg;
		std::vector<double>& ivw;
		std::vector<Eigen::MatrixXd> covar;
		int nvar;
		int n_studies;
		
		bool init_0;
		
	public:
		std::vector<int> kept_snps;
		
		svar_sumstat ss_meta;
		svar_sumstat ss_meta_0;
	
		// void triform_pval(const int& k, double& pval_hom, double& pval_het, double& pval_acat, double& pval_omni, const std::vector<svar_sumstat>& vss = ss, const svar_sumstat& mss = ss_meta){
		void triform_pval(const int& k, double& pval_hom, double& pval_het, double& pval_acat, double& pval_omni){
			
			if( ss_meta.V(k) <= 0 || ss_meta.V_0(k) <= 0 || ss_meta.V(k)/ss_meta.V_0(k) < 1.0 - global_opts::RSQ_PRUNE ){
				pval_hom = -99; pval_het = -99;
				pval_acat = -99; pval_omni = -99;
				return;
			}
			
			// Calculate homogeneous-effect p-value.
			pval_hom = ss_meta.single_snp_pval(k);
			
			double numer_het = 0.0, denom_het = 0.0, df_tot = 0.0, q_acat = 0.0, ns = 0.0;
			for( const auto& ss_i : ss ){
				
				// Skip study if the variant is monomorphic.
				if( ss_i.V_0(k) > 0  ){
				
					// Skip variant if it fails filters in any study.
					if( ss_i.V(k) <= 0 || ss_i.V(k)/ss_i.V_0(k) < 1.0 - global_opts::RSQ_PRUNE ){
						pval_hom = -99; pval_het = -99;
						pval_acat = -99; pval_omni = -99;
						return;
					}
				
					double U2_V_i = ss_i.U(k) * ss_i.U(k) / ss_i.V(k);
					double fstat = (ss_i.DF - 1)*U2_V_i/(ss_i.SSR - U2_V_i);
					
					if( !std::isnan(fstat) && fstat > 0 && ss_i.DF > 2 ){
					
						numer_het += U2_V_i;
						denom_het += ss_i.SSR;
						df_tot += ss_i.DF;
						
						q_acat += qcauchy(pf(fstat, 1.0, ss_i.DF - 1, true));
						
						ns++;
					}else{
						pval_hom = -99; pval_het = -99;
						pval_acat = -99; pval_omni = -99;
						return;
					}
				}
			}

			// Calculate heterogeneous-effect p-value.
			denom_het -= numer_het;
			numer_het /= ns;
			denom_het /= (df_tot - ns);
			
			double fstat_het = numer_het/denom_het;
			if( !std::isnan(fstat_het) && fstat_het > 0 && df_tot - ns > 2 && ns > 0 ){
				pval_het = pf( fstat_het, ns, df_tot - ns, true );
			}else{
				pval_het = -99;
			}
			
			// Calculate ACAT p-value.
			if ( ns > 0 ){
				q_acat /= ns;
				pval_acat = pcauchy(q_acat);
			}else{
				pval_acat = -99;
			}
			
			// Calculate omnibus p-value
			double n_omni = 0.0, stat_omni = 0.0, min_pval = 1.0;
			
			if( global_opts::het_use_hom ){
				if( pval_hom >= 0 ){
					if( pval_hom < min_pval ){
						min_pval = pval_hom;
					}
					stat_omni += qcauchy(pval_hom);
					n_omni++;
				}
			}
			
			if( global_opts::het_use_het ){
				if( pval_het >= 0 ){
					if( pval_het < min_pval ){
						min_pval = pval_het;
					}
					stat_omni += qcauchy(pval_het);
					n_omni++;
				}
			}
			
			if( global_opts::het_use_acat ){
				if( pval_acat >= 0 ){
					if( pval_acat < min_pval ){
						min_pval = pval_acat;
					}
					stat_omni += qcauchy(pval_acat);
					n_omni++;
				}
			}
			
			if( n_omni > 0 ){
				pval_omni = pcauchy(stat_omni/n_omni);
				
				// Switch to Bonferroni if something went wrong
				if( pval_omni > n_omni*min_pval ){
					pval_omni = n_omni*min_pval;
				}
			}else{
				pval_omni = -99;
			}
			
			return;
		};
		
		void triform_pval(double& pval_hom, double& pval_het, double& pval_acat, double& pval_omni, const std::vector<ss_lm_single>& vss, const ss_lm_single& mss){
			
			// Calculate homogeneous-effect p-value.
			pval_hom = mss.pval;
			
			double numer_het = 0.0, denom_het = 0.0, df_tot = 0.0, q_acat = 0.0, ns = 0.0;
			for( const auto& ss_i : vss ){
				
				double U2_V_i = ss_i.U * ss_i.U / ss_i.V;
				double fstat = (ss_i.df - 1)*U2_V_i/(ss_i.SSR - U2_V_i);
				
				if( !std::isnan(fstat) && fstat > 0 && ss_i.df > 2 ){
				
					numer_het += U2_V_i;
					denom_het += ss_i.SSR;
					df_tot += ss_i.df;
					
					q_acat += qcauchy(pf(fstat, 1.0, ss_i.df - 1, true));
					
					ns++;
				}else{
					pval_hom = -99; pval_het = -99;
					pval_acat = -99; pval_omni = -99;
					return;
				}
			}

			// Calculate heterogeneous-effect p-value.
			denom_het -= numer_het;
			numer_het /= ns;
			denom_het /= (df_tot - ns);
			
			double fstat_het = numer_het/denom_het;
			if( fstat_het > 0 && df_tot - ns > 2 && ns > 0 ){
				pval_het = pf( fstat_het, ns, df_tot - ns, true );
			}else{
				pval_het = -99;
			}
			
			// Calculate ACAT p-value.
			if ( ns > 0 ){
				q_acat /= ns;
				pval_acat = pcauchy(q_acat);
			}else{
				pval_acat = -99;
			}
			
			// Calculate omnibus p-value
			double n_omni = 0.0, stat_omni = 0.0, min_pval = 1.0;
			
			if( global_opts::het_use_hom ){
				if( pval_hom >= 0 ){
					if( pval_hom < min_pval ){
						min_pval = pval_hom;
					}
					stat_omni += qcauchy(pval_hom);
					n_omni++;
				}
			}
			
			if( global_opts::het_use_het ){
				if( pval_het >= 0 ){
					if( pval_het < min_pval ){
						min_pval = pval_het;
					}
					stat_omni += qcauchy(pval_het);
					n_omni++;
				}
			}
			
			if( global_opts::het_use_acat ){
				if( pval_acat >= 0 ){
					if( pval_acat < min_pval ){
						min_pval = pval_acat;
					}
					stat_omni += qcauchy(pval_acat);
					n_omni++;
				}
			}
			
			if( n_omni > 0 ){
				pval_omni = pcauchy(stat_omni/n_omni);
				
				// Switch to Bonferroni if something went wrong
				if( pval_omni > n_omni*min_pval ){
					pval_omni = n_omni*min_pval;
				}
			}else{
				pval_omni = -99;
			}
			
			return;
		};
		
		void omni_min_pval(int& k, double& min_omni_p, double& acat_omni_p){
			
			double MAX_P_VAR = -1.0, NUMER = 0.0, DENOM = 0.0, N_PVAL1 = 0.0;
			
			k = -1;
			min_omni_p = -99; acat_omni_p = -99;
			
			for(int i = 0; i < ss_meta.U.size(); i++){
				if( ss_meta.V(i) > 0 && ss_meta.V_0(i) > 0 ){
					if( ss_meta.V(i)/ss_meta.V_0(i) > 1.0 - global_opts::RSQ_PRUNE ){
						
						double hom_p, het_p, aca_p, omn_p;
						triform_pval(i, hom_p, het_p, aca_p, omn_p);
						
						if( omn_p >= 1){
							DENOM++;
							N_PVAL1++;
						}else if( omn_p > 0 && !std::isnan(omn_p) ){
							DENOM++;
							NUMER += qcauchy(omn_p);
							if( (omn_p < min_omni_p && omn_p > 0) || min_omni_p < 0 ){
								min_omni_p = omn_p;
								k = i;
							}
						}
					}
				}
			}
			// std::cout << U(k)/V(k) << "\t" << DF << "\t" << SSR << "\t" << P_VAR(k) << "\n";
			if( DENOM <= 0 || NUMER == 0.0 ){
				k = -1; min_omni_p = -99; acat_omni_p = -99;
			}else{
				double NUDGED_PVAL1 = DENOM >= 4 ? DENOM/(DENOM + 1) : 0.80;
				NUMER += N_PVAL1 * qcauchy( NUDGED_PVAL1 );
				acat_omni_p = pcauchy(NUMER/DENOM);
			}
			
			return;
		}
		
		meta_svar_sumstat(vcov_getter& vg_, std::vector<double>& ivw_) : init_0(false), vg(vg_), ivw(ivw_) {n_studies = ivw.size();};
		
		void add_sstat(const Eigen::VectorXd& U_, const Eigen::VectorXd& V_, const double& DF_, const double& SSR_, const double& SCALE_Y_SD_ = 1.0){
			// std::cout << DF_ << "\t" << SSR_ << "\n";
			ss.push_back( svar_sumstat(U_,V_,DF_,SSR_,SCALE_Y_SD_) );
			nvar = ss[0].U.size();
		};
		
		void update_meta_ss(){
			Eigen::VectorXd U_ = Eigen::VectorXd::Zero(nvar);
			Eigen::VectorXd V_ = Eigen::VectorXd::Zero(nvar);
			double DF_ = 0;
			double SSR_ = 0;
			double SCALE_Y_SD_ = 0;
			double DENOM = 0;
			for( int i = 0; i < n_studies; i++ ){
				U_ += ivw[i] * ss[i].U;
				V_ += ivw[i] * ss[i].V;
				DF_ += ss[i].DF;
				SSR_ += ivw[i] * ss[i].SSR;
				SCALE_Y_SD_ += ivw[i] * (ss[i].SCALE_Y_SD * ss[i].SCALE_Y_SD);
				DENOM += ivw[i];
			}
			U_ /= DENOM;
			V_ /= DENOM;
			SSR_ /= DENOM;
			SCALE_Y_SD_ = std::sqrt(SCALE_Y_SD_/DENOM);
			
			// std::cout << DF_ << "\n";
			
			ss_meta = svar_sumstat(U_,V_,DF_,SSR_,SCALE_Y_SD_);
			
			if( init_0 ){
				ss_meta.U_0 = ss_meta_0.U_0;
				ss_meta.V_0 = ss_meta_0.V_0;
			}else{
				ss_meta_0 = ss_meta;
				init_0 = true;
			}
		};
		
		void condition_on_het(const int& k);
		
		void drop_snp(const int& i){
			for( svar_sumstat& ss_i : ss ){
				ss_i.drop_snp(i);
			}
			update_meta_ss();
			kept_snps.erase(kept_snps.begin() + i);
		}
		
		void condition_on_subset(const std::vector<int>& vi){
			for( svar_sumstat& ss_i : ss ){
				ss_i.condition_on_subset(vi);
			}
		}
		
		ss_lm_single final_model(const int& c_i){
			std::vector<ss_lm_single> lm_singles;
			for( svar_sumstat& ss_i : ss ){
				lm_singles.push_back(ss_i.final_model(c_i));
			}
			return ss_lm_single(lm_singles);
		}
		
		ss_lm_single final_model_triform_pval(const int& k, double& pval_hom, double& pval_het, double& pval_acat, double& pval_omni){
			std::vector<ss_lm_single> lm_singles;
			for( svar_sumstat& ss_i : ss ){
				lm_singles.push_back(ss_i.final_model(k));
			}
			ss_lm_single lm_meta = ss_lm_single(lm_singles);
			triform_pval(pval_hom, pval_het, pval_acat, pval_omni, lm_singles, lm_meta);
			return lm_meta;
		};
				
		ss_lm_single marginal_triform_pval(const int& k, double& pval_hom, double& pval_het, double& pval_acat, double& pval_omni){
			std::vector<ss_lm_single> lm_singles;
			for( svar_sumstat& ss_i : ss ){
				lm_singles.push_back(ss_i.marginal_model(k));
			}
			ss_lm_single lm_meta = ss_lm_single(lm_singles);
			triform_pval(pval_hom, pval_het, pval_acat, pval_omni, lm_singles, lm_meta);
			return lm_meta;
		};
};

class lm_output
{
	public:
		std::vector<std::string> variable;
		
		std::vector<double> beta;
		std::vector<double> se;
		std::vector<double> pval;
		
		double df;
		std::string name;
		std::string info;
		
		void push_back(double,double,double);
		void print_coefs();
};


class forward_lm
{
	public:
		std::vector<std::string> variable;
		
		std::vector<double> beta;
		std::vector<double> se;
		std::vector<double> beta_0;
		std::vector<double> se_0;
		std::vector<double> pval_0;
		std::vector<double> pval_seq;
		std::vector<double> pval_joint;
		std::vector<double> pval_adj;
		
		std::vector<std::vector<int>> buddy_list;
		std::vector<std::vector<double>> corr_list; 
		
		std::vector<int> keep;
		std::vector<std::vector<int>> conditioned;
		
		double df;
		std::string name;
		std::string info;
		
		forward_lm(const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& df_0, const double& stdev, vcov_getter& vget, double pval_thresh);
		
		void check_joint_pvalues(int&, double&, const Eigen::VectorXd&, const Eigen::VectorXd&,const Eigen::VectorXd&, const Eigen::MatrixXd&, const double&, const double&);
		
		void push_back(double,double,double);
		void print_coefs();
};

lm_output lm_from_sumstats( const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& m, const double& stdev = 1.0, const Eigen::VectorXd& U_0 = vec0, const Eigen::MatrixXd& J_0 = mat0, const Eigen::MatrixXd& Cov = vec0, const bool& check_filter = true, const std::vector<bool>& exclude = std::vector<bool>(0));

#endif

