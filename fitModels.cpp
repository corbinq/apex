#include "fitModels.hpp"

using namespace std;

double get_neg_logLik_REML(const double& delta, Eigen::MatrixXd& X_tilde, Eigen::VectorXd& y_tilde, Eigen::VectorXd& lambda )
{
	
	double df_resid = X_tilde.rows() - X_tilde.cols();

	Eigen::VectorXd vals = 1.00 + delta*lambda.array();
	
	Eigen::DiagonalMatrix<double,Eigen::Dynamic> Di = vals.asDiagonal().inverse();
	
	Eigen::MatrixXd XtDX = X_tilde.transpose() * Di * X_tilde;
	Eigen::VectorXd XtDy = X_tilde.transpose() * Di * y_tilde;
	Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
	
	double sigma2 = (y_tilde.dot(Di * y_tilde) - XtDy.dot(b))/df_resid;
	
	double ll = -0.5*(
		df_resid*log(sigma2) + 
		vals.array().log().sum() + log(XtDX.determinant()) + 1.0
	);
	
	return ll;
}


void meta_svar_sumstat::condition_on_het(const int& k){
	kept_snps.push_back(k);
	for( int s = 0; s < ss.size(); s++ ){
		//Eigen::MatrixXd cov_mat = vg.Covar_perStudy(s, k);
		//cout << cov_mat << "\n";
		ss[s].condition_on( vector<int>(1, k), vg.Covar_perStudy(s, k) );
		// cout << "Conditioned " << s << "\n";
	}
	update_meta_ss();
}

lm_output lm_from_sumstats( const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& df_0, const double& stdev, const Eigen::VectorXd& U_0, const Eigen::MatrixXd& J_0, const Eigen::MatrixXd& Cov, const bool& check_filter, const vector<bool>& exclude ){
	
	lm_output out;
	
	Eigen::MatrixXd Ji;
	Eigen::VectorXd JiU_0;
	Eigen::VectorXd C_JiU_0;
	
	double m_0 = U_0.size();
	
	double SSE_0 = df_0;
	
	double df = df_0 - m_0 - 1;
	
	int i_m = 0;
	
	if( m_0 > 0 ){
		
		Ji = J_0.inverse();

		JiU_0 = Ji * U_0;

		C_JiU_0 = Cov * JiU_0;

		SSE_0 = df_0 - U_0.dot(JiU_0);
	}
	
	for( int i = 0; i < U.size(); ++i){
		bool skip = false;
		if( exclude.size() > 0 ) skip = exclude[i];
		if( skip && check_filter ){
			// cerr << "\nskip\n\n";
			out.push_back(-99, -99, -99);
		}else{
			double SCORE, VARSC, SSE_i;
			SSE_i = SSE_0;
			if( m_0 > 0 ){
				SCORE = U(i) - C_JiU_0(i);
				VARSC = V(i) - Cov.row(i)*Ji*Cov.row(i).transpose();
			}else{
				SCORE = U(i);
				VARSC = V(i);
			}
			if( (VARSC > 0 && (VARSC/V(i)) > 1 - global_opts::RSQ_PRUNE && V(i) > 0) || !check_filter ){
				double beta = stdev * SCORE/VARSC;
				
				double se = stdev * sqrt(SSE_i - SCORE*SCORE/VARSC) / sqrt(df*VARSC);
				double pval = -99;
				
				double PVAR = SCORE*SCORE/VARSC;
				double STAT = df*PVAR/(SSE_i - PVAR);
				
				if( SSE_i - PVAR > 0 && STAT > 0 && PVAR > 0 ){
					pval = pf( STAT, 1.0, df, true );
				}
				out.push_back(beta, se, pval);
			}else{
				// cerr << "\nWARNING: RSQ_VIF = " << (VARSC/V(i)) << ", VARSC = "<< VARSC << "\n\n";
				out.push_back(-99, -99, -99);
			}
		}
	}
	
	return out;
}

int which_min( const vector<double>& p, bool gt0 ){
	double mp = p[0];
	int wm = -1;
	for(int i = 0; i < p.size(); ++i){
		if( gt0 && mp < 0 ){
			mp = p[i];
		}
		if( p[i] <= mp ){
			if( !(gt0 && p[i] <= 0) ){
				wm = i;
				mp = p[i];
			}
		}
	}
	return wm;
}


void forward_lm::check_joint_pvalues(int& index_of_largest_pvalue, double& largest_pvalue, const Eigen::VectorXd& U, const Eigen::VectorXd& V, const Eigen::VectorXd& U_0, const Eigen::MatrixXd& J_0, const double& n, const double& m){
	
	index_of_largest_pvalue = -1;
	largest_pvalue = 0;
	
	int k_i = 0;
	vector<int> kept_snps = seq_int(keep.size());
	for(const int k : keep ){
		vector<int> kept_snps_not_k = kept_snps;
		kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );
	
		double current_pvalue;

		Eigen::VectorXd U_k = U(vector<int>(1,k));
		Eigen::VectorXd V_k = V(vector<int>(1,k));
			
		Eigen::VectorXd U_0k = U_0(kept_snps_not_k); 
		Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k); 
		Eigen::MatrixXd Cov_k = J_0(vector<int>(1,k_i), kept_snps_not_k);
		
		lm_output reg_k = lm_from_sumstats(U_k, V_k, n, m, 1.00, U_0k, J_0k, Cov_k, false);
		
		current_pvalue = reg_k.pval[0];

		if( current_pvalue > largest_pvalue || current_pvalue < 0 ){
			largest_pvalue = current_pvalue;
			index_of_largest_pvalue = k_i;
		}
		
		k_i++;
	}
}


forward_lm::forward_lm(const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& m, const double& stdev, vcov_getter& vget, double pval_thresh )
{

	Eigen::VectorXd U_0 = Eigen::VectorXd(0);
	Eigen::MatrixXd J_0 = Eigen::MatrixXd(0,0);
	Eigen::MatrixXd Cov = Eigen::VectorXd(0);
	
	if( U.size() != V.size() ){
		cerr << "U.size() != V.size()" << "\n";
		exit(1);
	}
	
	if( U.size() <= 0 ){
		cerr << "U.size() <= 0" << "\n";
		exit(1);
	}
	
	int n_var = U.size();
	int nk = 0;
	
	vector<bool> excl(n_var,false);
	
	lm_output reg0;
	
	double alpha_thresh = global_opts::LM_ALPHA;
	
	while( 1 )
	{
		
		int steps_taken = 0;
		
		lm_output reg = lm_from_sumstats(U, V, n, m, stdev, U_0, J_0, Cov, true, excl);
		// cerr << "Fit model.\n";

		if( nk == 0 && reg0.beta.size() == 0 ){
			reg0 = reg;
		}

		int wk = which_min(reg.pval, true);
		double adj_pval = ACAT_non_missing(reg.pval);
		
		//cout << nk << ":" << wk << "\t" << reg0.beta[wk]  << "\t" << reg0.se[wk]  << "\t" << reg.pval[wk]  << "\t" << adj_pval << "\n";
		
		double pval_check = -99;
		
		if( wk >= 0 ){
			if( global_opts::step_marginal ){
				pval_check = reg.pval[wk];
			}else{
				pval_check = adj_pval;
			}
		}
		
		// -----------------------------------
		// Forward step. 
		// -----------------------------------
		if(  (pval_check >= 0 && pval_check < alpha_thresh) || (keep.size() == 0 && wk >= 0 )  )
		{
			keep.push_back(wk);
			// excl[wk] = true;
			
			beta_0.push_back(reg0.beta[wk]);
			se_0.push_back(reg0.se[wk]);
			pval_0.push_back(reg0.pval[wk]);
			pval_seq.push_back(reg.pval[wk]);
			pval_adj.push_back(adj_pval);
			
			nk++;
			
			if( nk == 1 ){
				U_0 = Eigen::VectorXd(1);
				Cov = Eigen::MatrixXd(n_var,1);

			}else{
				U_0.conservativeResize(nk);
				Cov.conservativeResize(Eigen::NoChange, nk);
			}
			

			//get_var(keep, J_0);
			J_0 = vget.Var(keep);

			U_0(nk-1) = U(wk);
			
			if( pval_check <= alpha_thresh || keep.size() < 2 ){
				// Eigen::VectorXd new_cov; 
				// get_cov(wk, new_cov);
				Eigen::VectorXd new_cov = vget.Covar(wk);


				if( new_cov.size() != Cov.rows() )
				{
					cerr << "new_cov.size() != Cov.rows()" << "\n";
					exit(1);
				}
				
				vector<int> b_list;
				vector<double> r_list;
				
				for(int i = 0; i < n_var; ++i)
				{
					
					Cov(i, nk-1) = new_cov(i);
					
					if( global_opts::RSQ_BUDDY < 1.00 ){
						double corr = new_cov(i)*new_cov(i)/(V(i)*V(wk));
						
						// if( corr*corr > global_opts::RSQ_PRUNE ) excl[i] = true;
						
						if( corr* corr > global_opts::RSQ_BUDDY ){
							b_list.push_back(i);
							r_list.push_back(corr);
						}
					}
				}
				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.push_back(b_list);
					corr_list.push_back(r_list);
				}
			}
			
			steps_taken++;
			//cout << "forward: " << nk << " SNPs. Added SNP " << wk << "\n";
		}
		
		// -----------------------------------
		// Backward step. 
		// -----------------------------------
		if( global_opts::backward_step && nk > 1 ){

			double max_joint_pvalue = 0; 
			int k_rm = -1;
			
			check_joint_pvalues(k_rm, max_joint_pvalue, U, V, U_0, J_0, n, m);
			//cout << max_joint_pvalue << "\n";
			
			if( (max_joint_pvalue > alpha_thresh || max_joint_pvalue < 0 ) && k_rm < nk - 1 && k_rm >= 0 ){
				
				vector<int> kept_snps_not_k_rm = seq_int(nk);
				kept_snps_not_k_rm.erase( kept_snps_not_k_rm.begin() + k_rm );
		
				U_0 = (U_0(kept_snps_not_k_rm)).eval(); 
				J_0 = (J_0(kept_snps_not_k_rm, kept_snps_not_k_rm)).eval(); 
				Cov = (Cov(seq_int(n_var), kept_snps_not_k_rm)).eval();
					
				beta_0.erase( beta_0.begin() + k_rm );
				se_0.erase( se_0.begin() + k_rm );
				pval_0.erase( pval_0.begin() + k_rm );
				
				pval_seq.erase( pval_seq.begin() + k_rm );
				pval_adj.erase( pval_adj.begin() + k_rm );
				
				if( global_opts::RSQ_BUDDY < 1.00 ){
					buddy_list.erase(buddy_list.begin() + k_rm);
					corr_list.erase(corr_list.begin() + k_rm);
				}
				
				// excl[keep[k_rm]] = false;
				
				keep.erase( keep.begin() + k_rm );
				
				nk--;
				
				// cout << "backward\n";
				steps_taken++;
			}
		}
		// cout << steps_taken << "\n";
		
		if( steps_taken == 0 ){
			break;
		}
	}

	int k_i = 0;
	vector<int> kept_snps = seq_int(nk);
	
	for( const int k : keep ){
		
		vector<int> kept_snps_not_k = kept_snps;
		kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );
	
		if( kept_snps_not_k.size() == 0 ){
			
			beta.push_back(reg0.beta[k]);
			se.push_back(reg0.se[k]);
			pval_joint.push_back(reg0.pval[k]);

		}else{
			
			Eigen::VectorXd U_k = U(vector<int>(1,k));
			Eigen::VectorXd V_k = V(vector<int>(1,k));
				
			Eigen::VectorXd U_0k = U_0(kept_snps_not_k); 
			Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k); 
			Eigen::MatrixXd Cov_k = J_0(vector<int>(1,k_i), kept_snps_not_k);
			
			lm_output reg_k = lm_from_sumstats(U_k, V_k, n, m, stdev, U_0k, J_0k, Cov_k, false);
			
			beta.push_back(reg_k.beta[0]);
			se.push_back(reg_k.se[0]);
			pval_joint.push_back(reg_k.pval[0]);
		}

		conditioned.push_back(kept_snps);
		k_i++;
		
	}
	
}

void lm_output::push_back(double b, double s, double p)
{
	beta.push_back(b);
	se.push_back(s);
	pval.push_back(p);
}

void lm_output::print_coefs()
{
	for(int i = 0; i < pval.size(); ++i){
		cout << 
			beta[i] << "\t" << 
			se[i] << "\t" << 
			pval[i] << "\n";
	}
}

