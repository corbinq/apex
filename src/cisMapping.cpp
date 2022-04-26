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


#include "Main.hpp"

void run_cis_QTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	std::cerr << "Started cis-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);

	make_half_hat_matrix(X);

	if( !global_opts::low_mem ){
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = X.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = std::sqrt(g_data.var[i]);
		}
		std::cerr << "Done.\n";
	}
	
	std::cerr << "Calculating expression residuals...\n";
	make_resid_from_half_hat(Y, X);
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y, e_data.stdev);
	
	if( rknorm_r ){
		std::cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y);
		std::cerr << "Re-residualizing transformed residuals ...\n";
		make_resid_from_half_hat(Y, X);
		scale_and_center(Y);
	}
	
	// std::cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y.transposeInPlace();
	
	double n_samples = X.rows();
	double n_covar = X.cols();
	
	std::string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	std::string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
	BGZF* block_file;
	BGZF* bed_block_file;
	BGZF* long_file;
	
	if( !just_long ){
		
		block_file = bgzf_open(block_file_path.c_str(), "w");
		
		bed_block_file = bgzf_open(bed_block_file_path.c_str(), "w");
		
		write_to_bgzf("#chrom\tstart\tend\tgene\tegene_pval\tn_samples\tn_covar\tresid_sd\tn_cis_variants\n", bed_block_file);
		
	}
	
	bool write_long = (make_long || just_long);
	
	if( write_long ){
		long_file = bgzf_open(long_file_path.c_str(), "w");
		write_to_bgzf("#chrom\tpos\tref\talt\tgene\tbeta\tse\tpval\n", long_file);
	}
	
	int bl = 0;
	
	std::string iter_cerr_suffix = " cis-QTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y.rows() ? n_e : Y.rows()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y.rows() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (X.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
				{
					g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}

				out_mat = (Y.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
			}else{
				
				const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(bm.bcf_s[i], n_g);
				
				out_mat = (Y.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
				// out_mat = (Y.middleRows(bm.bed_s[i], n_e)* g_data.genotypes.middleCols(bm.bcf_s[i], n_g)).transpose().eval();
				
			}
			
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; ++jj, ++jm){
				
				std::stringstream block_line;
				
				std::vector<long double> pvals;
				std::vector<double> dist;
				double gene_tss = 0.5*(e_data.start[jj] + e_data.end[jj]);
		
				int pos_s = -1;
				int pos_e = -1;
		
				int v_s = -1;
				int v_e = -1;
				
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, ++im){
				
					std::stringstream long_line;
				
					if( g_data.chr[ii] == e_data.chr[jj] && e_data.start[jj] - g_data.pos[ii] < global_opts::cis_window_bp && g_data.pos[ii] - e_data.end[jj] < global_opts::cis_window_bp ){
						
						if( pos_s < 0 ) pos_s = g_data.pos[ii]; 
						pos_e = g_data.pos[ii];

						if( v_s < 0 ) v_s = ii; 
						v_e = ii;

						double g_yres_crossprod = out_mat(im, jm);
						if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
						
						if( !just_long ) block_line << "\t" << g_yres_crossprod;
						
						double beta, beta_se, zscore, log_pval_esnp;
						long double pval_esnp;

						if( g_data.var[ii] > 0 ){
							beta = g_yres_crossprod/g_data.var[ii];
							beta_se = std::sqrt( (n_samples - 1)/g_data.var[ii] - beta*beta)/std::sqrt(n_samples - n_covar - 1);
							zscore = beta/beta_se;
							log_pval_esnp = rmath::pf( zscore*zscore, 1.0, n_samples - n_covar - 1, false, true );
              pval_esnp = exp(log_pval_esnp);
						}else{
							beta = 0;
							beta_se = 0;
							zscore = 0;
              log_pval_esnp = 1;
							pval_esnp = -1;
						}

						dist.push_back( std::abs(g_data.pos[ii] - gene_tss) );
						pvals.push_back(pval_esnp);

						if( write_long ){
              std::string str_pval_esnp = log_to_string(log_pval_esnp);

							long_line << 
								clean_chrom(g_data.chr[ii]) << "\t" <<
								g_data.pos[ii] << "\t" <<
								g_data.ref[ii] << "\t" <<
								g_data.alt[ii] << "\t" <<
								e_data.gene_id[jj] << "\t" <<
								y_scale[jj] * e_data.stdev[jj] * beta << "\t" <<
								y_scale[jj] * e_data.stdev[jj] * beta_se  << "\t" <<
                str_pval_esnp << "\n";
						}
					}
					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						std::string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							std::to_string(pos_s) + "\t" + 
							std::to_string(pos_e) + "\t" + 
							std::to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						std::stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals, dist) << "\t" << 
							X.rows() << "\t" << 
							X.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			std::cerr << "\nERROR: " <<bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	std::cerr << "\n";

	if ( write_long ){
		bgzf_close(long_file);
		//build_tabix_index(long_file_path);
	}

	if( !just_long ){
		bgzf_close(block_file);
		bgzf_close(bed_block_file);
		build_tabix_index(block_file_path, 1);
		build_tabix_index(bed_block_file_path, 1);
	}

	return;
}

void run_cis_QTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long, const std::string& theta_path, const std::string& anchor_path)
{
	
	/*
	Eigen::SparseMatrix<double> Q;
	Eigen::VectorXd Q_lambda;

	subset_eigen(L, GRM_lambda, Q, Q_lambda);

	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
		
	// std::cerr << "Reordering trait and covariate matrices ...\n";
	// std::cerr << "Reordering genotypes ...\n";
	// g_data.genotypes = (Tr * g_data.genotypes).eval();
		
	
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_snps = g_data.n_variants;
	double n_covar = C.cols();
	
	if ( n_traits < 2 ) {
		std::cerr << "ERROR: apex cis requires 2 or more traits (" << n_traits << " specified).\n";
		abort();
	}

	// std::cerr << "Started cis-QTL analysis ...\n";
	
	if( e_data.is_residual ){
		std::cerr << "ERROR: LMM residuals are not currently supported in mode cis. \n";
		abort();
	}
	
	
	if( rknorm_y ){
		if( e_data.is_residual ){
			std::cerr << "IGNORING instruction to rank-normalize LMM residuals. \n";
		}else{
			std::cerr << "Rank-normalizing expression traits ... \n";
			rank_normalize(Y);
		}
	}
	
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	
	// std::cerr << "Calculating genotype-covariate covariance ... \n";
	// Eigen::MatrixXd CtG = (C.transpose() * g_data.genotypes).eval();
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtC = (Q.transpose() * C).eval();
	Eigen::MatrixXd QtG = (Q.transpose() * g_data.genotypes).eval();
	Eigen::MatrixXd QtY = (Q.transpose() * Y).eval();
	Eigen::MatrixXd CtY = (C.transpose() * Y).eval();
	Eigen::MatrixXd CtC = (C.transpose() * C).eval();
	Eigen::MatrixXd CtC_i = (CtC.inverse()).eval();
	
	std::cerr << "Rotating expression and covariates ... ";
	Y = (L.transpose() * Y).eval();
	Eigen::MatrixXd X = (L.transpose() * C).eval();
	std::cerr << "Done.\n";
	
	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	
	Eigen::MatrixXd V_mat;
	calculate_V_anchor_points(V_mat, g_data, C, hsq_vals, CtC, CtC_i, QtG, QtC, Q_lambda);
	QtG.resize(0,0);
	
	*/
	
	
	if( e_data.is_residual ){
		std::cerr << "NOTE: BED file contains null-model LMM trait residuals. \n";
	}
	
	
	Eigen::SparseMatrix<double> Q;
	Eigen::VectorXd Q_lambda;

	subset_eigen(L, GRM_lambda, Q, Q_lambda);

	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
	
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_snps = g_data.n_variants;
	double n_covar = C.cols();

	if( rknorm_y ){
		if( e_data.is_residual ){
			std::cerr << "IGNORING instruction to rank-normalize LMM residuals ... \n";
		}else{
			std::cerr << "Rank-normalizing expression traits ... \n";
			rank_normalize(Y);
		}
	}
	
	std::vector<double> y_scale;
	
	Eigen::MatrixXd QtY, CtY, QtC, QtG, CtC, CtC_i;
	bool gvar_precomputed = ( anchor_path != "" );
	
	if( !e_data.is_residual ){
		
		std::cerr << "Scaling expression traits ... \n";
	
		scale_and_center(Y, y_scale);
		
		QtY = (Q.transpose() * Y).eval();
		CtY = (C.transpose() * Y).eval();
		Y = (L.transpose() * Y).eval();
	}else{
		
		y_scale = std::vector<double>( (int) n_traits, 1.00 );
	}
	
	if( !gvar_precomputed ||  !e_data.is_residual ){
		std::cerr << "Calculating partial rotations ...\n";
		
		QtC = (Q.transpose() * C).eval();
		CtC = (C.transpose() * C).eval();
		
		if( !gvar_precomputed ){
			QtG = (Q.transpose() * g_data.genotypes).eval();
			CtC_i = (CtC.inverse()).eval();
		}
		
		std::cerr << "Done.\n";
	}

	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	
	Eigen::MatrixXd V_mat, X;
	
	if( !gvar_precomputed ){
		calculate_V_anchor_points(V_mat, g_data, C, hsq_vals, CtC, CtC_i, QtG, QtC, Q_lambda);
		QtG.resize(0,0);
		CtC_i.resize(0,0);
	}else{
		read_V_anchor_points(anchor_path, V_mat, g_data, hsq_vals);
	}
	
	
	// Set up files
	
	// std::string theta_file_path = global_opts::out_prefix + "." + "theta" + ".gz";
	std::string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	std::string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
	// BGZF* theta_file;
	BGZF* block_file;
	BGZF* bed_block_file;
	BGZF* long_file;
	
	// theta_file = bgzf_open(theta_file_path.c_str(), "w");
	
	if( !just_long ){
		
		block_file = bgzf_open(block_file_path.c_str(), "w");
		
		bed_block_file = bgzf_open(bed_block_file_path.c_str(), "w");
		
		write_to_bgzf("#chrom\tstart\tend\tgene\tegene_pval\tn_samples\tn_covar\tresid_sd\tn_cis_variants\n", bed_block_file);
		
	}
	
	bool write_long = (make_long || just_long);
	
	if( write_long ){
		long_file = bgzf_open(long_file_path.c_str(), "w");
		write_to_bgzf("#chrom\tpos\tref\talt\tgene\tbeta\tse\tpval\n", long_file);
	}
	
	// -------------------------------------
	// Begin fitting null models. 
	// -------------------------------------
	
	std::vector<double> phi_v(Y.cols());
	std::vector<double> hsq_v(Y.cols());
	std::vector<double> sigma_v(Y.cols());
	std::vector<double> SSR_v(Y.cols());
	
	e_data.stdev.resize(Y.cols());
	
	theta_data t_data;
	bool use_theta = false;
	
	if( theta_path != "" ){
		std::cerr << "Set null model for ";
		t_data.open(theta_path);
		use_theta = true;
	}else{
		std::cerr << "Fit null model for ";
	}
	
	X = (L.transpose() * C).eval();
	
	std::string iter_cerr_suffix = " traits out of " + std::to_string(Y.cols()) + " total";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	int last_j = 0;
	for(int j = 0; j < Y.cols(); j++ ){
		DiagonalXd Vi;
		double sigma2;
		double phi;
		double tau2;

		if( use_theta ){
			t_data.getTheta(j, e_data.gene_id[j], sigma2, tau2, phi, SSR_v[j], y_scale[j]);
			Vi = calc_Vi(phi, GRM_lambda);
		}else{
			
			if( e_data.is_residual ){
				std::cerr << "Error: --theta must be specified when --bed contains LMM residuals.\n";
				abort();
			}
			
			LMM_fitter fit(X, Y.col(j), GRM_lambda);
			fit.fit_REML();
			
			Vi = fit.Vi;
			sigma2 = fit.sigma2;
			phi = fit.phi;
			tau2 = phi*sigma2;
		}
		
		// DiagonalXd Vi = Eigen::VectorXd::Ones(Y.col(j).size()).asDiagonal();
		// double sigma2 = 1.00;
		// double phi = 0.05;
		
		double hsq = tau2 / (tau2 + sigma2);
		double scale = std::sqrt(tau2 + sigma2);
		
		if( !e_data.is_residual ){
			DiagonalXd Psi = calc_Psi(phi, Q_lambda);
			Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC )/(1.00 + phi);
			Eigen::VectorXd XtDy = X.transpose() * Vi * Y.col(j);
			Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
			Eigen::VectorXd y_hat = X * b;
			
			Eigen::VectorXd y_res = Y.col(j) - y_hat;
			
			SSR_v[j] = y_res.dot(Vi * Y.col(j))/sigma2;
			 
			Y.col(j) = (Vi * y_res/std::sqrt(sigma2)).eval();
		}
		
		// DiagonalXd Psi = calc_Psi(phi, Q_lambda);
		// Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC )/(1.00 + phi);
		// Eigen::VectorXd XtDy = X.transpose() * Vi * Y.col(j);
		// Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
		// Eigen::VectorXd y_hat = X * b;
		
		// Eigen::VectorXd y_res = Y.col(j) - y_hat;
		
		// SSR_v[j] = y_res.dot(Vi * Y.col(j))/sigma2;
		
		// Y.col(j) = (Vi * y_res/std::sqrt(sigma2)).eval();
		
		
		phi_v[j] = phi;
		hsq_v[j] = hsq;
		sigma_v[j] = sigma2;
		
		e_data.stdev[j] = std::sqrt(sigma2);
		
		// thinned_iter_cerr(last_j, j+1, iter_cerr_suffix, 20);
		print_iter_cerr(j, j+1, iter_cerr_suffix);
	}
	//print_iter_cerr(last_j, Y.cols(), iter_cerr_suffix);
	std::cerr << "\n";
	
	if( !e_data.is_residual ){
		Y = (L * Y).eval();
	}
	
	// -------------------------------------
	// Start cis-QTL analysis
	// -------------------------------------
	
	int bl = 0;
	
	iter_cerr_suffix = " cis-QTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y.cols() ? n_e : Y.cols()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y.cols() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			
			// if( global_opts::low_mem ){
				// g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				// Eigen::SparseMatrix<double>& G = g_data.genotypes;
				// Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval();
				// for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; si++, ii++)
				// {
					// g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
				// }
			// }
			
			int s_slice = global_opts::low_mem ? 0 : bm.bcf_s[i];
			int n_slice = global_opts::low_mem ? g_data.genotypes.cols() : n_g;
			
			if( s_slice < 0 ){
				std::cerr << "\ns_slice = " << s_slice << "\n";
				continue;
			}else if( s_slice >= g_data.genotypes.cols() ){
				s_slice = g_data.genotypes.cols() - 1;
			}
			if( s_slice + n_slice >= g_data.genotypes.cols() ){
				n_slice = g_data.genotypes.cols() - s_slice;
			}
			
			if( n_slice <= 0 ){
				continue;
			}
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_slice, n_slice);
			
			for( int jj = bm.bed_s[i]; jj < bm.bed_s[i] + n_e; jj++){
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
				std::vector<double> dist;
				double gene_tss = 0.5*(e_data.start[jj] + e_data.end[jj]);
				
				int idx_s = -1;
				int idx_n = 0;
				
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, im++ ){
					if( g_data.chr[ii] == e_data.chr[jj] && e_data.start[jj] - g_data.pos[ii] < global_opts::cis_window_bp && g_data.pos[ii] - e_data.end[jj] < global_opts::cis_window_bp ){
						if( idx_s < 0 ) idx_s = im;
						idx_n++;
					}else{
						if( idx_s >= 0 ){
							break;
						}
					}
				}
				
				if( idx_s < 0 ){
					//std::cerr << "\nidx_s = " << idx_s << "\n";
					continue;
				}else if( idx_s >= G.cols() ){
					idx_s = G.cols() - 1;
				}
				if( idx_s + idx_n > G.cols() ){
					idx_n = G.cols() - idx_s;
				}
				if( n_slice <= 0 ){
					continue;
				}
				
				int v_s = idx_s + s_slice;
				int v_e = v_s + idx_n - 1;
				
				int pos_s = g_data.pos[v_s];
				int pos_e = g_data.pos[v_e];
				
				const Eigen::SparseMatrix<double>& G_slice = G.middleCols(idx_s, idx_n);

				const double& sigma2 = sigma_v[jj];
				const double& phi = phi_v[jj];
				
				double tau2 = phi*sigma2;
				double hsq = tau2 / (tau2 + sigma2);
				double S2 = tau2 + sigma2;
				
				// ------------------------------------------------
				
				// std::stringstream theta_line;
				
				// theta_line <<
					// clean_chrom(e_data.chr[jj]) << "\t" << 
					// e_data.start[jj] << "\t" << 
					// e_data.end[jj] << "\t" << 
					// e_data.gene_id[jj] << "\t" << 
					// sigma2 << "\t" << 
					// tau2 << "\t" << 
					// phi << "\n";
				
				// write_to_bgzf(theta_line.str().c_str(), theta_file);
				

				// ------------------------------------------------
				
				double ADJ = 1.00/(1.00 + phi);
				
				Eigen::VectorXd U_vec = G_slice.transpose() * Y.col(jj);
				const double& SSR_0 = SSR_v[jj];
				// ------------------------------------------------------------
				
				if(  v_s + U_vec.size() - 1 >=  V_mat.rows() ){
					std::cerr << "\n";
					std::cerr << " v_s = " << v_s << "\n";
					std::cerr << " U_vec.size() = " << U_vec.size() << "\n";
					std::cerr << " V_mat.rows() = " << V_mat.rows() << "\n";
					std::cerr << "\n";
					continue;
				}
				Eigen::VectorXd diagV = predV(V_mat.middleRows(v_s, U_vec.size()), hsq)/(1.00 + phi);
				
				e_data.stdev[jj] = std::sqrt(sigma2);
				
				for( int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;
				
					double g_yres_crossprod = U_vec(si);
					if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
					
					if( !just_long ) block_line << "\t" << g_yres_crossprod;
					
					double beta, beta_se, zscore, pval_esnp;
					
					if( g_data.var[ii] > 0.0 && diagV(si) > 0.0 ){
						beta = g_yres_crossprod/diagV(si);
						beta_se = std::sqrt( SSR_0/diagV(si) - beta*beta )/std::sqrt(n_samples - n_covar - 1);
						if( beta_se > 0 ){
							zscore = beta/beta_se;
							pval_esnp = pf( zscore*zscore, 1.0, n_samples - n_covar - 1, true );
						}else{
							zscore = 0;
							pval_esnp = -1;
						}
					}else{
						beta = 0;
						beta_se = 0;
						zscore = 0;
						pval_esnp = -1;
					}


					dist.push_back( std::abs(g_data.pos[ii] - gene_tss) );
					pvals.push_back(pval_esnp);

					if( write_long ){

						long_line << 
							clean_chrom(g_data.chr[ii]) << "\t" <<
							g_data.pos[ii] << "\t" <<
							g_data.ref[ii] << "\t" <<
							g_data.alt[ii] << "\t" <<
							e_data.gene_id[jj] << "\t" <<
							y_scale[jj] * e_data.stdev[jj]*beta << "\t" <<
							y_scale[jj] * e_data.stdev[jj]*beta_se  << "\t" <<
							pval_esnp << "\n";
					}

					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						std::string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							std::to_string(pos_s) + "\t" + 
							std::to_string(pos_e) + "\t" + 
							std::to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						std::stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals, dist) << "\t" << 
							X.rows() << "\t" << 
							X.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			std::cerr << "\nERROR: " << bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	std::cerr << "\n";

	// bgzf_close(theta_file);

	if ( write_long ){
		bgzf_close(long_file);
		//build_tabix_index(long_file_path);
	}

	if( !just_long ){
		bgzf_close(block_file);
		bgzf_close(bed_block_file);
		build_tabix_index(block_file_path, 1);
		build_tabix_index(bed_block_file_path, 1);
	}

	return;
}


void run_cis_QTL_analysis_eLMM(const int& n_fac, const int& n_fac_fe, bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long, Eigen::MatrixXd& Y_epc )
{
	
	Eigen::MatrixXd Q;
	Eigen::VectorXd Q_lambda;
	
	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
		
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_snps = g_data.n_variants;
	double n_covar = C.cols();

	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
		rank_normalize(Y_epc);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	scale_and_center(Y_epc);
	
	// bool resid_ePC = true;
	
	// if( resid_ePC ){
			
		// Eigen::MatrixXd U = get_half_hat_matrix(C);

		// std::cerr << "Calculating expression residuals...\n";
		// make_resid_from_half_hat(Y_epc, U);
	// }
	
	std::cerr << "Estimating latent factors ... \n";
	calc_eGRM_PCs(Q, Q_lambda, Y_epc, n_fac);
	std::cerr << "Done.\n";
	
	Y_epc.resize(0,0);
	
	if( n_fac_fe > 0 ){
		int n_cv = C.cols();
		C.conservativeResize(Eigen::NoChange, n_cv + n_fac_fe);
		for(int ii = 0; ii < n_fac_fe; ii++){
			C.col(n_cv + ii) = Q.col(ii);
		}
		n_covar = C.cols();
		
		// int n_re = n_fac - n_cv;
		
		// Q = (Q.rightCols(n_re)).eval();
		
		// Q_lambda = (Q_lambda.tail(n_re)).eval();
		// double max_lambda = Q_lambda.maxCoeff();
		// if( max_lambda <= 1.00 && max_lambda > 0.00 ){
			// Q_lambda /= max_lambda;
		// }
	}
	
	//Y = Y(Eigen::all, test_idx)
	
	// std::cerr << "Calculating genotype-covariate covariance ... \n";
	// Eigen::MatrixXd CtG = (C.transpose() * g_data.genotypes).eval();
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtC = (Q.transpose() * C).eval();
	Eigen::MatrixXd QtG = (Q.transpose() * g_data.genotypes).eval();
	Eigen::MatrixXd QtY = (Q.transpose() * Y).eval();
	Eigen::MatrixXd CtY = (C.transpose() * Y).eval();
	Eigen::MatrixXd CtC = (C.transpose() * C).eval();
	Eigen::MatrixXd CtC_i = (CtC.inverse()).eval();
	
	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	
	Eigen::MatrixXd V_mat;
	calculate_V_anchor_points_low_rank(V_mat, g_data, C, hsq_vals, CtC, CtC_i, QtG, QtC, Q_lambda);
	
	QtG.resize(0,0);
	
	std::string theta_file_path = global_opts::out_prefix + "." + "theta" + ".gz";
	std::string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	std::string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
	BGZF* theta_file;
	BGZF* block_file;
	BGZF* bed_block_file;
	BGZF* long_file;
	
	theta_file = bgzf_open(theta_file_path.c_str(), "w");
	
	if( !just_long ){
		
		block_file = bgzf_open(block_file_path.c_str(), "w");
		
		bed_block_file = bgzf_open(bed_block_file_path.c_str(), "w");
		
		write_to_bgzf("#chrom\tstart\tend\tgene\tegene_pval\tn_samples\tn_covar\tresid_sd\tn_cis_variants\n", bed_block_file);
		
	}
	
	bool write_long = (make_long || just_long);
	
	if( write_long ){
		long_file = bgzf_open(long_file_path.c_str(), "w");
		write_to_bgzf("#chrom\tpos\tref\talt\tgene\tbeta\tse\tpval\n", long_file);
	}
	
	// -------------------------------------
	// Begin fitting null models. 
	// -------------------------------------
	
	std::vector<double> phi_v(Y.cols());
	std::vector<double> hsq_v(Y.cols());
	std::vector<double> sigma_v(Y.cols());
	std::vector<double> SSR_v(Y.cols());
	
	// Eigen::MatrixXd PY( Y.rows(), Y.cols() );
	
	e_data.stdev.resize(Y.cols());
	
	theta_data t_data;
	bool use_theta = false;
	
	// if( theta_path != "" ){
		// std::cerr << "Set null model for ";
		// t_data.open(theta_path);
		// use_theta = true;
	// }else{
		std::cerr << "Fit null model for ";
	// }
	
	std::string iter_cerr_suffix = " traits out of " + std::to_string(Y.cols()) + " total";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	int last_j = 0;
	for(int j = 0; j < Y.cols(); j++ ){
		DiagonalXd Vi;
		double sigma2;
		double phi;
		double tau2;

		// if( use_theta ){
			// t_data.getTheta(j, e_data.gene_id[j], sigma2, tau2, phi);
			// Vi = calc_Vi(phi, GRM_lambda);
		// }else{
			// LMM_fitter fit(X, Y.col(j), GRM_lambda);
			
			double yty = Y.col(j).squaredNorm();
			double sample_size = Y.rows();
			
			LMM_fitter_low_rank fit(CtC, QtC, CtY.col(j), QtY.col(j), yty, Q_lambda, sample_size);
			
			fit.fit_REML();
			
			// Vi = fit.Vi;
			sigma2 = fit.sigma2;
			phi = fit.phi;
			tau2 = phi*sigma2;
		// }
		
		// DiagonalXd Vi = Eigen::VectorXd::Ones(Y.col(j).size()).asDiagonal();
		// double sigma2 = 1.00;
		// double phi = 0.05;
		
		double hsq = tau2 / (tau2 + sigma2);
		double scale = std::sqrt(tau2 + sigma2);
		
		DiagonalXd Psi = calc_Psi_low_rank(phi, Q_lambda);
		Eigen::MatrixXd XtDX = (CtC - QtC.transpose() * Psi * QtC );
		
		Eigen::VectorXd XtDy = QtC.transpose() * Psi * QtY.col(j);
		XtDy = (CtY.col(j) - XtDy).eval();
		Eigen::VectorXd b = XtDX.colPivHouseholderQr().solve(XtDy);
		Eigen::VectorXd y_fixeff = C * b;
		
		Eigen::VectorXd v1 = QtC * b;
		Eigen::VectorXd v2 = QtY.col(j) - v1;
		Eigen::VectorXd y_raneff = Q * Psi * v2;
		Eigen::VectorXd y_resid = Y.col(j) - y_raneff - y_fixeff;
		
		SSR_v[j] = (yty - QtY.col(j).dot(Psi * QtY.col(j)) -  XtDy.dot(b))/sigma2;
		
		Y.col(j) = (y_resid/std::sqrt(sigma2)).eval();
		
		phi_v[j] = phi;
		hsq_v[j] = hsq;
		sigma_v[j] = sigma2;
		
		e_data.stdev[j] = std::sqrt(sigma2);
		
		// thinned_iter_cerr(last_j, j+1, iter_cerr_suffix, 20);
		print_iter_cerr(j, j+1, iter_cerr_suffix);
	}
	//print_iter_cerr(last_j, Y.cols(), iter_cerr_suffix);
	std::cerr << "\n";
	
	// PY = (L * PY).eval();
	
	
	// -------------------------------------
	// Start cis-QTL analysis
	// -------------------------------------
	
	int bl = 0;
	
	iter_cerr_suffix = " cis-QTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y.cols() ? n_e : Y.cols()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y.cols() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			int s_slice = global_opts::low_mem ? 0 : bm.bcf_s[i];
			int n_slice = global_opts::low_mem ? g_data.genotypes.cols() : n_g;
			
			if( s_slice < 0 ){
				std::cerr << "\ns_slice = " << s_slice << "\n";
				continue;
			}else if( s_slice >= g_data.genotypes.cols() ){
				s_slice = g_data.genotypes.cols() - 1;
			}
			if( s_slice + n_slice >= g_data.genotypes.cols() ){
				n_slice = g_data.genotypes.cols() - s_slice;
			}
			
			if( n_slice <= 0 ){
				continue;
			}
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_slice, n_slice);
			
			for( int jj = bm.bed_s[i]; jj < bm.bed_s[i] + n_e; jj++){
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
				std::vector<double> dist;
				double gene_tss = 0.5*(e_data.start[jj] + e_data.end[jj]);
				
				int idx_s = -1;
				int idx_n = 0;
				
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, im++ ){
					if( g_data.chr[ii] == e_data.chr[jj] && e_data.start[jj] - g_data.pos[ii] < global_opts::cis_window_bp && g_data.pos[ii] - e_data.end[jj] < global_opts::cis_window_bp ){
						if( idx_s < 0 ) idx_s = im;
						idx_n++;
					}else{
						if( idx_s >= 0 ){
							break;
						}
					}
				}
				
				if( idx_s < 0 ){
					//std::cerr << "\nidx_s = " << idx_s << "\n";
					continue;
				}else if( idx_s >= G.cols() ){
					idx_s = G.cols() - 1;
				}
				if( idx_s + idx_n > G.cols() ){
					idx_n = G.cols() - idx_s;
				}
				if( n_slice <= 0 ){
					continue;
				}
				
				int v_s = idx_s + s_slice;
				int v_e = v_s + idx_n - 1;
				
				int pos_s = g_data.pos[v_s];
				int pos_e = g_data.pos[v_e];
				
				const Eigen::SparseMatrix<double>& G_slice = G.middleCols(idx_s, idx_n);

				const double& sigma2 = sigma_v[jj];
				const double& phi = phi_v[jj];
				
				double tau2 = phi*sigma2;
				double hsq = tau2 / (tau2 + sigma2);
				double S2 = tau2 + sigma2;
				
				// ------------------------------------------------
				
				// std::cout << e_data.gene_id[jj] << "\t" << fit.sigma2 << "\t" << tau2 << "\t" << hsq << "\n";
						
				std::stringstream theta_line;
				
				theta_line <<
					clean_chrom(e_data.chr[jj]) << "\t" << 
					e_data.start[jj] << "\t" << 
					e_data.end[jj] << "\t" << 
					e_data.gene_id[jj] << "\t" << 
					sigma2 << "\t" << 
					tau2 << "\t" << 
					phi << "\t" << SSR_v[jj] << "\n";
				
				write_to_bgzf(theta_line.str().c_str(), theta_file);
				
				
				// -------- Modified, similar to trans mode. -----------------
				Eigen::VectorXd U_vec = G_slice.transpose() * Y.col(jj);
				const double& SSR_0 = SSR_v[jj];
				// ------------------------------------------------------------
				
				if(  v_s + U_vec.size() - 1 >=  V_mat.rows() ){
					std::cerr << "\n";
					std::cerr << " v_s = " << v_s << "\n";
					std::cerr << " U_vec.size() = " << U_vec.size() << "\n";
					std::cerr << " V_mat.rows() = " << V_mat.rows() << "\n";
					std::cerr << "\n";
					continue;
				}
				Eigen::VectorXd diagV = predV(V_mat.middleRows(v_s, U_vec.size()), hsq);
				
				e_data.stdev[jj] = std::sqrt(sigma2);
				
				for( int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;
				
					double g_yres_crossprod = U_vec(si);
					if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
					
					if( !just_long ) block_line << "\t" << g_yres_crossprod;
					
					//std::cout << diagV(si) << "\t" << g_data.var[ii] << "\t" << g_yres_crossprod << "\n";
					
					double beta, beta_se, zscore, pval_esnp;
					
					if( g_data.var[ii] > 0.0 && diagV(si) > 0.0 ){
						beta = g_yres_crossprod/diagV(si);
						beta_se = std::sqrt( SSR_0/diagV(si) - beta*beta )/std::sqrt(n_samples - n_covar - 1);
						if( beta_se > 0 ){
							zscore = beta/beta_se;
							pval_esnp = pf( zscore*zscore, 1.0, n_samples - n_covar - 1, true );
						}else{
							zscore = 0;
							pval_esnp = -1;
						}
					}else{
						beta = 0;
						beta_se = 0;
						zscore = 0;
						pval_esnp = -1;
					}

					dist.push_back( std::abs(g_data.pos[ii] - gene_tss) );
					pvals.push_back(pval_esnp);

					if( write_long ){

						long_line << 
							clean_chrom(g_data.chr[ii]) << "\t" <<
							g_data.pos[ii] << "\t" <<
							g_data.ref[ii] << "\t" <<
							g_data.alt[ii] << "\t" <<
							e_data.gene_id[jj] << "\t" <<
							y_scale[jj] * e_data.stdev[jj]*beta << "\t" <<
							y_scale[jj] * e_data.stdev[jj]*beta_se  << "\t" <<
							pval_esnp << "\n";
					}

					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						std::string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							std::to_string(pos_s) + "\t" + 
							std::to_string(pos_e) + "\t" + 
							std::to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						std::stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals, dist) << "\t" << 
							C.rows() << "\t" << 
							C.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			std::cerr << "\nERROR: " << bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	std::cerr << "\n";

	bgzf_close(theta_file);

	if ( write_long ){
		bgzf_close(long_file);
		//build_tabix_index(long_file_path);
	}

	if( !just_long ){
		bgzf_close(block_file);
		bgzf_close(bed_block_file);
		build_tabix_index(block_file_path, 1);
		build_tabix_index(bed_block_file_path, 1);
	}

	return;
}


void run_cis_QTL_analysis_eFE(const int& n_fac, bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long, Eigen::MatrixXd& Y_epc )
{
	
	Eigen::MatrixXd Q;
	Eigen::VectorXd Q_lambda;
	
	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
		
	// std::cerr << "Reordering trait and covariate matrices ...\n";
	// std::cerr << "Reordering genotypes ...\n";
	// g_data.genotypes = (Tr * g_data.genotypes).eval();
		
	
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_snps = g_data.n_variants;
	double n_covar = C.cols();

	// std::cerr << "Started cis-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
		rank_normalize(Y_epc);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	scale_and_center(Y_epc);
	
	// make_half_hat_matrix(C);
	// make_resid_from_half_hat(Y_epc, C);
	// scale_and_center(Y_epc);
	
	// bool resid_ePC = true;
	
	// if( resid_ePC ){
			
		// Eigen::MatrixXd U = get_half_hat_matrix(X);

		// std::cerr << "Calculating expression residuals...\n";
		// Y_epc = resid_from_half_hat(Y_epc, U);
	// }
	
	std::cerr << "Estimating latent factors ... \n";
	calc_eGRM_PCs(Q, Q_lambda, Y_epc, n_fac);
	std::cerr << "Done.\n";
	
	Y_epc.resize(0,0);
	
	if( n_fac > 0 ){
		int n_cv = C.cols();
		C.conservativeResize(Eigen::NoChange, n_cv + n_fac);
		for(int ii = 0; ii < n_fac; ii++){
			C.col(n_cv + ii) = Q.col(ii);
			
			std::string col_ii = "factor_";
			col_ii += std::to_string(ii);
			c_data.rows.keep.push_back(col_ii);
		}
		n_covar = C.cols();
	}
	
	Q.resize(0,0);
	

	make_half_hat_matrix(C);

	if( !global_opts::low_mem ){
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = C.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = std::sqrt(g_data.var[i]);
		}
		std::cerr << "Done.\n";
	}
	
	std::cerr << "Calculating expression residuals...\n";
	make_resid_from_half_hat(Y, C);
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y, e_data.stdev);
	
	if( rknorm_r ){
		std::cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y);
		std::cerr << "Re-residualizing transformed residuals ...\n";
		make_resid_from_half_hat(Y, C);
		scale_and_center(Y);
	}
	
	// std::cout << Y_res.format(EigenTSV) << "\n";
	// return 0; X
	
	Y.transposeInPlace();
	
	n_samples = C.rows();
	n_covar = C.cols();
	
	std::string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	std::string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
	BGZF* block_file;
	BGZF* bed_block_file;
	BGZF* long_file;
	
	if( !just_long ){
		
		block_file = bgzf_open(block_file_path.c_str(), "w");
		
		bed_block_file = bgzf_open(bed_block_file_path.c_str(), "w");
		
		write_to_bgzf("#chrom\tstart\tend\tgene\tegene_pval\tn_samples\tn_covar\tresid_sd\tn_cis_variants\n", bed_block_file);
		
	}
	
	bool write_long = (make_long || just_long);
	
	if( write_long ){
		long_file = bgzf_open(long_file_path.c_str(), "w");
		write_to_bgzf("#chrom\tpos\tref\talt\tgene\tbeta\tse\tpval\n", long_file);
	}
	
	int bl = 0;
	
	std::string iter_cerr_suffix = " cis-QTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y.rows() ? n_e : Y.rows()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y.rows() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (C.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
				{
					g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}

				out_mat = (Y.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
			}else{
				
				const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(bm.bcf_s[i], n_g);
				
				out_mat = (Y.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
				// out_mat = (Y.middleRows(bm.bed_s[i], n_e)* g_data.genotypes.middleCols(bm.bcf_s[i], n_g)).transpose().eval();
				
			}
			
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; ++jj, ++jm){
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
				std::vector<double> dist;
				double gene_tss = 0.5*(e_data.start[jj] + e_data.end[jj]);
		
				int pos_s = -1;
				int pos_e = -1;
		
				int v_s = -1;
				int v_e = -1;
		
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, ++im){
				
					std::stringstream long_line;
				
					if( g_data.chr[ii] == e_data.chr[jj] && e_data.start[jj] - g_data.pos[ii] < global_opts::cis_window_bp && g_data.pos[ii] - e_data.end[jj] < global_opts::cis_window_bp ){
						
						if( pos_s < 0 ) pos_s = g_data.pos[ii]; 
						pos_e = g_data.pos[ii];

						if( v_s < 0 ) v_s = ii; 
						v_e = ii;

						double g_yres_crossprod = out_mat(im, jm);
						if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
						
						if( !just_long ) block_line << "\t" << g_yres_crossprod;
						
						double beta, beta_se, zscore, pval_esnp;
						
						if( g_data.var[ii] > 0 ){
							beta = g_yres_crossprod/g_data.var[ii];
							beta_se = std::sqrt( (n_samples - 1)/g_data.var[ii] - beta*beta)/std::sqrt(n_samples - n_covar - 1);
							zscore = beta/beta_se;
							pval_esnp = pf( zscore*zscore, 1.0, n_samples - n_covar - 1, true );
						}else{
							beta = 0;
							beta_se = 0;
							zscore = 0;
							pval_esnp = -1;
						}

						dist.push_back( std::abs(g_data.pos[ii] - gene_tss) );
						pvals.push_back(pval_esnp);

						if( write_long ){

							long_line << 
								clean_chrom(g_data.chr[ii]) << "\t" <<
								g_data.pos[ii] << "\t" <<
								g_data.ref[ii] << "\t" <<
								g_data.alt[ii] << "\t" <<
								e_data.gene_id[jj] << "\t" <<
								y_scale[jj] * e_data.stdev[jj] * beta << "\t" <<
								y_scale[jj] * e_data.stdev[jj] * beta_se  << "\t" <<
								pval_esnp << "\n";
						}
					}
					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						std::string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							std::to_string(pos_s) + "\t" + 
							std::to_string(pos_e) + "\t" + 
							std::to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						std::stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals, dist) << "\t" << 
							C.rows() << "\t" << 
							C.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			std::cerr << "\nERROR: " <<bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	std::cerr << "\n";

	if ( write_long ){
		bgzf_close(long_file);
		//build_tabix_index(long_file_path);
	}

	if( !just_long ){
		bgzf_close(block_file);
		bgzf_close(bed_block_file);
		build_tabix_index(block_file_path, 1);
		build_tabix_index(bed_block_file_path, 1);
	}

	return;
}



void save_fa_covariates(const int& n_fac, genotype_data& g_data, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r)
{
	
	Eigen::MatrixXd Q;
	Eigen::VectorXd Q_lambda;
	
	Eigen::MatrixXd& Y = e_data.data_matrix;
	Eigen::MatrixXd& C = c_data.data_matrix;
		
	// std::cerr << "Reordering trait and covariate matrices ...\n";
	// std::cerr << "Reordering genotypes ...\n";
	// g_data.genotypes = (Tr * g_data.genotypes).eval();
		
	
	double n_traits = Y.cols();
	double n_samples = Y.rows();
	double n_snps = g_data.n_variants;
	double n_covar = C.cols();

	// std::cerr << "Started cis-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);
	
	std::cerr << "Estimating latent factors ... \n";
	calc_eGRM_PCs(Q, Q_lambda, Y, n_fac);
	std::cerr << "Done.\n";
	
	// Y_epc.resize(0,0);
	
	if( n_fac > 0 ){
		int n_cv = C.cols();
		C.conservativeResize(Eigen::NoChange, n_cv + n_fac);
		for(int ii = 0; ii < n_fac; ii++){
			C.col(n_cv + ii) = Q.col(ii);
			
			std::string col_ii = "factor_";
			col_ii += std::to_string(ii);
			c_data.rows.keep.push_back(col_ii);
		}
		n_covar = C.cols();
	}
	
	if( c_data.cols.keep.size() == 0 ){
		c_data.cols.keep = e_data.ids.keep;
	}
	
	std::string out_cov_path = global_opts::out_prefix + ".cov.gz";
	
	c_data.write_table(out_cov_path);
	
}

