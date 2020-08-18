#include "GQT.hpp"

using namespace std;

void run_cis_eQTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	cerr << "Started cis-eQTL analysis ...\n";
	
	if( rknorm_y ){
		cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	cerr << "Scaling expression traits ... \n";
	scale_and_center(Y);

	Eigen::MatrixXd U = get_half_hat_matrix(X);

	if( !global_opts::low_mem ){
		
		cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
		cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = sqrt(g_data.var[i]);
		}
		cerr << "Done.\n";
	}
	
	cerr << "Calculating expression residuals...\n";
	Eigen::MatrixXd Y_res = resid_from_half_hat(Y, U);
	
	cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y_res, e_data.stdev);
	
	if( rknorm_r ){
		cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y_res);
		cerr << "Re-residualizing transformed residuals ...\n";
		Eigen::MatrixXd tmp = resid_from_half_hat(Y_res, U);
		Y_res = tmp;
		scale_and_center(Y_res);
	}
	
	// cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y_res.transposeInPlace();
	
	double n_samples = X.rows();
	double n_covar = X.cols();
	
	string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
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
	
	string iter_cerr_suffix = " cis-eQTL blocks out of " + to_string(bm.size()) + " total";
	cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y_res.rows() ? n_e : Y_res.rows()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y_res.rows() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
				{
					g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}

				out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
			}else{
				
				const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(bm.bcf_s[i], n_g);
				
				out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
				// out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* g_data.genotypes.middleCols(bm.bcf_s[i], n_g)).transpose().eval();
				
			}
			
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; ++jj, ++jm){
				
				stringstream block_line;
				
				vector<double> pvals;
		
				int pos_s = -1;
				int pos_e = -1;
		
				int v_s = -1;
				int v_e = -1;
		
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, ++im){
				
					stringstream long_line;
				
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
							beta_se = sqrt( (n_samples - 1)/g_data.var[ii] - beta*beta)/sqrt(n_samples - n_covar - 1);
							zscore = beta/beta_se;
							pval_esnp = pf( zscore*zscore, 1.0, n_samples - n_covar - 1, true );
						}else{
							beta = 0;
							beta_se = 0;
							zscore = 0;
							pval_esnp = -1;
						}

						pvals.push_back(pval_esnp);

						if( write_long ){

							long_line << 
								clean_chrom(g_data.chr[ii]) << "\t" <<
								g_data.pos[ii] << "\t" <<
								g_data.ref[ii] << "\t" <<
								g_data.alt[ii] << "\t" <<
								e_data.gene_id[jj] << "\t" <<
								e_data.stdev[jj]*beta << "\t" <<
								e_data.stdev[jj]*beta_se  << "\t" <<
								pval_esnp << "\n";
						}
					}
					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							to_string(pos_s) + "\t" + 
							to_string(pos_e) + "\t" + 
							to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals) << "\t" << 
							X.rows() << "\t" << 
							X.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			cerr << "\nERROR: " <<bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	cerr << "\n";

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


void run_cis_eQTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

	cerr << "Starting eigendecomposition of GRM ... ";
	
	Eigen::SelfAdjointEigenSolver <Eigen::SparseMatrix<double>> GRM_eig(GRM);
	
	if (GRM_eig.info() != Eigen::Success){
		cerr << "FATAL ERROR: GRM decomposition failed!\n";
		abort();
	}
	cerr << "Done!\n";
	
	Eigen::VectorXd GRM_lambda = GRM_eig.eigenvalues();
	Eigen::SparseMatrix<double> L = GRM_eig.eigenvectors().sparseView(1.00, 5e-3);
	L.makeCompressed();

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	cerr << "Started cis-eQTL analysis ...\n";
	
	if( rknorm_y ){
		cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	cerr << "Scaling expression traits ... \n";
	scale_and_center(Y);
	
	cerr << "Rotating expression and covariates ... ";
	Y = (L.transpose() * Y).eval();
	X = (L.transpose() * X).eval();
	cerr << "Done!\n";

	Eigen::MatrixXd U = get_half_hat_matrix(X);
	
	if( !global_opts::low_mem ){
		
		cerr << "Rotating genotypes ... ";
		g_data.genotypes = (L.transpose()*g_data.genotypes).eval();
		g_data.genotypes.makeCompressed();
		cerr << "Done!\n";
		
		cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
		cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = sqrt(g_data.var[i]);
		}
		cerr << "Done.\n";
	}
	
	/* 
	cerr << "Calculating expression residuals...\n";
	Eigen::MatrixXd Y_res = resid_from_half_hat(Y, U);
	
	cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y_res, e_data.stdev);
	
	if( rknorm_r ){
		// cerr << "Rank-normalizing expression residuals ...\n";
		// rank_normalize(Y_res);
		// cerr << "Re-residualizing transformed residuals ...\n";
		// Eigen::MatrixXd tmp = resid_from_half_hat(Y_res, U);
		// Y_res = tmp;
		// scale_and_center(Y_res);
		cerr << "FATAL ERROR: Rank-norm resid with LMM is not currently supported.\n";
		abort();
	}
	
	// cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y_res.transposeInPlace();
	*/
	
	double n_samples = X.rows();
	double n_covar = X.cols();
	
	string block_file_path = global_opts::out_prefix + "." + "cis_sumstats" + ".txt.gz";
	string bed_block_file_path = global_opts::out_prefix + "." + "cis_gene_table" + ".txt.gz";
	string long_file_path = global_opts::out_prefix + "." + "cis_long_table" + ".txt.gz";
	
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
	
	string iter_cerr_suffix = " cis-eQTL blocks out of " + to_string(bm.size()) + " total";
	cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	e_data.stdev.resize( Y.cols() );
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < Y.cols() ? n_e : Y.cols()-1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < Y.cols() && bm.bcf_s[i] < g_data.n_variants ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
						
				g_data.genotypes = (L.transpose()*g_data.genotypes).eval();
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
				{
					g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}

				// out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
			}
			
			int s_slice = global_opts::low_mem ? 0 : bm.bcf_s[i];
			int n_slice = global_opts::low_mem ? g_data.genotypes.cols() : n_g;
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_slice, n_slice);
			
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; ++jj, ++jm){
				
				stringstream block_line;
				
				vector<double> pvals;
		
				int pos_s = -1;
				int pos_e = -1;
		
				int v_s = -1;
				int v_e = -1;
				
				int idx_s = -1;
				int idx_n = 0;
				
				for( int ii = bm.bcf_s[i], im = 0; ii < bm.bcf_s[i] + n_g; ++ii, im++ ){
					if( g_data.chr[ii] == e_data.chr[jj] && e_data.start[jj] - g_data.pos[ii] < global_opts::cis_window_bp && g_data.pos[ii] - e_data.end[jj] < global_opts::cis_window_bp ){
						if( idx_s < 0 ) idx_s = im;
						
						if( pos_s < 0 ) pos_s = g_data.pos[ii]; 
						pos_e = g_data.pos[ii];

						if( v_s < 0 ) v_s = ii; 
						v_e = ii;

						idx_n++;
					}else{
						if( idx_s >= 0 ){
							break;
						}
					}
				}
		
				const Eigen::SparseMatrix<double>& G_slice = G.middleCols(idx_s, idx_n);
				
								
				LMM_fitter fit(X, Y.col(jj), GRM_lambda);
				fit.fit_REML();
				
				cout << e_data.gene_id[jj] << "\t" << fit.sigma2 << "\t" << fit.phi << "\n";

				
				Eigen::MatrixXd XtDXi = (X.transpose() * fit.Vi * X).inverse();
				Eigen::MatrixXd XtDG = X.transpose() * fit.Vi * G_slice;
				Eigen::MatrixXd XtDXi_XtDG = XtDXi * XtDG;
				
				Eigen::VectorXd y_res = Y.col(jj) - X * XtDXi * X.transpose() * fit.Vi * Y.col(jj);
				
				double SSR_0 = y_res.dot(fit.Vi * Y.col(jj))/fit.sigma2;
				
				Eigen::VectorXd U_vec = G_slice.transpose() * fit.Vi * y_res/sqrt(fit.sigma2);
				
				Eigen::VectorXd diagV(U_vec.size());
				
				for (int si = 0; si < U_vec.size(); si++){
					diagV(si) = G_slice.col(si).dot(fit.Vi * G_slice.col(si)) - XtDG.col(si).dot(XtDXi_XtDG.col(si));	
				}
				
				e_data.stdev[jj] = sqrt(fit.sigma2);
				
				for( int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					stringstream long_line;
				
					double g_yres_crossprod = U_vec(si);
					if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
					
					if( !just_long ) block_line << "\t" << g_yres_crossprod;
					
					double beta, beta_se, zscore, pval_esnp;
					
					if( g_data.var[ii] > 0 ){
						beta = g_yres_crossprod/diagV(si);
						beta_se = sqrt( SSR_0/diagV(si) - beta*beta)/sqrt(n_samples - n_covar - 1);
						zscore = beta/beta_se;
						pval_esnp = pf( zscore*zscore, 1.0, n_samples - n_covar - 1, true );
					}else{
						beta = 0;
						beta_se = 0;
						zscore = 0;
						pval_esnp = -1;
					}

					pvals.push_back(pval_esnp);

					if( write_long ){

						long_line << 
							clean_chrom(g_data.chr[ii]) << "\t" <<
							g_data.pos[ii] << "\t" <<
							g_data.ref[ii] << "\t" <<
							g_data.alt[ii] << "\t" <<
							e_data.gene_id[jj] << "\t" <<
							e_data.stdev[jj]*beta << "\t" <<
							e_data.stdev[jj]*beta_se  << "\t" <<
							pval_esnp << "\n";
					}

					if( write_long ) write_to_bgzf(long_line.str().c_str(), long_file);
				}
				
				if( !just_long ){
					
					if( pos_s > 0 && pos_e > 0 ){
						
						block_line << "\n";
						
						string block_line_left = 
							clean_chrom(e_data.chr[jj]) + "\t" + 
							to_string(pos_s) + "\t" + 
							to_string(pos_e) + "\t" + 
							to_string(v_s);
						
						write_to_bgzf(block_line_left, block_file);
						write_to_bgzf(block_line.str().c_str(), block_file);
						
					}
					
					if( pvals.size() > 0 ){
						
						stringstream bed_block_line;
						
						bed_block_line <<
							clean_chrom(e_data.chr[jj]) << "\t" << 
							e_data.start[jj] << "\t" << 
							e_data.end[jj] << "\t" << 
							e_data.gene_id[jj] << "\t" << 
							ACAT_non_missing(pvals) << "\t" << 
							X.rows() << "\t" << 
							X.cols() << "\t" << 
							e_data.stdev[jj] << "\t" <<
							v_e - v_s + 1 << "\n";
						
						write_to_bgzf(bed_block_line.str().c_str(), bed_block_file);
					}
				}
			}
			
		}else{
			cerr << "\nERROR: " << bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
		
		bl++;
	}
	cerr << "\n";

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

