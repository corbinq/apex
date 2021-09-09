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


#include "transMapping.hpp"


void parse_cis_signal_data(const std::string& fn, const std::string& region, Eigen::MatrixXd& X, std::vector<std::string>& samples, std::vector<std::string>& genes){
	data_parser dp;
	dp.add_header(samples, 5);
	dp.add_field(genes, 3);
	dp.add_matrix(X, true, 5);
	dp.parse_file(fn, region);
	return;
}

Eigen::VectorXd getVectorXd(const std::vector<double>& v, const int& s, const int& n){
	Eigen::VectorXd ev(n);
	for(int i = 0, ii = s; i < n; i++, ii++){
		ev(i) = v[ii];
	}
	return ev;
}


double u_stat_pval(const double& u_stat, const double& m, const double& n){
	double f_stat = u_stat * u_stat;
	f_stat = (n - m - 1)*f_stat/(n - 1 - f_stat);
	if( f_stat < 0.00 || std::isnan(f_stat) ){
		std::cerr << "Warning: F statistic = " << f_stat << "\n";
		return -1.00;
	}
	return pf(f_stat, 1, n - m - 1, true);
}

double usq_stat_pval(const double& usq_stat, const double& m, const double& n){
	double f_stat = (n - m - 1)*usq_stat/(n - 1 - usq_stat);
	if( f_stat < 0.00 || std::isnan(f_stat) ){
		std::cerr << "Warning: F statistic = " << f_stat << "\n";
		return -1.00;
	}
	return pf(f_stat, 1, n - m - 1, true);
}

void run_trans_QTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const int& chunk_size, const std::string& cis_fp, const std::string& cis_fp_region)
{

	struct Zsn {int s; int n; };

	Eigen::MatrixXd Z;
	std::vector<std::string> Z_genes;
	std::vector<std::string> Z_samples;
	std::vector<Zsn> Z_map;
	bool adj_cis = false;
	bool adj_cis_sloppy = global_opts::sloppy_covar;

	if( cis_fp != "" ){
		adj_cis = true;
		
		parse_cis_signal_data(cis_fp, "", Z, Z_samples, Z_genes);
		
		std::unordered_map<std::string, Zsn> Z_map0;
		
		std::string last_gene = Z_genes[0];
		int ni = 0;
		for(int i = 0; i < Z_genes.size(); i++){
			std::string& gene_i = Z_genes[i];
			if( gene_i == last_gene ){
				ni++;
			}else{
				Z_map0[last_gene] = {i - ni, ni};
				ni = 1;
				last_gene = gene_i;
			}
			if( i + 1 == Z_genes.size() ){
				Z_map0[gene_i] = {i - ni, ni};
			}
		}
		
		for( const std::string& gene_i : e_data.gene_id ){
			if( Z_map0.find(gene_i) == Z_map0.end() ){
				Z_map.push_back({-1,0});
			}else{
				// std::cerr << "\n\nMatch!\n\n";
				Z_map.push_back(Z_map0[gene_i]);
			}
		}
	}

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &U = c_data.data_matrix;
	
	make_half_hat_matrix(U);
	
	std::cerr << "Started trans-QTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	std::vector<double> y_scale;
	scale_and_center(Y, y_scale);

	if( !global_opts::low_mem ){
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";
		for( int i = 0; i < UtG.cols(); i++)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
		}
		std::cerr << "Done.\n";
	}
	
	std::cerr << "Calculating expression residuals...\n";
	make_resid_from_half_hat(Y, U);
	
	if( rknorm_r ){
		std::cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y);
		std::cerr << "Re-residualizing transformed residuals ...\n";
		make_resid_from_half_hat(Y, U);
		scale_and_center(Y);
	}
	
	if( adj_cis ){
		std::cerr << "Residualizing cis signals...\n";
		make_resid_from_half_hat(Z, U);
		std::cerr << "Re-residualizing expression ...\n";
		for( int i = 0; i < e_data.gene_id.size(); i++){
			const int& s_i = Z_map[i].s;
			const int& n_i = Z_map[i].n;
			if( n_i > 0 ){
				Eigen::MatrixXd Z_i = get_half_hat_matrix(Z.middleCols(s_i, n_i));
				Z.middleCols(s_i, n_i) = Z_i;
				
				Eigen::VectorXd y_res_z = resid_vec_from_half_hat(Y.col(i), Z_i);
				Y.col(i) = y_res_z;
			}
		}
	}
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y, e_data.stdev);
	
	Y.transposeInPlace();
	
	double n_samples = U.rows();
	double n_covar = U.cols();
	
	int n_genes = Y.rows();
	
	std::vector<double> gene_max_val(n_genes, 0.0);
	std::vector<int> gene_max_idx(n_genes, 0);
	
	std::string block_file_path = global_opts::out_prefix + "." + "trans_sumstats" + ".txt.gz";
	std::string gene_file_path = global_opts::out_prefix + "." + "trans_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "trans_long_table" + ".txt.gz";
	
	BGZF* block_file;
	BGZF* gene_file;
	BGZF* long_file;
	
	gene_file = bgzf_open(gene_file_path.c_str(), "w");
	write_to_bgzf("#chrom\tpos\tref\talt\tgene_chrom\tgene_id\tbeta\tse\tpval\n", gene_file);
	
	long_file = bgzf_open(long_file_path.c_str(), "w");
	write_to_bgzf("#chrom\tpos\tref\talt\tgene_chrom\tgene_id\tbeta\tse\tpval\n", long_file);
	
	int bl = 0;
	
	int n_blocks = ceil(g_data.n_variants/chunk_size);
	
	std::string iter_cerr_suffix = " genotype blocks out of " + std::to_string(n_blocks) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	double F_crit = qf(global_opts::LM_ALPHA, 1, n_samples - n_covar - 1, true);
	double P_crit = std::sqrt(F_crit * (n_samples - 1)/( F_crit + n_samples - n_covar - 1));
	
	for( ; ; bl++ ){
		
		int s_g = bl * chunk_size;
		int n_g = chunk_size;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		n_g = n_g < g_data.n_variants - s_g ? n_g : g_data.n_variants - s_g;
		
		if( s_g >= g_data.n_variants || n_g <= 0 ){
			break; 
		}
		
		if( n_g > 0 && s_g < g_data.n_variants ){
			
			Eigen::MatrixXd StdScore;
			Eigen::MatrixXd Denom;
			Eigen::VectorXd dV;
			
			int s_G, n_G;
			if( global_opts::low_mem ){
				
				g_data.read_genotypes(sr, hdr, s_g, n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval(); 
				
				for( int si = s_g, ii = 0; si < s_g + n_g; si++, ii++)
				{
					g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);
					//cout << g_data.var[i] << "\n";
				}
				s_G = 0;
				n_G = g_data.genotypes.cols();
			}else{
				s_G = s_g;
				n_G = n_g;
			}
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_G, n_G);
			
			dV = getVectorXd(g_data.var, s_g, n_g);
			
			if( (!adj_cis) || adj_cis_sloppy ){
				StdScore = dV.cwiseSqrt().asDiagonal().inverse() * (Y * G).transpose().eval();
			}else{
				StdScore = (Y * G).transpose();
				Eigen::MatrixXd GtZ = G.transpose() * Z;
				
				Denom.resize( StdScore.rows(), StdScore.cols() );
				
				for( int gi = 0; gi < e_data.gene_id.size(); gi++){
					const int& s_i = Z_map[gi].s;
					const int& n_i = Z_map[gi].n;
					if( n_i > 0 ){
						Eigen::VectorXd dV_adj = dV - GtZ.middleCols(s_i, n_i).cwiseAbs2().rowwise().sum();
						Denom.col(gi) = dV - GtZ.middleCols(s_i, n_i).cwiseAbs2().rowwise().sum();
					}else{
						Denom.col(gi) = dV;
					}
				}
				StdScore = StdScore.cwiseQuotient(Denom.cwiseSqrt());
			}
			
			std::stringstream long_line;
			
			for(int j = 0; j < StdScore.cols(); j++){
				
				double n_covar_j = n_covar;
				
				if( adj_cis && !adj_cis_sloppy ){
					n_covar_j += Z_map[j].n;
				}
				
				for(int i = 0; i < StdScore.rows(); i++){
					const auto& val = StdScore(i,j);
					
					if( std::abs(val) > std::abs(gene_max_val[j]) ){
						gene_max_val[j] = val;
						gene_max_idx[j] = s_g + i;
					}
					
					if( val > P_crit || val < - P_crit ){
						
						double V;
						if( adj_cis && !adj_cis_sloppy ){
							V = Denom(i,j);
							if( V/dV(i) < 0.001 || V <= 0 || std::isnan(V) ){
								continue;
							}
						}else{
							V = dV(i);
						}
							
						const double& scale = e_data.stdev[j];
						double U = val * std::sqrt(V);
						double beta = U/V;
						double beta_se = std::sqrt( ((n_samples - 1)/V- beta*beta)/(n_samples - n_covar_j - 1) );
						

						double pval_esnp = u_stat_pval(val, n_covar_j, n_samples);
						
						int ii = s_g + i;
						if( g_data.flipped[ii] ) beta *= (-1.00); 

						long_line << 
							clean_chrom(g_data.chr[ii]) << "\t" <<
							g_data.pos[ii] << "\t" <<
							g_data.ref[ii] << "\t" <<
							g_data.alt[ii] << "\t" <<
							clean_chrom(e_data.chr[j]) << "\t" <<
							e_data.gene_id[j] << "\t" <<
							y_scale[j]*scale*beta << "\t" <<
							y_scale[j]*scale*beta_se  << "\t" <<
							pval_esnp << "\n";
					}
				}
			}
			
			write_to_bgzf(long_line.str().c_str(), long_file);
			
		}else{
			std::cerr << "\nERROR: " <<bl << ", " << s_g << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
	}
	std::cerr << "\n";

	bgzf_close(long_file);
	
	for( int j = 0; j < n_genes; j++){
		
		const int& ii = gene_max_idx[j];
		const int& val = gene_max_val[j];

		const double& V = g_data.var[ii];
		const double& scale = e_data.stdev[j];
		double U = val * std::sqrt(V);
		double beta = U/V;
		double beta_se = std::sqrt( ((n_samples - 1)/V- beta*beta)/(n_samples - n_covar - 1) );

		if( g_data.flipped[ii] ) beta *= (-1.00);

		std::stringstream gene_line;
		gene_line << 
			clean_chrom(g_data.chr[ii]) << "\t" <<
			g_data.pos[ii] << "\t" <<
			g_data.ref[ii] << "\t" <<
			g_data.alt[ii] << "\t" <<
			clean_chrom(e_data.chr[j]) << "\t" <<
			e_data.gene_id[j] << "\t" <<
			scale*beta << "\t" <<
			scale*beta_se  << "\t" <<
			u_stat_pval(val, n_covar, n_samples) << "\n";
		
		write_to_bgzf(gene_line.str().c_str(), gene_file);	
	}
	
	
	bgzf_close(gene_file);
	
	return;
}


void run_trans_QTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr, genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& GRM_lambda, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const int& chunk_size, const std::string& theta_path, const std::string& anchor_path)
{
	
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
	
	std::vector<double> gene_max_val((int) n_traits, 0.0);
	std::vector<int> gene_max_idx((int) n_traits, 0);
	
	std::string block_file_path = global_opts::out_prefix + "." + "trans_sumstats" + ".txt.gz";
	std::string gene_file_path = global_opts::out_prefix + "." + "trans_gene_table" + ".txt.gz";
	std::string long_file_path = global_opts::out_prefix + "." + "trans_long_table" + ".txt.gz";
	
	BGZF* block_file;
	BGZF* gene_file;
	BGZF* long_file;
	
	gene_file = bgzf_open(gene_file_path.c_str(), "w");
	write_to_bgzf("#chrom\tpos\tref\talt\tgene_chrom\tgene_id\tbeta\tse\tpval\n", gene_file);
	
	long_file = bgzf_open(long_file_path.c_str(), "w");
	write_to_bgzf("#chrom\tpos\tref\talt\tgene_chrom\tgene_id\tbeta\tse\tpval\n", long_file);
	
	
	std::vector<double> phi_v(Y.cols());
	std::vector<double> hsq_v(Y.cols());
	std::vector<double> sigma_v(Y.cols());
	std::vector<double> SSR_v(Y.cols());
	
	e_data.stdev.resize(Y.cols());
	
	theta_data t_data;
	bool use_theta = ( theta_path != "" );
	
	if( use_theta ){
		std::cerr << "Set null model for ";
		t_data.open(theta_path);
	}else{
		std::cerr << "Fit null model for ";
	}
	
	X = (L.transpose() * C).eval();
	
	std::string iter_cerr_suffix = " traits out of " + std::to_string(Y.cols()) + " total";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	int last_j = 0;
	for(int j = 0; j < (int) n_traits; j++ ){
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
			
			// Y now stores the rotated residuals. 
			Y.col(j) = (Vi * y_res/std::sqrt(sigma2)).eval();
		}
		
		phi_v[j] = phi;
		hsq_v[j] = hsq;
		sigma_v[j] = sigma2;
		
		e_data.stdev[j] = std::sqrt(sigma2);
		
		thinned_iter_cerr(last_j, j+1, iter_cerr_suffix, 20);
	}
	print_iter_cerr(last_j, Y.cols(), iter_cerr_suffix);
	
	if( !e_data.is_residual ){
		Y = (L * Y).eval();
	}
	
	Eigen::MatrixXd VBeta = getVBeta(hsq_v, phi_v, hsq_vals.size());
	
	std::cerr << "\n";
	
	int bl = 0;
	
	int n_blocks = ceil(g_data.n_variants/chunk_size);
	
	iter_cerr_suffix = " genotype blocks out of " + std::to_string(n_blocks) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	double F_crit = qf(global_opts::LM_ALPHA, 1, n_samples - n_covar - 1, true);
	double P_crit = std::sqrt(F_crit * (n_samples - 1)/( F_crit + n_samples - n_covar - 1));
	double Psq_crit = P_crit*P_crit;
	

	for( ; ; bl++ ){
		
		int s_g = bl * chunk_size;
		int n_g = chunk_size;
		n_g = n_g < g_data.n_variants ? n_g : g_data.n_variants-1;
		n_g = n_g < g_data.n_variants - s_g ? n_g : g_data.n_variants - s_g;
		
		if( s_g >= g_data.n_variants || n_g <= 0 ){
			break; 
		}
		
		if( n_g > 0 && s_g < g_data.n_variants ){
			
			Eigen::MatrixXd StdScore;
			Eigen::VectorXd dV;
			
			if( global_opts::low_mem ){
				std::cerr << "Low mem trans LMM not supported\n.";
				abort();
			}
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_g, n_g);
			
			Eigen::MatrixXd U_b = G.transpose() * Y;
			Eigen::MatrixXd V_b = V_mat.middleRows(s_g, n_g) * VBeta;
			
			Eigen::VectorXd scale_vec(U_b.cols());
			for(int i = 0; i < U_b.cols(); i++){
				scale_vec(i) = (n_samples - 1)/SSR_v[i];
			}
			
			Eigen::MatrixXd StdScore2 = (U_b.cwiseAbs2() * scale_vec.asDiagonal()).cwiseQuotient(V_b);
			
			std::stringstream long_line;
			
			for(int j = 0; j < U_b.cols(); j++){
				for(int i = 0; i < U_b.rows(); i++){
					
					const double& val2 = StdScore2(i,j);
					
					if( val2 > gene_max_val[j] ){
						gene_max_val[j] = val2;
						gene_max_idx[j] = s_g + i;
					}
					
					if( val2 > Psq_crit ){
						const double& V = V_b(i,j);
						const double& scale = e_data.stdev[j];
						const double& U = U_b(i,j);
						double beta = U/V;
						double beta_se = std::sqrt( (SSR_v[j]/V- beta*beta)/(n_samples - n_covar - 1) );
						
						double pval_esnp = usq_stat_pval(val2, n_covar, n_samples);
						
						int ii = s_g + i;

						if( g_data.flipped[ii] ) beta *= (-1.00);

						long_line << 
							clean_chrom(g_data.chr[ii]) << "\t" <<
							g_data.pos[ii] << "\t" <<
							g_data.ref[ii] << "\t" <<
							g_data.alt[ii] << "\t" <<
							clean_chrom(e_data.chr[j]) << "\t" <<
							e_data.gene_id[j] << "\t" <<
							y_scale[j]*scale*beta << "\t" <<
							y_scale[j]*scale*beta_se  << "\t" <<
							pval_esnp << "\n";
					}
				}
			}
			
			write_to_bgzf(long_line.str().c_str(), long_file);
		}else{
			std::cerr << "\nERROR: " <<bl << ", " << s_g << ", " << n_g << "\n"; 
			abort();
		}
		
		print_iter_cerr(bl, bl+1, iter_cerr_suffix);
	}
	std::cerr << "\n";

	bgzf_close(long_file);
	bgzf_close(gene_file);

	return;
}



