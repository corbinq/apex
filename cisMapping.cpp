#include "GQT.hpp"


void run_cis_eQTL_analysis(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	std::cerr << "Started cis-eQTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	scale_and_center(Y);

	Eigen::MatrixXd U = get_half_hat_matrix(X);

	if( !global_opts::low_mem ){
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
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
	Eigen::MatrixXd Y_res = resid_from_half_hat(Y, U);
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y_res, e_data.stdev);
	
	if( rknorm_r ){
		std::cerr << "Rank-normalizing expression residuals ...\n";
		rank_normalize(Y_res);
		std::cerr << "Re-residualizing transformed residuals ...\n";
		Eigen::MatrixXd tmp = resid_from_half_hat(Y_res, U);
		Y_res = tmp;
		scale_and_center(Y_res);
	}
	
	// std::cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y_res.transposeInPlace();
	
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
	
	std::string iter_cerr_suffix = " cis-eQTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
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
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
		
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

Eigen::VectorXd calc_Py(const Eigen::MatrixXd& U, const Eigen::MatrixXd& C, const Eigen::MatrixXd& Q, const Eigen::VectorXd& Psi, const Eigen::VectorXd& Y,const Eigen::VectorXd& Sq, const Eigen::VectorXd& Sc ){
	Eigen::VectorXd Py = 
		Y - Q * Psi.asDiagonal() * Sq - 
		(
			C - Q * Psi.asDiagonal() * U.transpose()
		) * (
			Eigen::MatrixXd::Identity(Sq.size(),Sq.size()) - 
			U * Psi.asDiagonal() * U.transpose()
		).ldlt().solve(
			Sc - U * Psi.asDiagonal() * Sq
		);
	return Py;
}

typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalXd;

DiagonalXd calc_Psi(const double& delta, const Eigen::VectorXd& lambdas){
	Eigen::VectorXd numer = delta * lambdas;
	Eigen::VectorXd denom = delta * lambdas;
	denom.array() += 1.00;
	return (numer.cwiseQuotient(denom)).asDiagonal().inverse();	
}

void run_cis_eQTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{
	
	double delta_thresh = 0.02;
	double drop0_thresh = 5e-3;

	Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic> Tr(GRM.rows()); 

	std::vector<int> unrelateds;
	
	int i = 0;
	for( int j = 0; j < GRM.cols(); j++ ){
		if( i < relateds.size() ){
			if( relateds[i] == j ){
				i++;
			}else{
				unrelateds.push_back(j);
			}
		}else{
			unrelateds.push_back(j);
		}
	}
	i = 0;
	for( const int& ii : relateds ){
		Tr.indices()[ii] = i;
		i++;
	}
	for( const int& ii : unrelateds ){
		Tr.indices()[ii] = i;
		i++;
	}

	int nrel = relateds.size();
	int nind = unrelateds.size();
	
	std::cerr << "Found "<< nrel << " related individuals ... \n";
	
	std::cerr << "Reordering related GRM blocks ... ";
	
	// A better way to do this would be  
	//    grm = grm.twistedBy(Tr);
	// Oddly, this fails on my system. 
	
	GRM = (Tr * GRM).eval();
	GRM = (GRM * Tr.transpose()).eval();
	
	GRM.makeCompressed();

	std::cerr << "Done.\n";

	Eigen::SparseMatrix<double> GRM_rel = GRM.topLeftCorner(nrel,nrel);
	GRM_rel.makeCompressed();

	std::cerr << "GRM eigendecomposition ... ";
	
	Eigen::SelfAdjointEigenSolver <Eigen::SparseMatrix<double>> GRM_eig(GRM_rel);
	
	if (GRM_eig.info() != Eigen::Success){
		std::cerr << "FATAL ERROR: GRM decomposition failed!\n";
		abort();
	}
	std::cerr << "Done.\n";
	
	Eigen::VectorXd GRM_lambda = GRM_eig.eigenvalues();
	
	GRM_lambda.conservativeResize(nrel + nind);
	for(int i = nrel; i < nrel + nind; i++){
		GRM_lambda(i) = 1;
	}
	
	
	Eigen::SparseMatrix<double> L = GRM_eig.eigenvectors().sparseView(drop0_thresh, 1.0 - std::numeric_limits<double>::epsilon());
	
	L.conservativeResize(nrel + nind,nrel + nind);
	
	for(int i = nrel; i < nrel + nind; i++){
		L.coeffRef(i,i) = 1.00;
	}
	
	L.makeCompressed();
	
	
	/*
	int ii = 0;
	int ni = g_data.ids.keep.size();
	for( const auto & cn : g_data.ids.keep ){
		std::cout << cn;
		ii++;
		if( ii < ni ){
			std::cout << "\t";
		}
	}
	std::cout << "\n";
	
	std::cout << GRM_lambda(0);
	for( int i = 1; i < GRM_lambda.size(); i++ ){
			std::cout << "\t" << GRM_lambda(i);
	}
	std::cout << "\n";
	
	std::cout << Eigen::MatrixXd(L).format(EigenTSV) << "\n";
	abort();
	*/
	
	std::vector<int> kp;
	for( int i = 0; i < GRM_lambda.size(); i++ ){
		if( std::abs( GRM_lambda(i) - 1 ) > delta_thresh ){
			kp.push_back(i);
		}
	}
	
	// Permutation matrix to extract selected eigenvectors.
	using td = Eigen::Triplet<double>;
	Eigen::SparseMatrix<double> PC_sel(GRM_lambda.size(),kp.size());
	std::vector<td> PC_trips;
	
	// Selected eigenvalues.
	Eigen::VectorXd Q_lambda(kp.size());
	
	i = 0;
	for( const int& k : kp){
		PC_trips.push_back(td(k,i,1.00));
		Q_lambda(i) =  GRM_lambda(k);
		i++;
	}
	PC_sel.setFromTriplets(PC_trips.begin(), PC_trips.end());
	
	std::cerr << "Selected " << kp.size() << " eigenvectors.\n";
	

	Eigen::SparseMatrix<double> Q = (L * PC_sel).eval();
	Q.makeCompressed();

	std::cerr << "Reordering trait and covariate matrices ...\n";
	
	Eigen::MatrixXd Y = (Tr * e_data.data_matrix).eval();
	Eigen::MatrixXd C = (Tr * c_data.data_matrix).eval();
	
	// Eigen::MatrixXd C = get_half_hat_matrix(C_raw);
	
	int n_traits = Y.cols();
	int n_samples = Y.rows();
	int n_snps = g_data.n_variants;
	int n_covar = C.cols();

	std::cerr << "Started cis-eQTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	scale_and_center(Y);
	
	std::cerr << "Calculating partial rotations ...\n";
	Eigen::MatrixXd QtC = Q.transpose() * C;	
	Eigen::MatrixXd QtY = Q.transpose() * Y;
	Eigen::MatrixXd CtY = C.transpose() * Y;
	Eigen::MatrixXd CtC = C.transpose() * C;
	
	std::cerr << "Rotating expression and covariates ... ";
	Eigen::MatrixXd Yr = L.transpose() * Y;
	Eigen::MatrixXd Cr = L.transpose() * C;
	std::cerr << "Done.\n";

	std::vector<double> GtG_diag(n_snps);
	
	std::vector<double> hsq_vals{0.0, 0.5, 1.0};
	int n_hsq = hsq_vals.size();
	
	std::vector<Eigen::MatrixXd> M_list;
	std::vector<std::vector<double>> V_list(n_hsq);

	for( const double& hsq : hsq_vals ){
		M_list.push_back(
			(CtC - (QtC.transpose() * calc_Psi(hsq, Q_lambda) * QtC)).inverse()
		);
	}
	

	if( !global_opts::low_mem ){
		std::cerr << "Genotype Q-GRM product ... \n";
		Eigen::SparseMatrix<double> QtG = (Q.transpose() * g_data.genotypes).eval();
	
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd CtG = C.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";

		for( int i = 0; i < n_snps; i++)
		{
			GtG_diag[i] = g_data.genotypes.col(i).squaredNorm();
			g_data.var[i] = GtG_diag[i] - CtG.col(i).squaredNorm();
			//SD_vec[i] = std::sqrt(g_data.var[i]);
		}
		std::cerr << "Done.\n";
		
		
		
		for(int j = 0; j < n_hsq; j++){
			std::cerr << "Calculating genotype variance basis (" << j+1 << "/" << n_hsq << ") ... \n";

			Eigen::SparseMatrix<double> Psi_a = calc_Psi(hsq_vals[j], Q_lambda) * QtG;
			Eigen::MatrixXd Res = CtG - QtC.transpose() * Psi_a;
			//Eigen::Matrix MRes = M_list[j] * Res;
			V_list[j].resize(n_snps);
			
			for( int i = 0; i < n_snps; i++)
			{
				V_list[j][i] = (
					GtG_diag[i] - QtG.col(i).dot(Psi_a.col(i)) - Res.col(i).dot(M_list[j] * Res.col(i))
				);
			}
		}
		
		std::cerr << "Done.\n";
	}
	
	/* 
	std::cerr << "Calculating expression residuals...\n";
	Eigen::MatrixXd Y_res = resid_from_half_hat(Y, U);
	
	std::cerr << "Scaling expression residuals ...\n";
	scale_and_center(Y_res, e_data.stdev);
	
	if( rknorm_r ){
		// std::cerr << "Rank-normalizing expression residuals ...\n";
		// rank_normalize(Y_res);
		// std::cerr << "Re-residualizing transformed residuals ...\n";
		// Eigen::MatrixXd tmp = resid_from_half_hat(Y_res, U);
		// Y_res = tmp;
		// scale_and_center(Y_res);
		std::cerr << "FATAL ERROR: Rank-norm resid with LMM is not currently supported.\n";
		abort();
	}
	
	// std::cout << Y_res.format(EigenTSV) << "\n";
	// return 0;
	
	Y_res.transposeInPlace();
	*/
	
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
	
	std::string iter_cerr_suffix = " cis-eQTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
	print_iter_cerr(1, 0, iter_cerr_suffix);
	
	e_data.stdev.resize( n_traits );
	
	for( int i = 0; i < bm.size(); i++ ){
		
		int n_e = bm.bed_e[i] - bm.bed_s[i] + 1;
		n_e = n_e < n_traits ? n_e : n_traits - 1;
		
		int n_g = bm.bcf_e[i] - bm.bcf_s[i] + 1;
		n_g = n_g < n_snps ? n_g : n_snps - 1;
		
		if( n_g > 0 && n_e > 0 && bm.bed_s[i] < n_traits && bm.bcf_s[i] < n_snps ){
			
			Eigen::MatrixXd out_mat;
			
			if( global_opts::low_mem ){
						
				// g_data.genotypes = (L.transpose()*g_data.genotypes).eval();
				
				g_data.read_genotypes(sr, hdr, bm.bcf_s[i], n_g );
				
				Eigen::SparseMatrix<double>& G = g_data.genotypes;
				
				Eigen::VectorXd CtG_block_sqnm = (C.transpose() * G).colwise().squaredNorm().eval(); // G.middleCols(bm.bcf_s[i], n_g);
				
				for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; si++, ii++)
				{
					
					GtG_diag[si] = G.col(ii).squaredNorm();
					g_data.var[si] = GtG_diag[si] - CtG_block_sqnm(ii);
				}

				// out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();
				
			}
			
			int s_slice = global_opts::low_mem ? 0 : bm.bcf_s[i];
			int n_slice = global_opts::low_mem ? g_data.genotypes.cols() : n_g;
			
			const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(s_slice, n_slice);
			
			
			for( int jj = bm.bed_s[i], jm = 0; jj < bm.bed_s[i] + n_e; ++jj, ++jm){
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
		
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
				
								
				LMM_fitter fit(Cr, Yr.col(jj), GRM_lambda);
				fit.fit_REML();
				
				
				double tau2 = fit.phi*fit.sigma2;
				double hsq = tau2 / (tau2 + fit.sigma2);
				
				std::cout << e_data.gene_id[jj] << "\t" << fit.sigma2 << "\t" << fit.phi << "\t" << tau2 << "\t" << hsq << "\n";

				
				Eigen::MatrixXd CtDCi = (Cr.transpose() * fit.Vi * Cr).inverse();
				Eigen::MatrixXd CtDG = Cr.transpose() * fit.Vi * G_slice;
				Eigen::MatrixXd CtDCi_CtDG = CtDCi * CtDG;
				
				Eigen::VectorXd y_res = Yr.col(jj) - Cr * CtDCi * Cr.transpose() * fit.Vi * Yr.col(jj);
				
				//Eigen::VectorXd y_res = calc_Py(U, C, Q, Psi, Y.col(jj), QtY.col(jj), CtY.col(jj));
				
				
				double SSR_0 = y_res.dot(fit.Vi * Yr.col(jj))/fit.sigma2;
				
				Eigen::VectorXd U_vec = G_slice.transpose() * fit.Vi * y_res/std::sqrt(fit.sigma2);
				
				Eigen::VectorXd diagV(U_vec.size());
				
				for (int si = 0; si < U_vec.size(); si++){
					diagV(si) = G_slice.col(si).dot(fit.Vi * G_slice.col(si)) - CtDG.col(si).dot(CtDCi_CtDG.col(si));	
				}
				
				e_data.stdev[jj] = std::sqrt(fit.sigma2);
				
				for( int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;
				
					double g_yres_crossprod = U_vec(si);
					if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
					
					if( !just_long ) block_line << "\t" << g_yres_crossprod;
					
					double beta, beta_se, zscore, pval_esnp;
					
					if( g_data.var[ii] > 0 ){
						beta = g_yres_crossprod/diagV(si);
						beta_se = std::sqrt( SSR_0/diagV(si) - beta*beta)/std::sqrt(n_samples - n_covar - 1);
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
							ACAT_non_missing(pvals) << "\t" << 
							n_samples << "\t" << 
							n_covar << "\t" << 
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



/*
void run_cis_eQTL_analysis_LMM_original(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

	std::cerr << "Starting eigendecomposition of GRM ... ";
	
	Eigen::SelfAdjointEigenSolver <Eigen::SparseMatrix<double>> GRM_eig(GRM);
	
	if (GRM_eig.info() != Eigen::Success){
		std::cerr << "FATAL ERROR: GRM decomposition failed!\n";
		abort();
	}
	std::cerr << "Done!\n";
	
	Eigen::VectorXd GRM_lambda = GRM_eig.eigenvalues();
	Eigen::SparseMatrix<double> L = GRM_eig.eigenvectors().sparseView(1.00, 5e-3);
	L.makeCompressed();

	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	std::cerr << "Started cis-eQTL analysis ...\n";
	
	if( rknorm_y ){
		std::cerr << "Rank-normalizing expression traits ... \n";
		rank_normalize(Y);
	}
	std::cerr << "Scaling expression traits ... \n";
	scale_and_center(Y);
	
	std::cerr << "Rotating expression and covariates ... ";
	Y = (L.transpose() * Y).eval();
	X = (L.transpose() * X).eval();
	std::cerr << "Done!\n";

	Eigen::MatrixXd U = get_half_hat_matrix(X);
	
	if( !global_opts::low_mem ){
		
		std::cerr << "Rotating genotypes ... ";
		g_data.genotypes = (L.transpose()*g_data.genotypes).eval();
		g_data.genotypes.makeCompressed();
		std::cerr << "Done!\n";
		
		std::cerr << "Calculating genotype-covariate covariance...\n";
		
		Eigen::MatrixXd UtG = U.transpose() * g_data.genotypes;
		std::cerr << "Calculating genotype residual variances ...";
		//Eigen::VectorXd SD_vec(UtG.cols());
		for( int i = 0; i < UtG.cols(); ++i)
		{
			
			g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();
			//SD_vec[i] = std::sqrt(g_data.var[i]);
		}
		std::cerr << "Done.\n";
	}
	
	
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
	
	std::string iter_cerr_suffix = " cis-eQTL blocks out of " + std::to_string(bm.size()) + " total";
	std::cerr << "Processed ";
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
				
				std::stringstream block_line;
				
				std::vector<double> pvals;
		
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
				
				std::cout << e_data.gene_id[jj] << "\t" << fit.sigma2 << "\t" << fit.phi << "\n";

				
				Eigen::MatrixXd XtDXi = (X.transpose() * fit.Vi * X).inverse();
				Eigen::MatrixXd XtDG = X.transpose() * fit.Vi * G_slice;
				Eigen::MatrixXd XtDXi_XtDG = XtDXi * XtDG;
				
				Eigen::VectorXd y_res = Y.col(jj) - X * XtDXi * X.transpose() * fit.Vi * Y.col(jj);
				
				double SSR_0 = y_res.dot(fit.Vi * Y.col(jj))/fit.sigma2;
				
				Eigen::VectorXd U_vec = G_slice.transpose() * fit.Vi * y_res/std::sqrt(fit.sigma2);
				
				Eigen::VectorXd diagV(U_vec.size());
				
				for (int si = 0; si < U_vec.size(); si++){
					diagV(si) = G_slice.col(si).dot(fit.Vi * G_slice.col(si)) - XtDG.col(si).dot(XtDXi_XtDG.col(si));	
				}
				
				e_data.stdev[jj] = std::sqrt(fit.sigma2);
				
				for( int ii = v_s, im = idx_s, si = 0; im < idx_s + idx_n; ii++, im++, si++){
				
					std::stringstream long_line;
				
					double g_yres_crossprod = U_vec(si);
					if( g_data.flipped[ii] ) g_yres_crossprod *= -1.0;
					
					if( !just_long ) block_line << "\t" << g_yres_crossprod;
					
					double beta, beta_se, zscore, pval_esnp;
					
					if( g_data.var[ii] > 0 ){
						beta = g_yres_crossprod/diagV(si);
						beta_se = std::sqrt( SSR_0/diagV(si) - beta*beta)/std::sqrt(n_samples - n_covar - 1);
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
			std::cerr << "\nERROR: " << bl << "; " << bm.bed_s[i] << ", " << n_e << "; " << bm.bcf_s[i] << ", " << n_g << "\n"; 
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

*/
