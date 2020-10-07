# 1 "yax/src/cisMapping.cpp.c"
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

  for( int i = 0; i < UtG.cols(); ++i)
  {

   g_data.var[i] = g_data.genotypes.col(i).squaredNorm() - UtG.col(i).squaredNorm();

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

    Eigen::VectorXd UtG_block_sqnm = (U.transpose() * G).colwise().squaredNorm().eval();

    for( int si = bm.bcf_s[i], ii = 0; si < bm.bcf_s[i] + n_g; ++si, ++ii)
    {
     g_data.var[si] = G.col(ii).squaredNorm() - UtG_block_sqnm(ii);

    }

    out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();

   }else{

    const Eigen::SparseMatrix<double>& G = g_data.genotypes.middleCols(bm.bcf_s[i], n_g);

    out_mat = (Y_res.middleRows(bm.bed_s[i], n_e)* G).transpose().eval();



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
        e_data.stdev[jj]*beta_se << "\t" <<
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

 }

 if( !just_long ){
  bgzf_close(block_file);
  bgzf_close(bed_block_file);
  build_tabix_index(block_file_path, 1);
  build_tabix_index(bed_block_file_path, 1);
 }

 return;
}

void run_cis_eQTL_analysis_LMM(bcf_srs_t*& sr, bcf_hdr_t*& hdr,genotype_data& g_data, table& c_data, bed_data& e_data, Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, block_intervals& bm, const bool& rknorm_y, const bool& rknorm_r, const bool& make_sumstat, const bool& make_long, const bool& just_long)
{

 PermutXd Tr;
 Eigen::SparseMatrix<double> Q;
 Eigen::VectorXd Q_lambda;
 Eigen::SparseMatrix<double> L;
 Eigen::VectorXd GRM_lambda;

 GRM_decomp(GRM, relateds, Tr, L, GRM_lambda, Q, Q_lambda);

 std::cerr << "Reordering trait and covariate matrices ...\n";

 Eigen::MatrixXd Y = (Tr * e_data.data_matrix).eval();
 Eigen::MatrixXd C = (Tr * c_data.data_matrix).eval();

 std::cerr << "Reordering genotypes ...\n";

 g_data.genotypes = (Tr * g_data.genotypes).eval();


 double n_traits = Y.cols();
 double n_samples = Y.rows();
 double n_snps = g_data.n_variants;
 double n_covar = C.cols();



 if( rknorm_y ){
  std::cerr << "Rank-normalizing expression traits ... \n";
  rank_normalize(Y);
 }
 std::cerr << "Scaling expression traits ... \n";
 scale_and_center(Y);




 std::cerr << "Calculating partial rotations ...\n";
 Eigen::MatrixXd QtC = (Q.transpose() * C).eval();
 Eigen::MatrixXd QtG = (Q.transpose() * g_data.genotypes).eval();
 Eigen::MatrixXd QtY = (Q.transpose() * Y).eval();
 Eigen::MatrixXd CtY = (C.transpose() * Y).eval();
 Eigen::MatrixXd CtC = (C.transpose() * C).eval();
 Eigen::MatrixXd CtC_i = (CtC.inverse()).eval();

 std::cerr << "Rotating expression and covariates ... ";
 Eigen::MatrixXd Y_raw = (Y).eval();
 Y = (L.transpose() * Y).eval();
 Eigen::MatrixXd X = (L.transpose() * C).eval();
 std::cerr << "Done.\n";

 std::vector<double> hsq_vals{0.0, 0.5, 1.0};

 Eigen::MatrixXd V_mat;
 calcVBasis(V_mat, g_data, C, hsq_vals, CtC, CtC_i, QtG, QtC, Q_lambda);

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
# 382 "yax/src/cisMapping.cpp.c"
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
    int v_e = v_s + n_slice - 1;

    int pos_s = g_data.pos[v_s];
    int pos_e = g_data.pos[v_e];

    const Eigen::SparseMatrix<double>& G_slice = G.middleCols(idx_s, idx_n);
    const Eigen::MatrixXd& QtG_slice = QtG.middleCols(s_slice, n_slice).middleCols(idx_s, idx_n);


    LMM_fitter fit(X, Y.col(jj), GRM_lambda);
    fit.fit_REML();

    const double& sigma2 = fit.sigma2;
    const double& phi = fit.phi;
    double tau2 = phi*sigma2;
    double hsq = tau2 / (tau2 + sigma2);
    double S2 = tau2 + sigma2;





    std::stringstream theta_line;

    theta_line <<
     clean_chrom(e_data.chr[jj]) << "\t" <<
     e_data.start[jj] << "\t" <<
     e_data.end[jj] << "\t" <<
     e_data.gene_id[jj] << "\t" <<
     sigma2 << "\t" <<
     tau2 << "\t" <<
     phi << "\n";

    write_to_bgzf(theta_line.str().c_str(), theta_file);
# 477 "yax/src/cisMapping.cpp.c"
    double ADJ = 1.00/(1.00 + phi);

    DiagonalXd Psi = calc_Psi(phi, Q_lambda);
    Eigen::MatrixXd XtDXi = ((CtC - QtC.transpose() * Psi * QtC ).inverse())/ADJ;
# 518 "yax/src/cisMapping.cpp.c"
    Eigen::VectorXd y_res = Y.col(jj) - X * XtDXi * X.transpose() * fit.Vi * Y.col(jj);
    double SSR_0 = y_res.dot(fit.Vi * Y.col(jj))/sigma2;

    Eigen::VectorXd Py = (L * (fit.Vi * y_res/std::sqrt(sigma2)));

    Eigen::VectorXd U_vec = G_slice.transpose() * Py;
# 541 "yax/src/cisMapping.cpp.c"
    if( v_s + U_vec.size() - 1 >= V_mat.rows() ){
     std::cerr << "\n";
     std::cerr << " v_s = " << v_s << "\n";
     std::cerr << " U_vec.size() = " << U_vec.size() << "\n";
     std::cerr << " V_mat.rows() = " << V_mat.rows() << "\n";
     std::cerr << "\n";
     continue;
    }
    Eigen::VectorXd diagV = predV(V_mat.middleRows(v_s, U_vec.size()), hsq)/(1.00 + fit.phi);

    e_data.stdev[jj] = std::sqrt(fit.sigma2);

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

     pvals.push_back(pval_esnp);

     if( write_long ){

      long_line <<
       clean_chrom(g_data.chr[ii]) << "\t" <<
       g_data.pos[ii] << "\t" <<
       g_data.ref[ii] << "\t" <<
       g_data.alt[ii] << "\t" <<
       e_data.gene_id[jj] << "\t" <<
       e_data.stdev[jj]*beta << "\t" <<
       e_data.stdev[jj]*beta_se << "\t" <<
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

 bgzf_close(theta_file);

 if ( write_long ){
  bgzf_close(long_file);

 }

 if( !just_long ){
  bgzf_close(block_file);
  bgzf_close(bed_block_file);
  build_tabix_index(block_file_path, 1);
  build_tabix_index(bed_block_file_path, 1);
 }

 return;
}
