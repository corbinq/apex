# 1 "yax/src/fitModels.cpp.c"
#include "fitModels.hpp"
# 28 "yax/src/fitModels.cpp.c"
static const double min_hsq_thresh = 1e-5;

DiagonalXd calc_Vi(const double& phi, const Eigen::VectorXd& lambdas){
 Eigen::VectorXd vals = 1.00 + phi*lambdas.array();
 return vals.asDiagonal().inverse();
}

DiagonalXd calc_Psi(const double& phi, const Eigen::VectorXd& lambdas){

 double hsq = phi/(1.00 + phi);

 if( hsq <= min_hsq_thresh ){
  return Eigen::VectorXd::Zero(lambdas.size()).asDiagonal();
 }
 Eigen::VectorXd numer = hsq * lambdas;
 Eigen::VectorXd denom = hsq * lambdas;
 denom.array() += 1.00;
 return (numer.cwiseQuotient(denom)).asDiagonal();
}

Eigen::MatrixXd getVBeta(const std::vector<double>& hsq_v, const std::vector<double>& phi_v, const int& p){

 int q = hsq_v.size();

 Eigen::MatrixXd Beta(p, q);

 for(int i = 0; i < p; i++){
  for( int j = 0; j < q; j++){
   if( i == 0 ){
    Beta(i, j) = 1.00/(1.00 + phi_v[j]);
   }else{
    Beta(i, j) = std::pow(hsq_v[j], (double) i)/(1.00 + phi_v[j]);
   }
  }
 }
 return Beta;
}

Eigen::MatrixXd getPredParams( const std::vector<double>& vals ){
 int n = vals.size();
 Eigen::MatrixXd X(n,n);
 for(int i = 0; i < n; i++){
  for(int j = 0; j < n; j++){
   if ( i == 0 && j == 0 ){
    X(i,j) = 1.0;
   }else{
    X(i,j) = std::pow(vals[i], (double) j );
   }
  }
 }
 return (X.transpose() * X).ldlt().solve(X.transpose());
}

Eigen::VectorXd predV( const Eigen::MatrixXd& vv, const double& hsq ){
 Eigen::VectorXd out = vv.col(0);

 Eigen::VectorXd Beta(vv.cols());

 Beta(0) = 1.00;

 for(int i = 1; i < vv.cols(); i++){
  Beta(i) = std::pow(hsq, (double) i);
 }
 return vv * Beta;
}

void calcVBasis( Eigen::MatrixXd& V_mat, genotype_data& g_data, const Eigen::MatrixXd& C, const std::vector<double>& hsq_vals, const Eigen::MatrixXd& CtC, const Eigen::MatrixXd& CtC_i, const Eigen::MatrixXd& QtG, Eigen::MatrixXd& QtC, const Eigen::VectorXd Q_lambda ){

 int n_hsq = hsq_vals.size();
 int n_snps = g_data.genotypes.cols();

 std::vector<Eigen::MatrixXd> M_list;
 V_mat = Eigen::MatrixXd(n_snps, n_hsq);

 for( const double& hsq : hsq_vals ){
  M_list.push_back(
   (CtC - (QtC.transpose() * calc_Psi(hsq, Q_lambda) * QtC)).inverse()
  );
 }

 std::cerr << "Calculating genotype-covariate covariance ... \n";
 Eigen::MatrixXd CtG = (C.transpose() * g_data.genotypes).eval();

 std::vector<double> GtG_diag(n_snps);

 std::cerr << "Calculating genotype residual variances ...";

 for( int i = 0; i < n_snps; i++)
 {
  GtG_diag[i] = g_data.genotypes.col(i).squaredNorm();
  g_data.var[i] = GtG_diag[i] - CtG.col(i).dot(CtC_i * CtG.col(i));
  V_mat(i,0) = g_data.var[i];
 }
 std::cerr << "Done.\n";

 for(int j = 1; j < n_hsq; j++){
  std::cerr << "Calculating genotype variance basis (" << j+1 << "/" << n_hsq << ") ... \n";

  Eigen::MatrixXd Psi_a = calc_Psi(hsq_vals[j], Q_lambda) * QtG;
  Eigen::MatrixXd Res = CtG - QtC.transpose() * Psi_a;

  for( int i = 0; i < n_snps; i++)
  {
   V_mat(i,j) = (
    GtG_diag[i] - QtG.col(i).dot(Psi_a.col(i)) - Res.col(i).dot(M_list[j] * Res.col(i))
   );
  }
 }

 V_mat = (V_mat * (getPredParams(hsq_vals).transpose())).eval();

 std::cerr << "Done.\n";

 return;
}


void GRM_decomp( Eigen::SparseMatrix<double>& GRM, const std::vector<int>& relateds, PermutXd& Tr, Eigen::SparseMatrix<double>& L, Eigen::VectorXd& L_lambda, Eigen::SparseMatrix<double>& Q, Eigen::VectorXd& Q_lambda ){

 double delta_thresh = 0.01;
 double drop0_thresh = 1e-3;

 std::vector<int> unrelateds;
 Tr = PermutXd(GRM.rows());

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

 L_lambda = GRM_eig.eigenvalues();

 L_lambda.conservativeResize(nrel + nind);
 for(int i = nrel; i < nrel + nind; i++){
  L_lambda(i) = 1;
 }


 L = GRM_eig.eigenvectors().sparseView(drop0_thresh, 1.0 - std::numeric_limits<double>::epsilon());

 L.conservativeResize(nrel + nind,nrel + nind);

 for(int i = nrel; i < nrel + nind; i++){
  L.coeffRef(i,i) = 1.00;
 }

 L.makeCompressed();

 std::vector<int> kp;
 for( int i = 0; i < L_lambda.size(); i++ ){
  if( std::abs( L_lambda(i) - 1 ) > delta_thresh ){
   kp.push_back(i);
  }
 }


 using td = Eigen::Triplet<double>;
 Eigen::SparseMatrix<double> PC_sel(L_lambda.size(),kp.size());
 std::vector<td> PC_trips;


 Q_lambda = Eigen::VectorXd(kp.size());

 i = 0;
 for( const int& k : kp){
  PC_trips.push_back(td(k,i,1.00));
  Q_lambda(i) = L_lambda(k) - 1.00;
  i++;
 }
 PC_sel.setFromTriplets(PC_trips.begin(), PC_trips.end());

 std::cerr << "Selected " << kp.size() << " eigenvectors.\n";


 Q = (L * PC_sel).eval();
 Q.makeCompressed();

 return;
}


void meta_svar_sumstat::condition_on_het(const int& k){
 kept_snps.push_back(k);
 for( int s = 0; s < ss.size(); s++ ){


  ss[s].condition_on( std::vector<int>(1, k), vg.Covar_perStudy(s, k) );

 }
 update_meta_ss();
}

lm_output lm_from_sumstats( const Eigen::VectorXd& U, const Eigen::VectorXd& V, const double& n, const double& df_0, const double& stdev, const Eigen::VectorXd& U_0, const Eigen::MatrixXd& J_0, const Eigen::MatrixXd& Cov, const bool& check_filter, const std::vector<bool>& exclude ){

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

    double se = stdev * std::sqrt(SSE_i - SCORE*SCORE/VARSC) / std::sqrt(df*VARSC);
    double pval = -99;

    double PVAR = SCORE*SCORE/VARSC;
    double STAT = df*PVAR/(SSE_i - PVAR);

    if( SSE_i - PVAR > 0 && STAT > 0 && PVAR > 0 ){
     pval = pf( STAT, 1.0, df, true );
    }
    out.push_back(beta, se, pval);
   }else{

    out.push_back(-99, -99, -99);
   }
  }
 }

 return out;
}

int which_min( const std::vector<double>& p, bool gt0 ){
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
 std::vector<int> kept_snps = seq_int(keep.size());
 for(const int k : keep ){
  std::vector<int> kept_snps_not_k = kept_snps;
  kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );

  double current_pvalue;

  Eigen::VectorXd U_k = U(std::vector<int>(1,k));
  Eigen::VectorXd V_k = V(std::vector<int>(1,k));

  Eigen::VectorXd U_0k = U_0(kept_snps_not_k);
  Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k);
  Eigen::MatrixXd Cov_k = J_0(std::vector<int>(1,k_i), kept_snps_not_k);

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
  std::cerr << "U.size() != V.size()" << "\n";
  exit(1);
 }

 if( U.size() <= 0 ){
  std::cerr << "U.size() <= 0" << "\n";
  exit(1);
 }

 int n_var = U.size();
 int nk = 0;

 std::vector<bool> excl(n_var,false);

 lm_output reg0;

 double alpha_thresh = global_opts::LM_ALPHA;

 while( 1 )
 {

  int steps_taken = 0;

  lm_output reg = lm_from_sumstats(U, V, n, m, stdev, U_0, J_0, Cov, true, excl);


  if( nk == 0 && reg0.beta.size() == 0 ){
   reg0 = reg;
  }

  int wk = which_min(reg.pval, true);
  double adj_pval = ACAT_non_missing(reg.pval);



  double pval_check = -99;

  if( wk >= 0 ){
   if( global_opts::step_marginal ){
    pval_check = reg.pval[wk];
   }else{
    pval_check = adj_pval;
   }
  }




  if( (pval_check >= 0 && pval_check < alpha_thresh) || (keep.size() == 0 && wk >= 0 ) )
  {
   keep.push_back(wk);


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



   J_0 = vget.Var(keep);

   U_0(nk-1) = U(wk);

   if( pval_check <= alpha_thresh || keep.size() < 2 ){


    Eigen::VectorXd new_cov = vget.Covar(wk);


    if( new_cov.size() != Cov.rows() )
    {
     std::cerr << "new_cov.size() != Cov.rows()" << "\n";
     exit(1);
    }

    std::vector<int> b_list;
    std::vector<double> r_list;

    for(int i = 0; i < n_var; ++i)
    {

     Cov(i, nk-1) = new_cov(i);

     if( global_opts::RSQ_BUDDY < 1.00 ){
      double corr = new_cov(i)/std::sqrt(V(i)*V(wk));



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

  }




  if( global_opts::backward_step && nk > 1 ){

   double max_joint_pvalue = 0;
   int k_rm = -1;

   check_joint_pvalues(k_rm, max_joint_pvalue, U, V, U_0, J_0, n, m);


   if( (max_joint_pvalue > alpha_thresh || max_joint_pvalue < 0 ) && k_rm < nk - 1 && k_rm >= 0 ){

    std::vector<int> kept_snps_not_k_rm = seq_int(nk);
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



    keep.erase( keep.begin() + k_rm );

    nk--;


    steps_taken++;
   }
  }


  if( steps_taken == 0 ){
   break;
  }
 }

 int k_i = 0;
 std::vector<int> kept_snps = seq_int(nk);

 for( const int k : keep ){

  std::vector<int> kept_snps_not_k = kept_snps;
  kept_snps_not_k.erase( kept_snps_not_k.begin() + k_i );

  if( kept_snps_not_k.size() == 0 ){

   beta.push_back(reg0.beta[k]);
   se.push_back(reg0.se[k]);
   pval_joint.push_back(reg0.pval[k]);

  }else{

   Eigen::VectorXd U_k = U(std::vector<int>(1,k));
   Eigen::VectorXd V_k = V(std::vector<int>(1,k));

   Eigen::VectorXd U_0k = U_0(kept_snps_not_k);
   Eigen::MatrixXd J_0k = J_0(kept_snps_not_k, kept_snps_not_k);
   Eigen::MatrixXd Cov_k = J_0(std::vector<int>(1,k_i), kept_snps_not_k);

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
  std::cout <<
   beta[i] << "\t" <<
   se[i] << "\t" <<
   pval[i] << "\n";
 }
}
