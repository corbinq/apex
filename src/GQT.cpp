# 1 "yax/src/GQT.cpp.c"
#include "GQT.hpp"
#include <string>
# 12 "yax/src/GQT.cpp.c"
std::string prefix = "";
bool use_low_mem = false;
double rsq_buddy = 2.0;
double rsq_prune = 0.80;
double pval_thresh = 5e-5;
int window_size = 1000000;
std::vector<std::string> target_genes;
bool use_ivw_1 = false;
bool use_ds = false;
bool trim_gene_ids = false;
bool stepwise_backward_step = false;
bool t_hom = false;
bool t_het = false;
bool t_acat = false;
bool stepwise_marginal_thresh = false;





int cis(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int trans(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int meta(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int store(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);

using mode_fun = std::function<int(const std::string &, std::vector<std::string>::const_iterator, std::vector<std::string>::const_iterator)>;

std::string help_string =
"\n"
"  GQT: Toolkit for xQTL analysis\n"
"     (c) 2019-2020 Corbin Quick and Li Guan.\n"
"\n"
"  Usage and options:\n"
"     ./gqt [mode] --help       Print help menu for [mode].\n"
"\n"
"  Analysis modes:\n"
"     ./gqt cis {OPTIONS}       cis-xQTL analysis from\n"
"                                 individual-level data.\n"
"\n"
"     ./gqt trans {OPTIONS}     trans-xQTL analysis from\n"
"                                 individual-level data.\n"
"\n"
"     ./gqt meta {OPTIONS}      Single and multi-variant\n"
"                                 xQTL meta-analysis from\n"
"                                 sumstat and vcov files.\n"
"\n"
"     ./gqt store {OPTIONS}     Store vcov (LD) data for\n"
"                                 xQTL meta-analysis or\n"
"                                 data sharing.\n"
"\n"
"  Contact: corbinq@gmail.com\n";





int parseModeArgs(args::ArgumentParser& p, std::vector<std::string>::const_iterator& args_begin, std::vector<std::string>::const_iterator& args_end){

    try
    {
        p.ParseArgs(args_begin, args_end);
    }
 catch (const args::Completion& e)
    {
        std::cout << e.what();
        return 0;
    }
    catch (args::Help)
    {
  restore_cursor();
        std::cout << p;
        exit(0);
    }
    catch (args::ParseError e)
    {
  restore_cursor();
        std::cerr << e.what() << "\n";
        std::cerr << p;
        exit(1);
    }
    catch (args::ValidationError e)
    {
  restore_cursor();
        std::cerr << e.what() << "\n";
        std::cerr << p;
        exit(1);
    }
}





int main(int argc, char* argv[])
{





 auto lam_kill =
      [] (int i) { restore_cursor(); std::cerr << "\nKilled.\n" << "\n"; exit(0); };

    signal(SIGINT, lam_kill);
    signal(SIGABRT, lam_kill);
    signal(SIGTERM, lam_kill);


 at_quick_exit (restore_cursor);
 atexit (restore_cursor);






    std::unordered_map<std::string, mode_fun> map{
        {"cis", cis},
        {"trans", trans},
        {"meta", meta},
        {"store", store}
 };


 args::ArgumentParser p0("gqt: GWAS/QTL Toolkit.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help0(p0, "help", "Display this help menu", {'h', "help"});

 p0.Prog(argv[0]);

 std::string mode_summary = "\ncis: cis-xQTL analysis.\n\ntrans: trans-xQTL analysis.\nmeta: xQTL meta-analysis.\n\nstore: store xQTL vcov (LD) data.\n\n";

 args::MapPositional<std::string, mode_fun> mode(p0, "mode", mode_summary, map);
 mode.KickOut(true);

 const std::vector<std::string> args(argv + 1, argv + argc);

 try
    {
        auto next = p0.ParseArgs(args);
        if (mode)
        {
   return args::get(mode)(argv[0], next, std::end(args));
        } else
        {
            std::cout << help_string;
        }
    }
    catch (args::Help)
    {
        std::cout << help_string;
        return 0;
    }
    catch (args::Error e)
    {
        std::cerr << "\nUnknown command line argument(s).\nPrinting help menu:\n" << std::endl;
        std::cerr << help_string;
        return 1;
    }
    return 0;
}


int cis(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){

 args::ArgumentParser p("gqt cis: cis-xQTL analysis.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
 args::CompletionFlag completion(p, {"complete"});

 p.Prog(progname);

 args::Group analysis_args(p, "Analysis options");
  args::Flag fit_null(analysis_args, "", "Estimate and store LMM null model parameters.", {"fit-null"});
  args::ValueFlag<std::string> theta_arg(analysis_args, "", "Use stored LMM null model parameters.", {"theta-file"});

 args::Group cis_args(p, "Output options");
  args::Flag make_long(cis_args, "", "Write cis-eQTL results in long table format.", {'l', "long"});
  args::Flag just_long(cis_args, "", "Only write long-table cis-eQTL results.", {'L', "just-long"});

 args::Group scale_args(p, "Scale and transform options");
  args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
  args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
  args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
  args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});

 args::Group input_args(p, "Input files");
  args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
  args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});

  args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for eQTL analysis.", {'b', "bed", "expression"});
  args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
  args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
  args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
# 215 "yax/src/GQT.cpp.c"
 args::Group opt_args(p, "General options");
  args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
  args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
  args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
  args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
  args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-eQTL or gene-based analysis.", {'w', "window"});
  args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});

  args::ValueFlag<double> pval_arg(opt_args, "", "P-value threshold for stepwise procedures.", {"pvalue"});
  args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});





 parseModeArgs(p, beginargs, endargs);





 std::string theta_path = args::get(theta_arg);

 prefix = args::get(out_arg);
 std::string e_path = args::get(bed_arg);
 std::string g_path = args::get(bcf_arg);
 std::string c_path = args::get(cov_arg);






 std::string grm_path = args::get(grm_arg);
 std::string kin_path = args::get(kin_arg);
 double grm_scale = 1.00;

 if( kin_path != "" ){
  if( grm_path != "" ){
   std::cerr << "ERROR: Specify --kin or --grm, but not both.\n";
   abort();
  }
  grm_path = kin_path;
  grm_scale = 2.00;
 }






 trim_gene_ids = (bool) trim_ids;

 std::string region = args::get(region_arg);
 std::string gtds = args::get(gtds_arg);
 target_genes = split_string(args::get(gene_arg), ',');


 if( gtds == "" || gtds == "GT" ){
  use_ds = false;
 }else if( gtds == "DS" ){
  use_ds = true;
 }else{
  std::cerr << "Invalid specification --field " << gtds << "\n";
  std::cerr << "Valid options are \"GT\" or \"DS\". Exiting. \n";
  return 1;
 }

 std::string window_size_s = args::get(window_arg);

 if( window_size_s == "" ){
  window_size = 1000000;
 }else{
  window_size = stoi(window_size_s);
  std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
 }

 if( window_size <= 0 ) window_size = 1000000;






 use_low_mem = (bool) low_mem;

 int nthreads = args::get(threads_arg);
 if( nthreads >= 1 ){
  omp_set_num_threads(nthreads);
  Eigen::setNbThreads(nthreads);
 }
 std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";

 global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_step, t_hom, t_het, t_acat, stepwise_marginal_thresh);

 if( prefix == "" )
 {
  restore_cursor();
 }else
 {
  hide_cursor();
 }

 if( prefix == "" ){
  std::cerr << "Error: Output prefix not specified. Try --help to see options.\n";
  return 0;
 }

 int max_var = 100000000;

 int n_var = 0;

 std::vector<int> variants_per_chrom;

 genotype_data g_data;
 table c_data;
 bed_data e_data;

 std::vector<std::string> bcf_chroms = get_chroms(g_path, variants_per_chrom);
 std::vector<std::string> bed_chroms = get_chroms(e_path);
 std::vector<std::string> keep_chroms = intersect_ids(bcf_chroms,bed_chroms);

 for(int i = 0; i < bcf_chroms.size(); i++ ){
  if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
   n_var += variants_per_chrom[i];
  }
 }


 for(const auto& c : keep_chroms){
  std::cerr << c << ",";
 }
 std::cerr << "\b present in both bcf and bed file.\n";
 std::cerr << n_var << " total variants on selected chromosomes.\n\n";

 bcf_srs_t *sr = bcf_sr_init();

 if( region != "" ){

  std::cerr << "Setting region to " << region << " in bcf file ... \n";
  bcf_sr_set_regions(sr, region.c_str(), 0);

 }else{
  std::string region_string = "";
  for(std::string& chr : keep_chroms){
   region_string += ( region_string=="" ? "" : "," ) + chr;
  }
  bcf_sr_set_regions(sr, region_string.c_str(), 0);
 }

 bcf_sr_add_reader(sr, g_path.c_str());
 bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);


 g_data.read_bcf_header(hdr);
 std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";


 if( c_path == "" ){
  std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
  std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
  std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
  std::cerr << "    --rankNormal-resid for trait residuals.\n";
 }else{
  c_data.readHeader(c_path.c_str());
  std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
 }


 e_data.readBedHeader(e_path.c_str());
 std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";

 std::vector<std::string> intersected_samples;
 if( c_path == "" ){
  intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
 }else{
  intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
 }



 std::vector<std::string> intersected_samples_gto = g_data.ids.file;
 for(int i = 0; i < intersected_samples_gto.size(); ){
  if( has_element(intersected_samples, intersected_samples_gto[i]) ){
   i++;
  }else{
   intersected_samples_gto.erase(intersected_samples_gto.begin() + i);
  }
 }

 intersected_samples = intersected_samples_gto;

 std::cerr << "Found " << intersected_samples.size() << " samples in common across all three files.\n\n";


 g_data.ids.setKeepIDs(intersected_samples);
 g_data.n_samples = intersected_samples.size();

 e_data.ids.setKeepIDs(intersected_samples);
 if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);

 std::vector<std::string> keep_regions = keep_chroms;
 if( region != "" ){
  keep_regions.clear();
  keep_regions.push_back(region);
 }



 if( c_path == "" ){
   c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
 }else{
  c_data.readFile(c_path.c_str());
  std::cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";

  if( !no_scale_x ){
   scale_and_center(c_data.data_matrix);
  }
  appendInterceptColumn(c_data.data_matrix);
 }

 e_data.readBedFile(e_path.c_str(),keep_regions);

 std::cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";



 g_data.read_bcf_variants(sr, hdr, n_var, !low_mem, !low_mem);

 if( g_data.chr.size() == 0 ){
  std::cerr << "\nNo variants present in specified region(s). Exiting.\n\n";
  return 0;
 }

 if( low_mem || fit_null ){
  clear_line_cerr();
  std::cerr << "Processed variant data for " << n_var << " variants.\n\n";
  g_data.genotypes.resize(0,0);
 }else{
  std::cerr << "\rFreezing genotype data for " << n_var << " variants ... \r";
  g_data.freeze_genotypes();
  clear_line_cerr();
  std::cerr << "Processed genotype data for " << n_var << " variants.\n\n";
 }

 block_intervals bm;
 bm.make_blocks(e_data, g_data, window_size);

 Eigen::MatrixXd &Y = e_data.data_matrix;
 Eigen::MatrixXd &X = c_data.data_matrix;

 Eigen::SparseMatrix<double> GRM;
 std::vector<int> relateds;
 if( grm_path != "" ){
  read_sparse_GRM(grm_path, GRM, intersected_samples, grm_scale, 3, relateds);
 }

 if( fit_null ){
  if( grm_path == "" ){
   std::cerr << "Error: GRM is required to fit null models.\n";
   return 1;
  }
  fit_LMM_null_models(c_data, e_data, GRM, relateds, rknorm_y, rknorm_r);
  std::cerr << "\nNull model estimation complete. \nSpecify --theta-file {prefix}.theta.gz to re-use estimates for analysis.\n";
  return 0;
 }

 if( grm_path == "" ){
  run_cis_eQTL_analysis(sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r, true, make_long, just_long);
 }else{
  run_cis_eQTL_analysis_LMM(sr, hdr, g_data, c_data, e_data, GRM, relateds, bm, rknorm_y, rknorm_r, true, make_long, just_long);
 }

    return 0;
};


int trans(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){





 args::ArgumentParser p("gqt trans: trans-xQTL analysis.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
 args::CompletionFlag completion(p, {"complete"});

 p.Prog(progname);

 args::Group analysis_args(p, "Analysis options");
  args::Flag fit_null(analysis_args, "", "Estimate and store LMM null model parameters.", {"fit-null"});
  args::ValueFlag<std::string> theta_arg(analysis_args, "", "Use stored LMM null model parameters.", {"theta-file"});

 args::Group scale_args(p, "Scale and transform options");
  args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
  args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
  args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
  args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});

 args::Group input_args(p, "Input files");
  args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
  args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});

  args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for eQTL analysis.", {'b', "bed", "expression"});
  args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
  args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
  args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
# 533 "yax/src/GQT.cpp.c"
 args::Group opt_args(p, "General options");
  args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
  args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
  args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
  args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
  args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-eQTL or gene-based analysis.", {'w', "window"});
  args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});
  args::ValueFlag<int> blocks_arg(opt_args, "100", "Number of variants per trans-xQTL block.", {"block-size"});
  args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});





 parseModeArgs(p, beginargs, endargs);





 std::string theta_path = args::get(theta_arg);

 prefix = args::get(out_arg);
 std::string e_path = args::get(bed_arg);
 std::string g_path = args::get(bcf_arg);
 std::string c_path = args::get(cov_arg);






 std::string grm_path = args::get(grm_arg);
 std::string kin_path = args::get(kin_arg);
 double grm_scale = 1.00;

 if( kin_path != "" ){
  if( grm_path != "" ){
   std::cerr << "ERROR: Specify --kin or --grm, but not both.\n";
   abort();
  }
  grm_path = kin_path;
  grm_scale = 2.00;
 }






 trim_gene_ids = (bool) trim_ids;

 std::string region = args::get(region_arg);
 std::string gtds = args::get(gtds_arg);
 target_genes = split_string(args::get(gene_arg), ',');

 int block_size = args::get(blocks_arg);

 if( block_size <= 0 ){
  block_size = 100;
 }


 use_ds;
 if( gtds == "" || gtds == "GT" ){
  use_ds = false;
 }else if( gtds == "DS" ){
  use_ds = true;
 }else{
  std::cerr << "Invalid specification --field " << gtds << "\n";
  std::cerr << "Valid options are \"GT\" or \"DS\". Exiting. \n";
  return 1;
 }



 std::string window_size_s = args::get(window_arg);

 if( window_size_s == "" ){
  window_size = 1000000;
 }else{
  window_size = stoi(window_size_s);
  std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
 }





 use_low_mem = (bool) low_mem;

 int nthreads = args::get(threads_arg);
 if( nthreads >= 1 ){
  omp_set_num_threads(nthreads);
  Eigen::setNbThreads(nthreads);
 }
 std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";

 global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_step, t_hom, t_het, t_acat, stepwise_marginal_thresh);

 if( prefix == "" )
 {
  restore_cursor();
 }else
 {
  hide_cursor();
 }

 if( prefix == "" ){
  std::cerr << "Error: Output prefix not specified. Try --help to see options.\n";
  return 0;
 }

 int max_var = 100000000;

 int n_var = 0;

 std::vector<int> variants_per_chrom;

 genotype_data g_data;
 table c_data;
 bed_data e_data;

 std::vector<std::string> bcf_chroms = get_chroms(g_path, variants_per_chrom);
 std::vector<std::string> bed_chroms = get_chroms(e_path);
 std::vector<std::string> keep_chroms = intersect_ids(bcf_chroms,bed_chroms);

 for(int i = 0; i < bcf_chroms.size(); i++ ){
  if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
   n_var += variants_per_chrom[i];
  }
 }


 for(const auto& c : keep_chroms){
  std::cerr << c << ",";
 }
 std::cerr << "\b present in both bcf and bed file.\n";
 std::cerr << n_var << " total variants on selected chromosomes.\n\n";

 bcf_srs_t *sr = bcf_sr_init();

 if( region != "" ){

  std::cerr << "Setting region to " << region << " in bcf file ... \n";
  bcf_sr_set_regions(sr, region.c_str(), 0);

 }else{
  std::string region_string = "";
  for(std::string& chr : keep_chroms){
   region_string += ( region_string=="" ? "" : "," ) + chr;
  }
  bcf_sr_set_regions(sr, region_string.c_str(), 0);
 }

 bcf_sr_add_reader(sr, g_path.c_str());
 bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);


 g_data.read_bcf_header(hdr);
 std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";


 if( c_path == "" ){
  std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
  std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
  std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
  std::cerr << "    --rankNormal-resid for trait residuals.\n";
 }else{
  c_data.readHeader(c_path.c_str());
  std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
 }


 e_data.readBedHeader(e_path.c_str());
 std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";

 std::vector<std::string> intersected_samples;
 if( c_path == "" ){
  intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
 }else{
  intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
 }



 std::vector<std::string> intersected_samples_gto = g_data.ids.file;
 for(int i = 0; i < intersected_samples_gto.size(); ){
  if( has_element(intersected_samples, intersected_samples_gto[i]) ){
   i++;
  }else{
   intersected_samples_gto.erase(intersected_samples_gto.begin() + i);
  }
 }

 intersected_samples = intersected_samples_gto;

 std::cerr << "Found " << intersected_samples.size() << " samples in common across all three files.\n\n";


 g_data.ids.setKeepIDs(intersected_samples);
 g_data.n_samples = intersected_samples.size();

 e_data.ids.setKeepIDs(intersected_samples);
 if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);

 std::vector<std::string> keep_regions = keep_chroms;
 if( region != "" ){
  keep_regions.clear();
  keep_regions.push_back(region);
 }



 if( c_path == "" ){
   c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
 }else{
  c_data.readFile(c_path.c_str());
  std::cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";

  if( !no_scale_x ){
   scale_and_center(c_data.data_matrix);
  }
  appendInterceptColumn(c_data.data_matrix);
 }

 e_data.readBedFile(e_path.c_str(), bed_chroms);
 std::cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";

 g_data.read_bcf_variants(sr, hdr, n_var, !low_mem, !low_mem);

 if( g_data.chr.size() == 0 ){
  std::cerr << "\nNo variants present in specified region(s). Exiting.\n\n";
  return 0;
 }

 if( low_mem || fit_null ){
  clear_line_cerr();
  std::cerr << "Processed variant data for " << n_var << " variants.\n\n";
  g_data.genotypes.resize(0,0);
 }else{
  std::cerr << "\rFreezing genotype data for " << n_var << " variants ... \r";
  g_data.freeze_genotypes();
  clear_line_cerr();
  std::cerr << "Processed genotype data for " << n_var << " variants.\n\n";
 }

 block_intervals bm;
 Eigen::MatrixXd &Y = e_data.data_matrix;
 Eigen::MatrixXd &X = c_data.data_matrix;

 Eigen::SparseMatrix<double> GRM;
 std::vector<int> relateds;
 if( grm_path != "" ){
  read_sparse_GRM(grm_path, GRM, intersected_samples, grm_scale, 3, relateds);
 }

 bool make_sumstat = true;
 bool make_long = true;

 if( fit_null ){
  if( grm_path == "" ){
   std::cerr << "Error: GRM is required to fit null models.\n";
   return 1;
  }
  fit_LMM_null_models(c_data, e_data, GRM, relateds, rknorm_y, rknorm_r);
  std::cerr << "Analysis complete. Specify --null-params {theta-file} to re-use null model estimates.\n";
  return 0;
 }

 if( grm_path == "" ){
  run_trans_eQTL_analysis(sr, hdr, g_data, c_data, e_data, rknorm_y, rknorm_r, make_sumstat, make_long, block_size);
 }else{
  run_trans_eQTL_analysis_LMM(sr, hdr, g_data, c_data, e_data, GRM, relateds, rknorm_y, rknorm_r, make_sumstat, make_long, block_size, theta_path);
 }

    return 0;
};


int meta(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){





 args::ArgumentParser p("gqt meta: Meta-analysis of xQTL studies.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
 args::CompletionFlag completion(p, {"complete"});

 p.Prog(progname);

 args::Group inp_args(p, "Input and output file paths");
  args::ValueFlag<std::string> meta_arg(inp_args, "", "Sumstat file prefixes for meta-analysis (comma-separated).", {"sumstats"});
  args::ValueFlag<std::string> out_arg(inp_args, "", "Prefix for output files.", {'o', "prefix", "out"});

 args::Group meta_args(p, "Analysis mode");
  args::Flag meta_svar(meta_args, "", "Perform single-variant meta-analysis.", {"meta"});
  args::Flag meta_stepwise(meta_args, "", "Perform stepwise (conditional) meta-analysis.", {"stepwise"});

 args::Group analysis_args(p, "Analysis options");
  args::ValueFlag<std::string> test_arg(analysis_args, "het,alt", "Meta-analysis test forms (comma-separated combination of 'het,hom,alt').", {'t', "tests"});
  args::Flag ivw_1(analysis_args, "", "Calculate IVW weights under the alternative hypothesis.", {"ivw1"});
  args::Flag het_meta(analysis_args, "", "Allow heterogeneous genotype effects across studies in stepwise meta-analysis.", {"het"});
  args::ValueFlag<double> rsq_arg(analysis_args, "", "Maximum rsq threshold.", {"rsq"});
  args::ValueFlag<double> buddy_arg(analysis_args, "", "Print LD buddies (specify rsq threshold).", {"buddies"});
  args::Flag use_marginal_pval(analysis_args, "", "Apply threshold to marginal rather than ACAT p-value in stepwise algorithm (no adjustment for no. variants).", {"marginal"});
  args::ValueFlag<double> pval_arg(analysis_args, "", "P-value threshold for stepwise procedures.", {"pvalue"});
  args::Flag backward_step(analysis_args, "", "Perform backward step in stepwise selection.", {"backward"});






 args::Group opt_args(p, "General options");
  args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
  args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
  args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-eQTL or gene-based analysis.", {'w', "window"});
  args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});

  args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});

 args::Group vcov_args(p, "Querying vcov files");
  args::Flag query_vcov(vcov_args, "", "Query meta-analysis vcov data (file list from --sumstats). Can be used with --region.", {"query"});
  args::ValueFlag<int> query_pos_arg(vcov_args, "pos", "Only get vcov data for specified variant (against all others in region).", {"pos"});
  args::Flag query_u_vcov(vcov_args, "", "Just return uncentered vcov data from bin file (default: centered v-cov matrix).", {"just-bin"});







 parseModeArgs(p, beginargs, endargs);





 use_ivw_1 = (bool) ivw_1;
 stepwise_backward_step = (bool) backward_step;
 stepwise_marginal_thresh = (bool) use_marginal_pval;

 std::string tests = args::get(test_arg);

 if( tests == "" ){
  tests = "hom,acat";
 }
 for( char& c_i : tests){
  c_i = tolower(c_i);
 }
 for( const std::string& test_i : split_string(tests, ',') ){
  if( test_i == "hom" ){
   t_hom = true;
  }else if( test_i == "het" ){
   t_het = true;
  }else if( test_i == "acat" ){
   t_acat = true;
  }else{
   std::cerr << "Error: Unknown test '" << test_i << "'.\n";
   std::cerr << "Valid options are het,hom,acat (comma-separated combinations)\n";
   abort();
  }
 }
 if( !t_hom && !t_het && !t_acat ){
  std::cerr << "Error: No valid test options specified.\n";
  std::cerr << "Valid options are het,hom,acat (comma-separated combinations)\n";
  abort();
 }





 prefix = args::get(out_arg);





 trim_gene_ids = (bool) trim_ids;

 std::string region = args::get(region_arg);
 target_genes = split_string(args::get(gene_arg), ',');

 std::string window_size_s = args::get(window_arg);

 if( window_size_s == "" ){
  window_size = 1000000;
 }else{
  window_size = stoi(window_size_s);
  std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
 }

 if( window_size <= 0 ) window_size = 1000000;


 std::string meta_prefix = args::get(meta_arg);
 int query_pos = args::get(query_pos_arg);

 rsq_buddy = args::get(buddy_arg);

 rsq_prune = args::get(rsq_arg);
 pval_thresh = args::get(pval_arg);

 rsq_buddy = rsq_buddy <= 0 ? 2 : rsq_buddy;
 rsq_prune = rsq_prune <= 0 ? 0.8 : rsq_prune;
 pval_thresh = pval_thresh <= 0 ? 5e-5 : pval_thresh;





 int nthreads = args::get(threads_arg);
 if( nthreads >= 1 ){
  omp_set_num_threads(nthreads);
  Eigen::setNbThreads(nthreads);
 }
 std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";

 global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_step, t_hom, t_het, t_acat, stepwise_marginal_thresh);

 if( prefix == "" ){
  if( meta_svar || meta_stepwise ){
   restore_cursor();
   std::cerr << "Error: Output prefix not specified. Try --help to see options.\n";
   return 1;
  }
 }

 if( meta_prefix == "" ){
  restore_cursor();
  std::cerr << "Error: No summary statistics files specified. Try --help to see options.\n";
  return 1;
 }
 hide_cursor();

 std::vector<std::string> meta_prefixes = split_string(meta_prefix, ',');

 if( meta_prefixes.size() <= 0 ){
  restore_cursor();
  std::cerr << "Error: No summary statistics files specified (length 0). Try --help to see options.\n";
  return 1;
 }

 cis_meta_data meta_dt(meta_prefixes, region);
# 989 "yax/src/GQT.cpp.c"
 if( meta_svar ){
  meta_dt.meta_analyze();
 }else if( meta_stepwise ){
  if( het_meta ){
   meta_dt.conditional_analysis_het();
  }else{
   meta_dt.conditional_analysis();
  }
 }

    return 0;
};


int store(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){





 args::ArgumentParser p("gqt store: Store vcov (LD; variance-covariance) data.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
 args::CompletionFlag completion(p, {"complete"});

 p.Prog(progname);

 args::Group make_sumstat_args(p, "Generating summary statistics");
  args::Flag make_vcov(make_sumstat_args, "", "Generate variance-covariance files.", {"make-vcov"});
  args::Flag make_sumstat(make_sumstat_args, "", "Generate GWAS summary statistics.", {"make-sumstats"});

 args::Group scale_args(p, "Scale and transform options");
  args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
  args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
  args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
  args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});

 args::Group input_args(p, "Input files");
  args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
  args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});

  args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for eQTL analysis.", {'b', "bed", "expression"});
  args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
  args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
  args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
# 1044 "yax/src/GQT.cpp.c"
 args::Group opt_args(p, "General options");
  args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
  args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
  args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
  args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
  args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs (LD stored in 2*size sliding window).", {'w', "window"});

  args::ValueFlag<int> blocks_arg(opt_args, "100", "Number of variants per trans-xQTL block.", {"block-size"});
  args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});





 parseModeArgs(p, beginargs, endargs);





 prefix = args::get(out_arg);
 std::string e_path = args::get(bed_arg);
 std::string g_path = args::get(bcf_arg);
 std::string c_path = args::get(cov_arg);






 std::string grm_path = args::get(grm_arg);
 std::string kin_path = args::get(kin_arg);
 double grm_scale = 1.00;

 if( kin_path != "" ){
  if( grm_path != "" ){
   std::cerr << "ERROR: Specify --kin or --grm, but not both.\n";
   abort();
  }
  grm_path = kin_path;
  grm_scale = 2.00;
 }






 trim_gene_ids = (bool) trim_ids;

 std::string region = args::get(region_arg);
 std::string gtds = args::get(gtds_arg);



 use_ds;
 if( gtds == "" || gtds == "GT" ){
  use_ds = false;
 }else if( gtds == "DS" ){
  use_ds = true;
 }else{
  std::cerr << "Invalid specification --field " << gtds << "\n";
  std::cerr << "Valid options are \"GT\" or \"DS\". Exiting. \n";
  return 1;
 }



 std::string window_size_s = args::get(window_arg);

 if( window_size_s == "" ){
  window_size = 1000000;
 }else{
  window_size = stoi(window_size_s);
  std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
 }

 if( window_size <= 0 ) window_size = 1000000;





 use_low_mem = (bool) low_mem;

 int nthreads = args::get(threads_arg);
 if( nthreads >= 1 ){
  omp_set_num_threads(nthreads);
  Eigen::setNbThreads(nthreads);
 }
 std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";

 global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_step, t_hom, t_het, t_acat, stepwise_marginal_thresh);

 if( prefix == "" ){
  restore_cursor();
  std::cerr << "Error: Output prefix not specified. Try --help to see options.\n";
  return 0;
 }
 hide_cursor();

 int max_var = 100000000;

 int n_var = 0;

 std::vector<int> variants_per_chrom;

 genotype_data g_data;
 table c_data;
 bed_data e_data;

 std::vector<std::string> bcf_chroms = get_chroms(g_path, variants_per_chrom);

 std::vector<std::string> keep_chroms = bcf_chroms;

 std::vector<std::string> bed_chroms;

 if( e_path != "" ){
  get_chroms(e_path);
  keep_chroms = intersect_ids(bcf_chroms,bed_chroms);
 }

 for(int i = 0; i < bcf_chroms.size(); i++ ){
  if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
   n_var += variants_per_chrom[i];
  }
 }


 for(const auto& c : keep_chroms){
  std::cerr << c << ",";
 }
 if( e_path != "" ){
  std::cerr << "\b present in both bcf and bed file.\n";
 }else{
  std::cerr << "\b present in bcf file.\n";
 }
 std::cerr << n_var << " total variants on selected chromosomes.\n\n";

 bcf_srs_t *sr = bcf_sr_init();

 if( region != "" ){

  std::cerr << "Setting region to " << region << " in bcf file ... \n";
  bcf_sr_set_regions(sr, region.c_str(), 0);

 }else{
  std::string region_string = "";
  for(std::string& chr : keep_chroms){
   region_string += ( region_string=="" ? "" : "," ) + chr;
  }
  bcf_sr_set_regions(sr, region_string.c_str(), 0);
 }

 bcf_sr_add_reader(sr, g_path.c_str());
 bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);


 g_data.read_bcf_header(hdr);
 std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";


 if( c_path == "" ){
  std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
  std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
  std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
  std::cerr << "    --rankNormal-resid for trait residuals.\n";
 }else{
  c_data.readHeader(c_path.c_str());
  std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
 }

 if( e_path != "" ){

  e_data.readBedHeader(e_path.c_str());
  std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";
 }

 std::vector<std::string> intersected_samples = g_data.ids.file;
 if( e_path != "" ){
  intersected_samples = intersect_ids(intersected_samples, e_data.ids.file);
 }
 if( c_path != "" ){
  intersected_samples = intersect_ids(intersected_samples, c_data.cols.file);
 }



 std::vector<std::string> intersected_samples_gto = g_data.ids.file;
 for(int i = 0; i < intersected_samples_gto.size(); ){
  if( has_element(intersected_samples, intersected_samples_gto[i]) ){
   i++;
  }else{
   intersected_samples_gto.erase(intersected_samples_gto.begin() + i);
  }
 }

 intersected_samples = intersected_samples_gto;

 std::cerr << "Found " << intersected_samples.size() << " samples in common across all files.\n\n";


 g_data.ids.setKeepIDs(intersected_samples);
 g_data.n_samples = intersected_samples.size();

 if( e_path != "" ) e_data.ids.setKeepIDs(intersected_samples);
 if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);

 std::vector<std::string> keep_regions = keep_chroms;
 if( region != "" ){
  keep_regions.clear();
  keep_regions.push_back(region);
 }



 if( c_path == "" ){
   c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
 }else{
  c_data.readFile(c_path.c_str());
  std::cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";

  if( !no_scale_x ){
   scale_and_center(c_data.data_matrix);
  }
  appendInterceptColumn(c_data.data_matrix);
 }

 if( e_path != "" ){
  e_data.readBedFile(e_path.c_str(),keep_regions);
 }
 std::cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";



 g_data.read_bcf_variants(sr, hdr, n_var, !low_mem, !low_mem);

 if( g_data.chr.size() == 0 ){
  std::cerr << "\nNo variants present in specified region(s). Exiting.\n\n";
  return 0;
 }

 if( low_mem ){
  clear_line_cerr();
  std::cerr << "Processed variant data for " << n_var << " variants.\n\n";
  g_data.genotypes.resize(0,0);
 }else{
  std::cerr << "\rFreezing genotype data for " << n_var << " variants ... \r";
  g_data.freeze_genotypes();
  clear_line_cerr();
  std::cerr << "Processed genotype data for " << n_var << " variants.\n\n";
 }

 block_intervals bm;

 std::cerr << "Generating LD blocks ... ";
 bm.make_blocks(e_data, g_data, window_size);
 std::cerr << "Done.\n\n";

 Eigen::MatrixXd &Y = e_data.data_matrix;
 Eigen::MatrixXd &X = c_data.data_matrix;

 Eigen::SparseMatrix<double> GRM;
 std::vector<int> relateds;
 if( grm_path != "" ){
  read_sparse_GRM(grm_path, GRM, intersected_samples, grm_scale, 3, relateds);
 }

 std::cerr << "Making variance-covariance (vcov) files ...\n";
 write_vcov_files(g_data, c_data);
 std::cerr << "\nDone!\n";

    return 0;
};
