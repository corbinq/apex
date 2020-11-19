/*  YAX:

    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of YAX.

    YAX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of YAX.
*/

#include "Main.hpp"
#include <string>

// ------------------------------------
//  Default global options. 
// ------------------------------------

// TO_DO_PL: re-organize global options.
// TO_DO_PL: split modes into separate source files.
// TO_DO_PM: sample and variant filters.

std::string prefix = "";
bool use_low_mem = false; 
double rsq_buddy = 2.0;
double rsq_prune = 0.80;
double pval_thresh = 5e-5;
int n_ePCs = 0;
int max_signals = 10;
int window_size = 1000000;
std::vector<std::string> target_genes;
bool use_ivw_1 = false;
bool use_ds = false;
bool trim_gene_ids = false;
double stepwise_backward_thresh = 1.00;
bool t_hom = false;
bool t_het = false;
bool t_acat = false;
bool stepwise_marginal_thresh = false;

// ------------------------------------
//  Analysis modes
// ------------------------------------

int cis(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int trans(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int factor(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int meta(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int store(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);

using mode_fun = std::function<int(const std::string &, std::vector<std::string>::const_iterator, std::vector<std::string>::const_iterator)>;

std::string help_string = 
"\n" 
"  YAX: Toolkit for xQTL analysis\n" 
"     (c) 2019-2020 Corbin Quick and Li Guan.\n" 
"\n" 
"  Usage and options:\n" 
"     ./yax [mode] --help       Print help menu for [mode].\n" 
"\n" 
"  Analysis modes:\n" 
"     ./yax cis {OPTIONS}       Run cis-xQTL analysis.\n" 
"\n" 
"     ./yax trans {OPTIONS}     Run trans-xQTL analysis.\n" 
"\n"
"     ./yax factor {OPTIONS}    Estimate latent factors from\n" 
"                                 molecular trait data.\n" 
"\n" 
"     ./yax meta {OPTIONS}      Single and multi-variant\n" 
"                                 xQTL meta-analysis from\n" 
"                                 sumstat and vcov files.\n" 
"\n" 
"     ./yax store {OPTIONS}     Store vcov (LD) data for\n" 
"                                 xQTL meta-analysis or\n" 
"                                 data sharing.\n" 
"\n" 
"  Contact: corbinq@gmail.com\n";


// ------------------------------------
//  Wrapper for parsing within-mode command line options 
// ------------------------------------

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


// ------------------------------------
//  Main function (determine running mode)
// ------------------------------------

int main(int argc, char* argv[])
{

// ------------------------------------
//  Hide cursor while running, restore on exit or interrupt
// ------------------------------------

	auto lam_kill = 
      [] (int i) { restore_cursor(); std::cerr << "\nKilled.\n" << "\n"; exit(0); };
	
    signal(SIGINT, lam_kill);
    signal(SIGABRT, lam_kill);
    signal(SIGTERM, lam_kill);
    // signal(SIGTSTP, lam_wait);
	
	at_quick_exit (restore_cursor);
	atexit (restore_cursor);
	

// ------------------------------------
//  YAX main menu: Parse running mode.
// ------------------------------------

    std::unordered_map<std::string, mode_fun> map{
        {"cis", cis},
        {"trans", trans},
        {"meta", meta},
		{"factor", factor},
        {"store", store}
	};
	
	// Argument parsing using https://github.com/Taywee/args
	args::ArgumentParser p0("yax: GWAS/QTL Toolkit.", "Contact: corbinq@gmail.com.\n");
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
	
	args::ArgumentParser p("yax cis: cis-xQTL analysis.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
	args::CompletionFlag completion(p, {"complete"});

	p.Prog(progname);

	args::Group analysis_args(p, "Analysis options");
		args::Flag fit_null(analysis_args, "", "Estimate and store LMM null model parameters.", {"fit-null"});
		args::Flag stepwise(analysis_args, "", "Estimate and store conditionally independent cis signal genotypes.", {"stepwise"});
		args::ValueFlag<std::string> theta_arg(analysis_args, "", "Use stored LMM null model parameters.", {"theta-file"});
		args::ValueFlag<int> max_mod_arg(analysis_args, "", "Maximum model size in stepwise regression.", {"max-model"});
		args::ValueFlag<double> dtss_arg(analysis_args, "", "dTSS weight for eGene p-values.", {"dtss-weight"});
		args::ValueFlag<int> epc_arg(analysis_args, "", "Number of ePCs extracted from trait matrix.", {"epcs"});
		args::ValueFlag<int> epc_fe_arg(analysis_args, "", "Number of ePCs used as fixed effect covariates.", {"fe-epcs"});
		args::Flag use_egrm(analysis_args, "", "Use eGRM rather than fixed-effect ePCs.", {"use-egrm"});
		args::ValueFlag<std::string> loco_arg(analysis_args, "", "Leave-one-chr-out (LOCO) to calculate ePCs or eGRMs.", {"loco"});
		
	args::Group cis_args(p, "Output options");
		args::Flag make_long(cis_args, "", "Write cis-QTL results in long table format.", {'l', "long"});
		args::Flag just_long(cis_args, "", "Only write long-table cis-QTL results.", {'L', "just-long"});
	
	args::Group scale_args(p, "Scale and transform options");
		args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
		args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
		// args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
		args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});
	
	args::Group input_args(p, "Input files");
		args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
		args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});
		// args::ValueFlag<std::string> trait_arg(input_args, "", "Trait file path.", {'t', "trait-file"});
		args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for QTL analysis.", {'b', "bed", "expression"});
		args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
		args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
		args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
	
	/*
	args::Group subset_args(p, "Subsetting samples");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of samples to exclude (file path or comma-separated).", {"exclude-iids"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified samples (file path or comma-separated).", {"include-iids"});
	
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
		args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
		args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
		args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
		args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-QTL or gene-based analysis.", {'w', "window"});
		args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});
		// args::ValueFlag<std::string> ld_window_arg(opt_args, "1000000", "Window size in base pairs for LD files.", {'w', "window"});
		args::ValueFlag<double> pval_arg(opt_args, "", "P-value threshold for stepwise procedures.", {"pvalue"});
		args::ValueFlag<double> rsq_arg(opt_args, "", "Rsq threshold for variable selection.", {"rsq"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
	
	// ----------------------------------
	// Parse command line arguments 
	// ----------------------------------
	
	parseModeArgs(p, beginargs, endargs);
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	std::string theta_path = args::get(theta_arg);
	
	prefix = args::get(out_arg);
	std::string e_path = args::get(bed_arg);
	std::string g_path = args::get(bcf_arg);
	std::string c_path = args::get(cov_arg);
	
	
	// ----------------------------------
	// GRM paths and options
	// ----------------------------------
	
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
	
	
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
	double dtss_w = args::get(dtss_arg);
	
	global_opts::set_exp_weight(dtss_w);
	
	std::string loco = args::get(loco_arg);
	n_ePCs = args::get(epc_arg);
	int n_ePCs_FE = args::get(epc_fe_arg);
	
	trim_gene_ids = (bool) trim_ids;
	
	std::string region = args::get(region_arg);
	std::string gtds = args::get(gtds_arg);
	target_genes = split_string(args::get(gene_arg), ',');
	
	// Use imputed dosages rather than genotype hard calls
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
	
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	int max_signals_new = args::get(max_mod_arg);
	if( max_signals_new > 0 ){
		max_signals = max_signals_new;
	}
	global_opts::set_max_signals(max_signals);
	
	double pval_thresh_new = args::get(pval_arg);
	
	if( pval_thresh_new > 0 ){
		pval_thresh = pval_thresh_new;
	}
	
	double rsq_new = args::get(rsq_arg);
	
	if( rsq_new > 0 ){
		rsq_buddy = rsq_new;
		rsq_prune = rsq_new;
	}
	
	use_low_mem = (bool) low_mem;
	
	int nthreads = args::get(threads_arg);
	if( nthreads >= 1 ){
		omp_set_num_threads(nthreads);
		Eigen::setNbThreads(nthreads);
	}
	std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";
	
	global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_thresh, t_hom, t_het, t_acat, stepwise_marginal_thresh);
	
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
	
	// Show chromosomes present across files.
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
	
	// read header from bcf file
	g_data.read_bcf_header(hdr);
	std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";
	
	// read header from covariate file
	if( c_path == "" ){
		std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
		std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
		std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
		std::cerr << "    --rankNormal-resid for trait residuals.\n";
	}else{
		c_data.readHeader(c_path.c_str());
		std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
	}
	
	// read header from expression bed file
	e_data.readBedHeader(e_path.c_str());
	std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";
	
	std::vector<std::string> intersected_samples;
	if( c_path == "" ){
		intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
	}else{
		intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
	}
	
	// order of intersected samples should match genotype file
	
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
	
	// set to the intersection across all three files
	g_data.ids.setKeepIDs(intersected_samples);
	g_data.n_samples = intersected_samples.size();
	
	e_data.ids.setKeepIDs(intersected_samples);
	if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);
	
	std::vector<std::string> keep_regions = keep_chroms;
	if( region != "" ){
		keep_regions.clear();
		keep_regions.push_back(region);
	} 
	
	// now let's read the expression and covariate data .. 
	
	if( c_path == "" ){
		c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
	}else{
		c_data.readFile(c_path.c_str());
		std::cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";
			
		// if( !no_scale_x ){
			scale_and_center(c_data.data_matrix);
		// }
		appendInterceptColumn(c_data.data_matrix);
	}
		
	std::vector<int> test_idx;
	std::vector<int> epc_idx;
	
	bool egrm_loco = false;
	
	//if( n_ePCs > 0 ){
	//	e_data.readBedFile(e_path.c_str(),bed_chroms);
	//}else{
		e_data.readBedFile(e_path.c_str(),keep_regions);
	//}
	
	std::cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";
	
	// let's get the genotypes. 
	
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
	
	if( stepwise ){
		if( grm_path != "" ){
			std::cerr << "Error: GRM not currently supported for cis_signal.\n";
			return 1;
		}
		scan_signals(sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r);
		return 0;
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
		if( n_ePCs > 0 ){
			std::cerr << "Using " << n_ePCs << " ePCs.\n";
			std::cerr << "Reading full trait matrix to construct ePCs...\n";
			bed_data r_data;
			r_data.readBedHeader(e_path.c_str());
			r_data.ids.setKeepIDs(intersected_samples);
			
			std::vector<std::string> r_chroms = bed_chroms;
			
			std::vector<int> rm_chroms;
			for( int i = r_chroms.size() - 1; i >= 0; i--){
				if(r_chroms[i] == loco){
					rm_chroms.push_back(i);
				}
			}
			for( const int& i : rm_chroms){
				r_chroms.erase(r_chroms.begin()+i);
			}
			
			r_data.readBedFile(e_path.c_str(),r_chroms);
			
			std::cerr << "Processed expression for " << r_data.data_matrix.cols() << " genes across " << r_data.data_matrix.rows() << " samples.\n";
			
			if ( use_egrm ){
				run_cis_QTL_analysis_eLMM(n_ePCs, n_ePCs_FE, sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r, true, make_long, just_long, r_data.data_matrix);
			}else{
				run_cis_QTL_analysis_eFE(n_ePCs, sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r, true, make_long, just_long, r_data.data_matrix);
			}
		}else{
			std::cerr << "Using OLS (no GRM specified).\n";
			run_cis_QTL_analysis(sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r, true, make_long, just_long);
		}
	}else{
		if( n_ePCs > 0 ){
			std::cerr << "Error: GRM cannot currently be included with random ePCs.\n";
			return 1;
		}else{
			run_cis_QTL_analysis_LMM(sr, hdr, g_data, c_data, e_data, GRM, relateds, bm, rknorm_y, rknorm_r, true, make_long, just_long);
		}
	}
	
    return 0;
};


int trans(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){
	
	// ----------------------------------
	// Define command line flags
	// ----------------------------------
	
	args::ArgumentParser p("yax trans: trans-xQTL analysis.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
	args::CompletionFlag completion(p, {"complete"});
	
	p.Prog(progname);

	args::Group analysis_args(p, "Analysis options");
		args::Flag fit_null(analysis_args, "", "Estimate and store LMM null model parameters.", {"fit-null"});
		args::Flag save_resid(analysis_args, "", "Estimate and store LMM null model residuals.", {"save-resid"});
		args::ValueFlag<std::string> theta_arg(analysis_args, "", "Use stored LMM null model parameters.", {"theta-file"});
		args::ValueFlag<std::string> cis_arg(analysis_args, "", "Use stored cis signals as covariates in trans analysis.", {"cis-file"});
		args::ValueFlag<double> fpr_arg(analysis_args, "", "Nominal false positive rate (p-value threshold).", {"fpr"});

	args::Group scale_args(p, "Scale and transform options");
		args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
		args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
		args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
		args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});
	
	args::Group input_args(p, "Input files");
		args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
		args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});
		// args::ValueFlag<std::string> trait_arg(input_args, "", "Trait file path.", {'t', "trait-file"});
		args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for QTL analysis.", {'b', "bed", "expression"});
		args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
		args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
		args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
	
	/*
	args::Group subset_args(p, "Subsetting samples");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of samples to exclude (file path or comma-separated).", {"exclude-iids"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified samples (file path or comma-separated).", {"include-iids"});
	
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
		args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
		args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
		args::Flag sloppy(opt_args, "", "Use sloppy covariate adjustment (faster, but less powerful).", {"sloppy"});
		args::Flag write_anchors(opt_args, "", "Save var interpolation points.", {"write-var"});
		args::ValueFlag<std::string> region_arg(opt_args, "", "Subset genotypes to specified genomic region.", {'r', "region"});
		args::ValueFlag<std::string> bed_region_arg(opt_args, "", "Subset bed to specified genomic region.", {"bed-region"});
		args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-QTL or gene-based analysis.", {'w', "window"});
		args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});
		args::ValueFlag<int> blocks_arg(opt_args, "100", "Number of variants per trans-xQTL block.", {"block-size"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
	
	// ----------------------------------
	// Parse command line arguments 
	// ----------------------------------
	
	parseModeArgs(p, beginargs, endargs);
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	std::string theta_path = args::get(theta_arg);
	std::string cis_path = args::get(cis_arg);
	
	prefix = args::get(out_arg);
	std::string e_path = args::get(bed_arg);
	std::string g_path = args::get(bcf_arg);
	std::string c_path = args::get(cov_arg);
	
	global_opts::save_residuals( (bool) save_resid);
	
	
	// ----------------------------------
	// GRM paths and options
	// ----------------------------------
	
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
	
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
	trim_gene_ids = (bool) trim_ids;
	
	std::string region = args::get(region_arg);
	std::string bed_region = args::get(bed_region_arg);
	std::string gtds = args::get(gtds_arg);
	target_genes = split_string(args::get(gene_arg), ',');
	
	int block_size = args::get(blocks_arg);
	
	if( block_size <= 0 ){
		block_size = 100;
	}
	
	// Use imputed dosages rather than genotype hard calls
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
	
	//cerr << "\n" << region << "\n\n";
	
	std::string window_size_s = args::get(window_arg);
	
	if( window_size_s == "" ){
		window_size = 1000000;
	}else{
		window_size = stoi(window_size_s);
		std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
	}
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	double fpr = args::get(fpr_arg);
	
	if( fpr > 0 ){
		pval_thresh = fpr;
	}
	
	use_low_mem = (bool) low_mem;
	
	global_opts::set_lmm_options( (bool) write_anchors );
	
	int nthreads = args::get(threads_arg);
	if( nthreads >= 1 ){
		omp_set_num_threads(nthreads);
		Eigen::setNbThreads(nthreads);
	}
	std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";
	
	global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_thresh, t_hom, t_het, t_acat, stepwise_marginal_thresh);
	
	if( sloppy ){
		global_opts::use_sloppy_covar();
	}
	
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
	// std::vector<std::string> keep_chroms = intersect_ids(bcf_chroms,bed_chroms);
	
	std::vector<std::string> keep_chroms = bcf_chroms;
	
	for(int i = 0; i < bcf_chroms.size(); i++ ){
		if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
			n_var += variants_per_chrom[i];
		}
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
	
	// read header from bcf file
	g_data.read_bcf_header(hdr);
	std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";
	
	// read header from covariate file
	if( c_path == "" ){
		std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
		std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
		std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
		std::cerr << "    --rankNormal-resid for trait residuals.\n";
	}else{
		c_data.readHeader(c_path.c_str());
		std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
	}
	
	// read header from expression bed file
	e_data.readBedHeader(e_path.c_str());
	std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";
	
	std::vector<std::string> intersected_samples;
	if( c_path == "" ){
		intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
	}else{
		intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
	}
	
	// order of intersected samples should match genotype file
	
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
	
	// set to the intersection across all three files
	g_data.ids.setKeepIDs(intersected_samples);
	g_data.n_samples = intersected_samples.size();
	
	e_data.ids.setKeepIDs(intersected_samples);
	if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);
	
	std::vector<std::string> keep_regions = keep_chroms;
	if( region != "" ){
		keep_regions.clear();
		keep_regions.push_back(region);
	} 
	
	// now let's read the expression and covariate data .. 
	
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
	
	if( bed_region == "" ){
		e_data.readBedFile(e_path.c_str(), bed_chroms);
	}else{
		keep_regions.clear();
		keep_regions.push_back(bed_region);
		e_data.readBedFile(e_path.c_str(), keep_regions);
		keep_regions.clear();
	}

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
		run_trans_QTL_analysis(sr, hdr, g_data, c_data, e_data, rknorm_y, rknorm_r, make_sumstat, make_long, block_size, cis_path, bed_region);
	}else{
		if( cis_path != "" ){
			std::cerr << "Error: Options --grm {grm} and --cis-file {file} and currently incompatible.\n";
			return 1;
		}
		run_trans_QTL_analysis_LMM(sr, hdr, g_data, c_data, e_data, GRM, relateds, rknorm_y, rknorm_r, make_sumstat, make_long, block_size, theta_path);
	}
	
    return 0;
};


int factor(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){
	
	args::ArgumentParser p("yax factor: high-dimensional factor analysis.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});
	args::CompletionFlag completion(p, {"complete"});

	p.Prog(progname);

	args::Group analysis_args(p, "Analysis options");
		args::ValueFlag<int> epc_arg(analysis_args, "", "Number of latent common factors.", {"factors"});
		args::ValueFlag<std::string> iter_arg(analysis_args, "", "Number of factor analysis iterations (0 for PCA).", {"iter"});
		args::ValueFlag<std::string> pp_arg(analysis_args, "", "Factor analysis prior p.", {"prior-p"});
		args::ValueFlag<std::string> pt_arg(analysis_args, "", "Factor analysis prior tau.", {"prior-tau"});

	args::Group scale_args(p, "Scale and transform options");
		args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
		args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});

	args::Group input_args(p, "Input files");
		args::ValueFlag<std::string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
		args::ValueFlag<std::string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});
		// args::ValueFlag<std::string> trait_arg(input_args, "", "Trait file path.", {'t', "trait-file"});
		args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for QTL analysis.", {'b', "bed", "expression"});
	
	/*
	args::Group subset_args(p, "Subsetting samples");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of samples to exclude (file path or comma-separated).", {"exclude-iids"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified samples (file path or comma-separated).", {"include-iids"});
	
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
		args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
		args::ValueFlag<std::string> region_arg(opt_args, "", "Subset genotypes to specified genomic region.", {'r', "region"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
	
	// ----------------------------------
	// Parse command line arguments 
	// ----------------------------------
	
	parseModeArgs(p, beginargs, endargs);
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	prefix = args::get(out_arg);
	std::string e_path = args::get(bed_arg);
	std::string g_path = args::get(bcf_arg);
	std::string c_path = args::get(cov_arg);
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
	n_ePCs = args::get(epc_arg);
	int n_FA_iter = 3;
	double fa_p = 0.001;
	double fa_t = 1.00;
	std::string n_FA_iter_s = args::get(iter_arg);
	std::string fa_p_s = args::get(pp_arg);
	std::string fa_t_s = args::get(pt_arg);
	
	if( n_FA_iter_s != "" ){
		n_FA_iter = std::stoi(n_FA_iter_s);
	}
	if( fa_p_s != "" ){
		fa_p = std::stod(fa_p_s);
	}
	if( fa_t_s != "" ){
		fa_t = std::stod(fa_t_s);
	}
	
	global_opts::set_factor_par(n_FA_iter, fa_p, fa_t);
	
	std::string region = args::get(region_arg);
	trim_gene_ids = (bool) trim_ids;
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	int nthreads = args::get(threads_arg);
	if( nthreads >= 1 ){
		omp_set_num_threads(nthreads);
		Eigen::setNbThreads(nthreads);
	}
	std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";
	
	global_opts::process_global_opts(prefix, true, 0.8, 0.8, 0.8, 1000000, std::vector<std::string>(0), true, false, trim_gene_ids, 0.8, true, true, true, 0.8);
	
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
	
	// Show chromosomes present across files.
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
	
	bcf_hdr_t *hdr;
	
	if( g_path != "" ){
		
		bcf_sr_add_reader(sr, g_path.c_str());
		hdr = bcf_sr_get_header(sr, 0);
		
		// read header from bcf file
		g_data.read_bcf_header(hdr);
		std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";	
			
	}
	
	// read header from covariate file
	if( c_path == "" ){
		std::cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
		std::cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
		std::cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
		std::cerr << "    --rankNormal-resid for trait residuals.\n";
	}else{
		c_data.readHeader(c_path.c_str());
		std::cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
	}
	
	// read header from expression bed file
	e_data.readBedHeader(e_path.c_str());
	std::cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";
	
	std::vector<std::string> intersected_samples;
	if( c_path == "" && g_path == "" ){
		intersected_samples = e_data.ids.file;
	}else if( c_path == "" ){
		intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
	}else{
		intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
	}
	
	if( g_path != "" ){
		
		// order of intersected samples should match genotype file
		
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

		// set to the intersection across all three files
		g_data.ids.setKeepIDs(intersected_samples);
		g_data.n_samples = intersected_samples.size();
		
	}
	
	e_data.ids.setKeepIDs(intersected_samples);
	if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);
	
	std::vector<std::string> keep_regions = keep_chroms;
	if( region != "" ){
		keep_regions.clear();
		keep_regions.push_back(region);
	} 
	
	// now let's read the expression and covariate data .. 
	
	if( c_path == "" ){
		c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
	}else{
		c_data.readFile(c_path.c_str());
		std::cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";
			
		// if( !no_scale_x ){
			scale_and_center(c_data.data_matrix);
		// }
		appendInterceptColumn(c_data.data_matrix);
	}
		
	std::vector<int> test_idx;
	std::vector<int> epc_idx;
	
	//if( n_ePCs > 0 ){
	//	e_data.readBedFile(e_path.c_str(),bed_chroms);
	//}else{
		e_data.readBedFile(e_path.c_str(),keep_regions);
	//}
	
	std::cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";
	
	// let's get the genotypes. 
	
	
	save_fa_covariates(n_ePCs, g_data, c_data, e_data, rknorm_y, rknorm_r);
	
    return 0;
};


int meta(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs){

	// ----------------------------------
	// Define command line flags
	// ----------------------------------
	
	args::ArgumentParser p("yax meta: Meta-analysis of xQTL studies.", "Contact: corbinq@gmail.com.\n");
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
		args::ValueFlag<double> dtss_arg(analysis_args, "", "dTSS weight for eGene p-values.", {"dtss-weight"});
		args::ValueFlag<std::string> test_arg(analysis_args, "het,alt", "Meta-analysis test forms (comma-separated combination of 'het,hom,alt').", {'t', "tests"});
		args::Flag ivw_1(analysis_args, "", "Calculate IVW weights under the alternative hypothesis.", {"ivw1"});
		args::Flag het_meta(analysis_args, "", "Allow heterogeneous genotype effects across studies in stepwise meta-analysis.", {"het"});
		args::ValueFlag<int> max_mod_arg(analysis_args, "", "Maximum model size in stepwise regression.", {"max-model"});
		args::ValueFlag<double> rsq_arg(analysis_args, "", "Maximum rsq threshold.", {"rsq"});
		args::ValueFlag<double> buddy_arg(analysis_args, "", "Print LD buddies (specify rsq threshold).", {"buddies"});
		args::Flag use_marginal_pval(analysis_args, "", "Apply threshold to marginal rather than ACAT p-value in stepwise algorithm (no adjustment for no. variants).", {"marginal"});
		args::ValueFlag<double> pval_arg(analysis_args, "", "P-value threshold for stepwise procedures.", {"pvalue"});
		args::ValueFlag<double> backward_arg(analysis_args, "", "Backward step p-value threshold in stepwise selection.", {"backward"});
	/*
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
		args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
		args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-QTL or gene-based analysis.", {'w', "window"});
		args::ValueFlag<std::string> gene_arg(opt_args, "", "Restrict analysis to specified genes (gene name or comma-separated list).", {"gene"});
		// args::ValueFlag<std::string> ld_window_arg(opt_args, "1000000", "Window size in base pairs for LD files.", {'w', "window"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
	
	args::Group vcov_args(p, "Querying vcov files");
		args::Flag query_vcov(vcov_args, "", "Query meta-analysis vcov data (file list from --sumstats). Can be used with --region.", {"query"});
		args::ValueFlag<int> query_pos_arg(vcov_args, "pos", "Only get vcov data for specified variant (against all others in region).", {"pos"});
		args::Flag query_u_vcov(vcov_args, "", "Just return uncentered vcov data from bin file (default: centered v-cov matrix).", {"just-bin"});
		//args::Flag query_colnames(vcov_args, "", "Show SNP IDs as column names of the (default: centered v-cov matrix).", {"just-bin"});
	
	
	// ----------------------------------
	// Parse command line arguments 
	// ----------------------------------
	
	parseModeArgs(p, beginargs, endargs);
	
	// ----------------------------------
	// Process test options
	// ----------------------------------
	
	use_ivw_1 = (bool) ivw_1;
	// stepwise_backward_step = (bool) backward_step;
	
	stepwise_backward_thresh = args::get(backward_arg);
	
	if(stepwise_backward_thresh <= 0){
		stepwise_backward_thresh = 1.00;
	}
	
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
	if( !t_hom && !t_het && !t_acat  ){
		std::cerr << "Error: No valid test options specified.\n";
		std::cerr << "Valid options are het,hom,acat (comma-separated combinations)\n";
		abort();
	}
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	prefix = args::get(out_arg);
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
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
	
	//string query_prefix = args::get(query_arg);
	std::string meta_prefix = args::get(meta_arg);
	int query_pos = args::get(query_pos_arg);
	
	rsq_buddy = args::get(buddy_arg);
	
	rsq_prune = args::get(rsq_arg);
	pval_thresh = args::get(pval_arg);

	rsq_buddy = rsq_buddy <= 0 ? 2 : rsq_buddy;
	rsq_prune = rsq_prune <= 0 ? 0.8 : rsq_prune;
	pval_thresh = pval_thresh <= 0 ? 5e-5 : pval_thresh;
	
	double dtss_w = args::get(dtss_arg);
	global_opts::set_exp_weight(dtss_w);
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	int max_signals_new = args::get(max_mod_arg);
	if( max_signals_new > 0 ){
		max_signals = max_signals_new;
	}
	global_opts::set_max_signals(max_signals);
	
	int nthreads = args::get(threads_arg);
	if( nthreads >= 1 ){
		omp_set_num_threads(nthreads);
		Eigen::setNbThreads(nthreads);
	}
	std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";
	
	global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_thresh, t_hom, t_het, t_acat, stepwise_marginal_thresh);
	
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
	/*
	for(const std::string& tg : target_genes){
		// std::cout << stoi(tg) << "\n";
		// std::cout << query_u_vcov << "\n";
		meta_dt.get_vcov_gene((int) stoi(tg), (bool) !query_u_vcov);
	}
	return 0;
	*/
	
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
	
	// ----------------------------------
	// Define command line flags
	// ----------------------------------
	
	args::ArgumentParser p("yax store: Store vcov (LD; variance-covariance) data.", "Contact: corbinq@gmail.com.\n");
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
		// args::ValueFlag<std::string> trait_arg(input_args, "", "Trait file path.", {'t', "trait-file"});
		args::ValueFlag<std::string> bed_arg(input_args, "", "Expression file path for QTL analysis.", {'b', "bed", "expression"});
		args::ValueFlag<std::string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
		args::ValueFlag<std::string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
		args::ValueFlag<std::string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
	
	/*
	args::Group subset_args(p, "Subsetting samples");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of samples to exclude (file path or comma-separated).", {"exclude-iids"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified samples (file path or comma-separated).", {"include-iids"});
	
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<std::string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<std::string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::ValueFlag<int> threads_arg(opt_args, "", "No. threads (not to exceed no. available cores).", {"threads"});
		args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
		args::ValueFlag<std::string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
		args::ValueFlag<std::string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
		args::ValueFlag<std::string> window_arg(opt_args, "1000000", "Window size in base pairs (LD stored in 2*size sliding window).", {'w', "window"});
		// args::ValueFlag<std::string> ld_window_arg(opt_args, "1000000", "Window size in base pairs for LD files.", {'w', "window"});
		args::ValueFlag<int> blocks_arg(opt_args, "100", "Number of variants per trans-xQTL block.", {"block-size"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
	
	// ----------------------------------
	// Parse command line arguments 
	// ----------------------------------
	
	parseModeArgs(p, beginargs, endargs);
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	prefix = args::get(out_arg);
	std::string e_path = args::get(bed_arg);
	std::string g_path = args::get(bcf_arg);
	std::string c_path = args::get(cov_arg);
	
	
	// ----------------------------------
	// GRM paths and options
	// ----------------------------------
	
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
	
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
	trim_gene_ids = (bool) trim_ids;
	
	std::string region = args::get(region_arg);
	std::string gtds = args::get(gtds_arg);
	// target_genes = split_string(args::get(gene_arg), ',');
	
	// Use imputed dosages rather than genotype hard calls
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
	
	//cerr << "\n" << region << "\n\n";
	
	std::string window_size_s = args::get(window_arg);
	
	if( window_size_s == "" ){
		window_size = 1000000;
	}else{
		window_size = stoi(window_size_s);
		std::cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
	}
	
	if( window_size <= 0 ) window_size = 1000000;
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	use_low_mem = (bool) low_mem;
	
	int nthreads = args::get(threads_arg);
	if( nthreads >= 1 ){
		omp_set_num_threads(nthreads);
		Eigen::setNbThreads(nthreads);
	}
	std::cerr << "Using " << Eigen::nbThreads() << " threads.\n";
	
	global_opts::process_global_opts(prefix, use_low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, use_ivw_1, use_ds, trim_gene_ids, stepwise_backward_thresh, t_hom, t_het, t_acat, stepwise_marginal_thresh);
	
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
		bed_chroms = get_chroms(e_path);
		keep_chroms = intersect_ids(bcf_chroms,bed_chroms);
	}
	
	for(int i = 0; i < bcf_chroms.size(); i++ ){
		if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
			n_var += variants_per_chrom[i];
		}
	}
	
	// Show chromosomes present across files.
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
	
	// read header from bcf file
	g_data.read_bcf_header(hdr);
	std::cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";
	
	// read header from covariate file
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
		// read header from expression bed file
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
	
	// order of intersected samples should match genotype file
	
	std::vector<std::string> intersected_samples_gto = g_data.ids.file;
	
	std::unordered_set<std::string> iid_set(intersected_samples.begin(), intersected_samples.end());
	
	for(int i = 0; i < intersected_samples_gto.size(); ){
		if( has_element(intersected_samples, intersected_samples_gto[i]) ){
			i++;
		}else{
			intersected_samples_gto.erase(intersected_samples_gto.begin() + i);
		}
	}
	
	intersected_samples = intersected_samples_gto;
	
	std::cerr << "Found " << intersected_samples.size() << " samples in common across all files.\n\n";
	
	// set to the intersection across all three files
	g_data.ids.setKeepIDs(intersected_samples);
	g_data.n_samples = intersected_samples.size();
	
	if( e_path != "" ) e_data.ids.setKeepIDs(intersected_samples);
	if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);
	
	std::vector<std::string> keep_regions = keep_chroms;
	if( region != "" ){
		keep_regions.clear();
		keep_regions.push_back(region);
	} 
	
	// now let's read the expression and covariate data .. 
	
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
	
	// let's get the genotypes. 
	
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
	write_vcov_files(sr, hdr, g_data, c_data);
	std::cerr << "\nDone!\n";
	
    return 0;
};






