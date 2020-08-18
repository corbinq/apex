#include "GQT.hpp"

using namespace std;


const static Eigen::IOFormat EigenCSV(Eigen::StreamPrecision, Eigen::DontAlignCols, ",", "\n");
const static Eigen::IOFormat EigenTSV(Eigen::StreamPrecision, Eigen::DontAlignCols, "\t", "\n");


int main(int argc, char* argv[])
{

	// Argument parsing using https://github.com/Taywee/args
	
	args::ArgumentParser p("gqt: GWAS/QTL Toolkit.", "Contact: corbinq@gmail.com.\n");
    args::HelpFlag help(p, "help", "Display this help menu", {'h', "help"});

	args::Group make_sumstat_args(p, "Generating summary statistics");
		args::Flag make_vcov(make_sumstat_args, "", "Generate variance-covariance files.", {"make-vcov"});
		args::Flag make_sumstat(make_sumstat_args, "", "Generate GWAS summary statistics.", {"make-sumstats"});
	
	args::Group cis_args(p, "cis-eQTL analysis");
		args::ValueFlag<string> test_arg(cis_args, "het,alt", "Meta-analysis test forms (comma-separated combination of 'het,hom,alt').", {'t', "tests"});
		args::Flag make_cis_sumstat(cis_args, "", "Generate cis-eQTL summary statistics.", {"make-cis-sumstats"});
		args::Flag make_long(cis_args, "", "Write cis-eQTL results in long table format.", {'l', "long"});
		args::Flag just_long(cis_args, "", "Only write long-table cis-eQTL results.", {'L', "just-long"});
	
	args::Group scale_args(p, "Scale and transform options");
		args::Flag rknorm_y(scale_args, "", "Apply rank normal transform to trait values.", {"rankNormal"});
		args::Flag rknorm_r(scale_args, "", "Apply rank normal transform to residuals (can be used with rankNormal).", {"rankNormal-resid"});
		args::Flag no_scale_x(scale_args, "", "Do not scale and center covariates (otherwise done by default).", {"no-scale-cov"});
		args::Flag no_resid_geno(scale_args, "", "Do not residualize genotypes (not recommended).", { "no-resid-geno"});
	
	args::Group input_args(p, "Input files");
		args::ValueFlag<string> bcf_arg(input_args, "", "Genotype file path (vcf, vcf.gz, or bcf format).", {'v', "vcf", "bcf"});
		args::ValueFlag<string> cov_arg(input_args, "", "Covariate/trait file path.", { 'c', "cov"});
		// args::ValueFlag<string> trait_arg(input_args, "", "Trait file path.", {'t', "trait-file"});
		args::ValueFlag<string> bed_arg(input_args, "", "Expression file path for eQTL analysis.", {'b', "bed", "expression"});
		args::ValueFlag<string> grm_arg(input_args, "", "Sparse GRM file.", {"grm"});
		args::ValueFlag<string> kin_arg(input_args, "", "Sparse kinship file.", {"kin"});
		args::ValueFlag<string> gtds_arg(input_args, "", "Genotype field (\"GT\" by default, or \"DS\" for imputed dosages).", {"field"});
	
	/*
	args::Group subset_args(p, "Subsetting samples");
		args::ValueFlag<string> iid_e_arg(subset_args, "", "List of samples to exclude (file path or comma-separated).", {"exclude-iids"});
		args::ValueFlag<string> iid_i_arg(subset_args, "", "Only include specified samples (file path or comma-separated).", {"include-iids"});
	
	args::Group filter_args(p, "Filtering variants");
		args::ValueFlag<string> iid_e_arg(subset_args, "", "List of variants to exclude (file path).", {"exclude-snps"});
		args::ValueFlag<string> iid_i_arg(subset_args, "", "Only include specified variants (file path ).", {"include-snps"});
	*/
	
	args::Group opt_args(p, "General options");
		args::Flag low_mem(opt_args, "", "Lower memory usage.", {"low-mem"});
		args::ValueFlag<string> out_arg(opt_args, "", "Prefix for output files.", {'o', "prefix", "out"});
		args::ValueFlag<string> region_arg(opt_args, "", "Subset to specified genomic region.", {'r', "region"});
		args::ValueFlag<string> window_arg(opt_args, "1000000", "Window size in base pairs for cis-eQTL or gene-based analysis.", {'w', "window"});
		args::ValueFlag<string> gene_arg(opt_args, "", "Restrict analysis to specified gene(s).", {"gene"});
		// args::ValueFlag<string> ld_window_arg(opt_args, "1000000", "Window size in base pairs for LD files.", {'w', "window"});
		args::ValueFlag<double> rsq_arg(opt_args, "", "Maximum rsq threshold.", {"rsq"});
		args::ValueFlag<double> buddy_arg(opt_args, "", "Print LD buddies (specify rsq threshold).", {"buddies"});
		args::Flag use_marginal_pval(opt_args, "", "Apply threshold to marginal rather than ACAT p-value in stepwise algorithm (no adjustment for no. variants).", {"marginal"});
		args::ValueFlag<double> pval_arg(opt_args, "", "P-value threshold for stepwise procedures.", {"pvalue"});
		//args::ValueFlag<int> blocks_arg(opt_args, "2500", "Number of variants per block for cis-eQTL or gene-based analysis.", {"block-size"});
		args::Flag trim_ids(opt_args, "", "Trim version numbers from Ensembl gene IDs.", {"trim-ids"});
		args::Flag backward_step(opt_args, "", "Perform backward step in stepwise selection.", {"backward"});
	
	args::Group meta_args(p, "Analyzing summary statistics");
		args::ValueFlag<string> meta_arg(meta_args, "", "Sumstat file prefixes for meta-analysis (comma-separated).", {"sumstats"});
		args::Flag meta_svar(meta_args, "", "Perform single-variant meta-analysis.", {"meta"});
		args::Flag meta_stepwise(meta_args, "", "Perform stepwise (conditional) meta-analysis.", {"stepwise"});
		args::Flag ivw_1(meta_args, "", "Calculate IVW weights under the alternative hypothesis.", {"ivw1"});
		args::Flag het_meta(meta_args, "", "Allow heterogeneous genotype effects across studies in stepwise meta-analysis.", {"het"});
	
	args::Group vcov_args(p, "Querying vcov files");
		args::Flag query_vcov(vcov_args, "", "Query meta-analysis vcov data (file list from --sumstats). Can be used with --region.", {"query"});
		args::ValueFlag<int> query_pos_arg(vcov_args, "pos", "Only get vcov data for specified variant (against all others in region).", {"pos"});
		args::Flag query_gtg(vcov_args, "", "Just return uncentered vcov data from bin file (default: centered v-cov matrix).", {"just-bin"});
		//args::Flag query_colnames(vcov_args, "", "Show SNP IDs as column names of the (default: centered v-cov matrix).", {"just-bin"});
	
	args::ValueFlag<long> start_arg(p, "", "", {'s'});
	args::ValueFlag<long> n_arg(p, "", "", {'n'});
	
	
	// Hide the cursor while printing 
	// logs to cerr, and restor on exit

	auto lam_kill = 
      [] (int i) { restore_cursor(); cerr << "\nKilled.\n" << endl; exit(0); };
	
    signal(SIGINT, lam_kill);
    signal(SIGABRT, lam_kill);
    signal(SIGTERM, lam_kill);
    // signal(SIGTSTP, lam_wait);
	
	at_quick_exit (restore_cursor);
	atexit (restore_cursor);
	
	
	// Parsing the command line options
	
    try
    {
        p.ParseCLI(argc, argv);
    }
    catch (args::Help)
    {
		restore_cursor();
        cout << p;
        return 0;
    }
    catch (args::ParseError e)
    {
		restore_cursor();
        cerr << e.what() << endl;
        cerr << p;
        return 1;
    }
    catch (args::ValidationError e)
    {
		restore_cursor();
        cerr << e.what() << endl;
        cerr << p;
        return 1;
    }
	
	// ----------------------------------
	// Process test options
	// ----------------------------------
	
	bool t_hom = false, t_het = false, t_acat = false;
	string tests = args::get(test_arg);
	
	if( tests == "" ){
		tests = "hom,acat";
	}
	for( char& c_i : tests){
		c_i = tolower(c_i);
	}
	for( const string& test_i : split_string(tests, ',') ){
		if( test_i == "hom" ){
			t_hom = true;
		}else if( test_i == "het" ){
			t_het = true;
		}else if( test_i == "acat" ){
			t_acat = true;
		}else{
			cerr << "Error: Unknown test '" << test_i << "'.\n";
			cerr << "Valid options are het,hom,acat (comma-separated combinations)\n";
			abort();
		}
	}
	if( !t_hom && !t_het && !t_acat  ){
		cerr << "Error: No valid test options specified.\n";
		cerr << "Valid options are het,hom,acat (comma-separated combinations)\n";
		abort();
	}
	
	// ----------------------------------
	// I/O File Paths
	// ----------------------------------
	
	string prefix = args::get(out_arg);
	string e_path = args::get(bed_arg);
	string g_path = args::get(bcf_arg);
	string c_path = args::get(cov_arg);
	
	
	// ----------------------------------
	// GRM paths and options
	// ----------------------------------
	
	string grm_path = args::get(grm_arg);
	string kin_path = args::get(kin_arg);
	double grm_scale = 1.00;
	
	if( kin_path != "" ){
		if( grm_path != "" ){
			cerr << "ERROR: Specify --kin or --grm, but not both.\n";
			abort();
		}
		grm_path = kin_path;
		grm_scale = 2.00; 
	}
	
	
	// ----------------------------------
	// Input subsetting: Regions, genotype fields, target genes
	// ----------------------------------
	
	string region = args::get(region_arg);
	string gtds = args::get(gtds_arg);
	vector<string> target_genes = split_string(args::get(gene_arg), ',');
	
	// Use imputed dosages rather than genotype hard calls
	bool use_ds;
	if( gtds == "" || gtds == "GT" ){
		use_ds = false;
	}else if( gtds == "DS" ){
		use_ds = true;
	}else{
		cerr << "Invalid specification --field " << gtds << "\n";
		cerr << "Valid options are \"GT\" or \"DS\". Exiting. \n";
		return 1;
	}
	
	//cerr << "\n" << region << "\n\n";
	
	string window_size_s = args::get(window_arg);
	
	int window_size;
	
	if( window_size_s == "" ){
		window_size = 1000000;
	}else{
		window_size = stoi(window_size_s);
		cerr << "Set window size to " << window_size/1000000 << " Mbp.\n";
	}
	
	if( window_size <= 0 ) window_size = 1000000;
	
	//string query_prefix = args::get(query_arg);
	string meta_prefix = args::get(meta_arg);
	int query_pos = args::get(query_pos_arg);
	
	double rsq_buddy = args::get(buddy_arg);
	
	double rsq_prune = args::get(rsq_arg);
	double pval_thresh = args::get(pval_arg);

	rsq_buddy = rsq_buddy <= 0 ? 2 : rsq_buddy;
	rsq_prune = rsq_prune <= 0 ? 0.8 : rsq_prune;
	pval_thresh = pval_thresh <= 0 ? 5e-5 : pval_thresh;
	
	// ----------------------------------
	// Set global options
	// ----------------------------------
	
	global_opts::process_global_opts(prefix, low_mem, rsq_buddy, rsq_prune, pval_thresh, window_size, target_genes, (bool) ivw_1, use_ds, (bool) trim_ids, (bool) backward_step, t_hom, t_het, t_acat, (bool) use_marginal_pval);
	
	if( prefix == "" && meta_prefix == "" )
	{
		restore_cursor();
	}else
	{
		hide_cursor();
	}
	
	long start = args::get(start_arg);
	long num = args::get(n_arg);
	
	if( meta_prefix != "" ){
		
		if( prefix == "" ){
			if( meta_svar || meta_stepwise ){
				cerr << "Error: Output prefix not specified. Try --help to see options.\n";
				return 0;
			}
		}
		
		vector<string> meta_prefixes = split_string(meta_prefix, ',');
		
		if( meta_prefixes.size() > 0 ){
			
			cis_meta_data test(meta_prefixes, region);
			
			if( meta_svar ){
				test.meta_analyze();
			}else if( meta_stepwise ){
				if( het_meta ){
					test.conditional_analysis_het();
				}else{
					test.conditional_analysis();
				}
			}else if( query_vcov ){
				
				Eigen::MatrixXd test_out;
				
				start = 0;
				num = test.vc.pos.size();
				vector<int> qi;
				for(int i = 0; i < num; i++){
					qi.push_back(i);
				}
				
				if( query_pos <= 0 ){
					
					for(int i = 0; i < num; i++){
						if( i > 0 ) cout << "\t";
						cout << test.vc.chr[i] << "_" << test.vc.pos[i] << "_" << test.vc.ref[i] << "_" << test.vc.alt[i];
					}
					cout << "\n";
					
					if( query_gtg ){
						test_out = test.vc.getGtG(qi);
					}else{
						test_out = test.vc.getV(qi);
					}
					
					cout << test_out.format(EigenTSV) << endl;
				
				}else{
					
					lindex ldx(test.vc.pos);
					int query_index = ldx.index(query_pos);
					
					if( test.vc.pos[query_index] != query_pos ){
						cerr << "No match found for pos=" << query_pos << " (closest=" << test.vc.pos[query_index] << "\n";
						return 0;
					}
					
					cout << "ID\t" << test.vc.chr[query_index] << "_" << test.vc.pos[query_index] << "_" << test.vc.ref[query_index] << "_" << test.vc.alt[query_index] << "\n";
					
					if( query_gtg ){
						test_out = test.vc.getGtG(vector<int>(1,query_index),qi);
					}else{
						test_out = test.vc.getV(vector<int>(1,query_index),qi);
					}
					
					for(int i = 0; i < test_out.rows(); i++){
						double val = test_out(i,0);
						cout << test.vc.chr[i] << "_" << test.vc.pos[i] << "_" << test.vc.ref[i] << "_" << test.vc.alt[i];
						cout << "\t" << val << "\n";
					}
				}
				
				
			}
		}
		
		return 0; 
	}
	
	if( prefix == "" ){
		cerr << "Error: Output prefix not specified. Try --help to see options.\n";
		return 0;
	}
	
	int max_var = 100000000;
	
	int n_var = 0;
	
	vector<int> variants_per_chrom;
	
	genotype_data g_data;
	table c_data;
	bed_data e_data;
	
	vector<string> bcf_chroms = get_chroms(g_path, variants_per_chrom);

	for(int i = 0; i < bcf_chroms.size(); ++i ){
		cout << bcf_chroms[i] << "\t" << variants_per_chrom[i] << "\n";
	}
	
	vector<string> bed_chroms = get_chroms(e_path);
	
	for(int i = 0; i < bed_chroms.size(); ++i ){
		cout << bed_chroms[i] << "\n";
	}
	
	vector<string> keep_chroms = intersect_ids(bcf_chroms,bed_chroms);
	
	for(int i = 0; i < bcf_chroms.size(); ++i ){
		if( find(keep_chroms.begin(), keep_chroms.end(), bcf_chroms[i]) != keep_chroms.end() ){
			n_var += variants_per_chrom[i];
		}
	}
	
	// Show chromosomes present across files.
	for(const auto& c : keep_chroms){
		cerr << c << ",";
	}
	cerr << "\b present in both bcf and bed file.\n";
	cerr << n_var << " total variants on selected chromosomes.\n\n";
	
	bcf_srs_t *sr = bcf_sr_init();
	
	if( region != "" ){
		
		cerr << "Setting region to " << region << " in bcf file ... \n";
		bcf_sr_set_regions(sr, region.c_str(), 0);
		
	}else{
		string region_string = "";
		for(string& chr : keep_chroms){   
			region_string += ( region_string=="" ? "" : "," ) + chr;
		}
		bcf_sr_set_regions(sr, region_string.c_str(), 0);
	}
	
	bcf_sr_add_reader(sr, g_path.c_str());
	bcf_hdr_t *hdr = bcf_sr_get_header(sr, 0);
	
	// read header from bcf file
	g_data.read_bcf_header(hdr);
	cerr << "Found " << g_data.ids.file.size() << " samples in bcf file ... \n";
	
	// read header from covariate file
	if( c_path == "" ){
		cerr << "\nWARNING: No covariate file specified. That's usually a bad idea.\n";
		cerr << "    Covariates can be specified using --cov FILE. Use --rankNormal\n";
		cerr << "    to rank-normal (aka, inverse-normal) transform traits, and use\n";
		cerr << "    --rankNormal-resid for trait residuals.\n";
	}else{
		c_data.readHeader(c_path.c_str());
		cerr << "Found " << c_data.cols.file.size() << " samples in covariate file ... \n";
	}
	
	// read header from expression bed file
	e_data.readBedHeader(e_path.c_str());
	cerr << "Found " << e_data.ids.file.size() << " samples in expression bed file ... \n";
	
	vector<string> intersected_samples;
	if( c_path == "" ){
		intersected_samples = intersect_ids(g_data.ids.file, e_data.ids.file);
	}else{
		intersected_samples = intersect_ids(intersect_ids(g_data.ids.file, c_data.cols.file), e_data.ids.file);
	}
	
	// order of intersected samples should match genotype file
	
	vector<string> intersected_samples_gto = g_data.ids.file;
	for(int i = 0; i < intersected_samples_gto.size(); ){
		if( has_element(intersected_samples, intersected_samples_gto[i]) ){
			i++;
		}else{
			intersected_samples_gto.erase(intersected_samples_gto.begin() + i);
		}
	}
	
	intersected_samples = intersected_samples_gto;
	
	cerr << "Found " << intersected_samples.size() << " samples in common across all three files.\n\n";
	
	// set to the intersection across all three files
	g_data.ids.setKeepIDs(intersected_samples);
	g_data.n_samples = intersected_samples.size();
	
	e_data.ids.setKeepIDs(intersected_samples);
	if( c_path != "" ) c_data.cols.setKeepIDs(intersected_samples);
	
	vector<string> keep_regions = keep_chroms;
	if( region != "" ){
		keep_regions.clear();
		keep_regions.push_back(region);
	} 
	
	// now let's read the expression and covariate data .. 
	
	if( c_path == "" ){
		 c_data.data_matrix = Eigen::MatrixXd::Constant(intersected_samples.size(), 1, 1.0);
	}else{
		c_data.readFile(c_path.c_str());
		cerr << "Processed data for " << c_data.data_matrix.cols() << " covariates across " << c_data.data_matrix.rows() << " samples.\n";
			
		if( !no_scale_x ){
			scale_and_center(c_data.data_matrix);
		}
		appendInterceptColumn(c_data.data_matrix);
	}
	
	e_data.readBedFile(e_path.c_str(),keep_regions);
	cerr << "Processed expression for " << e_data.data_matrix.cols() << " genes across " << e_data.data_matrix.rows() << " samples.\n";
	
	// let's get the genotypes. 
	
	g_data.read_bcf_variants(sr, hdr, n_var, !low_mem, !low_mem);
	
	if( g_data.chr.size() == 0 ){
		cerr << "\n\nWhat the hell, man?!  No variants present in specified region(s). Exiting.\n\n";
		return 0;
	}
	
	if( low_mem ){
		clear_line_cerr();
		cerr << "Processed variant data for " << n_var << " variants.\n\n";
		g_data.genotypes.resize(0,0);
	}else{
		cerr << "\rFreezing genotype data for " << n_var << " variants ... \r";
		g_data.freeze_genotypes();
		clear_line_cerr();
		cerr << "Processed genotype data for " << n_var << " variants.\n\n";
	}
	
	cerr << "Generating eVariant-eGene blocks ... ";
	block_intervals bm(e_data, g_data, window_size);
	cerr << "Done.\n\n";
	
	Eigen::MatrixXd &Y = e_data.data_matrix;
	Eigen::MatrixXd &X = c_data.data_matrix;
	
	Eigen::SparseMatrix<double> GRM;
	if( grm_path != "" ){
		read_sparse_GRM(grm_path, GRM, intersected_samples, grm_scale, 3);
	}
	
	if( make_sumstat || make_long || just_long ){
		if( grm_path == "" ){
			run_cis_eQTL_analysis(sr, hdr, g_data, c_data, e_data, bm, rknorm_y, rknorm_r, make_sumstat, make_long, just_long);
		}else{
			run_cis_eQTL_analysis_LMM(sr, hdr, g_data, c_data, e_data, GRM, bm, rknorm_y, rknorm_r, make_sumstat, make_long, just_long);
		}
	}
	
	if( make_vcov ){
		cerr << "Making variance-covariance (vcov) files ...\n";
		write_vcov_files(g_data, c_data);
		cerr << "\nDone!\n";
	}
	
    return 0;
}

