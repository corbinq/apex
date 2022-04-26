#include <gtest/gtest.h>
#include "setOptions.hpp"
#include "cli.hpp"
#include <vector>
#include "processVCOV.hpp"
using namespace std;

/**
 * Test to ensure --legacy-vcov functions to produce files in older format.
 */
TEST(StoreTest, LegacyVcovTest) {
  global_opts::reset();
  string output_prefix = "data/test_output.legacy_test.on";
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    output_prefix,           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  // Run with legacy format
  vector<string> on_args = {"--legacy-vcov", "--bed", "data/legacy_test.bed.gz", "--vcf", "data/legacy_test.vcf.gz", "--field", "DS", "--prefix", output_prefix, "--window", "100"};
  store("tests", on_args.begin(), on_args.end());

  global_opts::reset();
  output_prefix = "data/test_output.legacy_test.off";
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    output_prefix,           // prefix
    false,        // use low mem
    2,            // rsq_buddy
    0.8,          // rsq_prune
    0.05,         // p-value threshold
    1000000,      // window size
    {},           // target genes
    '0',          // ivw method
    false,        // use_ds (dosage)
    false,        // trim gene ids
    0.05,         // stepwise_backward_thresh
    true,         // t_hom
    false,        // t_het
    true,         // t_acat
    true,         // stepwise_marginal_thresh
    false         // write out log p-values
  );

  // Run with latest format
  vector<string> off_args = {"--bed", "data/legacy_test.bed.gz", "--vcf", "data/legacy_test.vcf.gz", "--field", "DS", "--prefix", output_prefix, "--window", "100"};
  store("tests", off_args.begin(), off_args.end());

  // Open vcov files
  vcov_bin_file bin_new("data/test_output.legacy_test.on.vcov.bin");
  vcov_bin_file bin_leg("data/test_output.legacy_test.off.vcov.bin");

  // Check to make sure two files differ
  // This is a pathological test case where input data was specifically crafted to cause a difference in vcov files under
  // the two different formats.
  auto byte1_new = bin_new.get_data(1, 100);
  auto byte1_leg = bin_leg.get_data(1, 100);

  ASSERT_NE(byte1_new, byte1_leg);
}