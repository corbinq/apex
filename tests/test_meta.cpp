#include <gtest/gtest.h>
#include "setOptions.hpp"
#include "metaAnalysis.hpp"
#include <string>
#include <vector>
#include "meta_stepwise_reader.h"
using namespace std;

// Demonstrate some basic assertions.
TEST(MetaTest, ConditionalTest) {
  global_opts::reset();
  global_opts::set_global_region("");
  global_opts::set_exp_weight(0);
  global_opts::set_max_signals(10);
  global_opts::process_global_opts(
    "data/test_output.pval0",           // prefix
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
    true,          // stepwise_marginal_thresh
    true           // write out log p-values
  );

  // TODO: need actual files for this test case, currently using dummy data
  vector<string> meta_prefixes = {"data/pval0"};
  string region = "";
  cis_meta_data meta_dt(meta_prefixes, region);
  meta_dt.conditional_analysis();
  auto reader_truth = StepwiseReader("data/pval0.cis_meta.stepwise.tsv"); // wrong file but just scratch testing for now
  auto reader_test = StepwiseReader("data/test_output.pval0.cis_meta.stepwise.tsv");
  ASSERT_TRUE(reader_truth == reader_test);
}