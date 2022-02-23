#include "meta_stepwise_reader.h"
#include <fstream>
#include <string>
#include <iostream>
#include <regex>
#include <optional>
#include "reader_util.h"

using namespace std;

void StepwiseReader::load(const string &file) {
  if (!filepath_exists(file)) {
    throw std::runtime_error("Could not find file: " + file);
  }

  ifstream input_file(file);
  string line;
  auto line_separator = regex("[ \t]");
  auto field_separator = regex(",");

  // Line regexes
  auto regex_header = regex("#gene\tstudies\tsignal.*");

  bool header_done = false;
  bool parse_trait = false;
  while (getline(input_file, line)) {
    smatch match;
    if (!header_done) {
      if (regex_search(line, match, regex_header)) {
        header_done = true;
      }
    }
    else {
      // Begin parsing record row
      vector<string> tokens;
      copy(sregex_token_iterator(line.begin(), line.end(), line_separator, -1), sregex_token_iterator(), back_inserter(tokens));

      if (line.substr(0,1) == "#") { continue; }

      // Create record
      auto rec = make_shared<StepwiseRecord>();
      rec->gene = tokens.at(0);
      copy(sregex_token_iterator(tokens.at(1).begin(), tokens.at(1).end(), field_separator, -1), sregex_token_iterator(), back_inserter(rec->studies));
      rec->signal = tokens.at(2);
      rec->variant = tokens.at(3);
      rec->beta = extract_fp<double>(tokens.at(4));
      rec->se = extract_fp<double>(tokens.at(5));
      rec->pval_joint = extract_fp<long double>(tokens.at(6));
      rec->pval_marginal = extract_fp<long double>(tokens.at(7));
      rec->pval_stepwise = extract_fp<long double>(tokens.at(8));

      // Keys
//      string chrpos = rec->chrom + ":" + to_string(rec->pos);
//      string variant = rec->chrom + ":" + to_string(rec->pos) + "_" + rec->ref + "/" + rec->alt;

      // Insert
      records.push_back(rec);
      this->index.emplace(rec->gene, rec);
      this->index.emplace(rec->variant, rec);
    }
  }

  if (records.empty()) {
    throw std::runtime_error("No score test results were read from file: " + file);
  }
}

StepwiseReader::StepwiseReader(const string &file) {
  load(file);
}

shared_ptr<StepwiseRecord> StepwiseReader::get_record(const string &i) {
  auto it = this->index.find(i);
  if (it != this->index.end()) {
    return it->second;
  }
  else {
    auto null = shared_ptr<StepwiseRecord>();
    return null;
  }
}

bool StepwiseReader::operator==(const StepwiseReader& other) {
  if (records.size() != other.records.size()) {
    return false;
  }

  const size_t nrecords = records.size();
  for (size_t i = 0; i < nrecords; i++) {
    auto& rec_a = records[i];
    auto& rec_b = other.records[i];

    if (rec_a->gene != rec_b->gene)
      return false;
    if (rec_a->signal != rec_b->signal)
      return false;
    if (rec_a->variant != rec_b->variant)
      return false;
    if (!approx_equal(rec_a->beta, rec_b->beta))
      return false;
    if (!approx_equal(rec_a->se, rec_b->se))
      return false;
    if (!approx_equal(rec_a->pval_joint, rec_b->pval_joint, PVALUE_MAX_REL_DIFF))
      return false;
    if (!approx_equal(rec_a->pval_signal, rec_b->pval_signal, PVALUE_MAX_REL_DIFF))
      return false;
    if (!approx_equal(rec_a->pval_marginal, rec_b->pval_marginal, PVALUE_MAX_REL_DIFF))
      return false;
    if (!approx_equal(rec_a->pval_stepwise, rec_b->pval_stepwise, PVALUE_MAX_REL_DIFF))
      return false;
  }

  return true;
}