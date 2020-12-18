

# APEX: variance-covariance / LD storage guide
This page describes how to store variance-covariance matrices, which capture covariate-adjusted linkage disequilibrium (LD), in APEX. These files are required for multiple-variant analysis from summary statistics (without individual-level data), including meta-analysis.  Once installed, you can quickly get started by running `./apex store --help`. <br />

## Overview
We recommend running `apex store` separately on each chromosome, which can be accomplished by specifying `--region chr1`.  Each chromosomal vcov file is indexed by chromosomal position and byte offset, allowing fast access to variance-covariance data within a specific region or a specific variant.  [See here](/apex/doc/benchmarking/#meta-analysis) for information on expected vcov file sizes. <br />

`apex store` should be run with the same set of input files, output prefix, and other options used to store association summary statistics with `apex cis`.  This ensures that both vcov and sumstat files include the same set of individuals, and the same set of covariates, and that `apex meta` and `apexR` are able to easily link sumstat and vcov files.  <br />

Note that APEX vcov files capture covariated-adjusted LD (the covariance of genotype residuals), which can differ substantially from unadjusted LD in structured samples.  Results from `apex meta` and the `apexR` R package are designed to be (nearly) numerically equivalent to results from individual-level data -- this requires covariate-adjusted (rather than raw) LD.  For detailed descriptions of input file formats, please see the [input file documentation page](/apex/doc/input_files/). <br />

##### Table of Contents  

 1. [Running APEX Store](#running-apex-store)
 2. [Command line arguments](#command-line-arguments)

 [*Return to APEX main page.*](/apex/)

## Running APEX Store
**Example command:** <br />
 `./apex store --vcf {vcf} --bed {trait-file} --cov {covariate-file} --prefix {out-name}` <br />

 **Output files.** The above command generates 3 output files, `{out-name}.vcov.bin`, `{out-name}.vcov.idx.gz`, `{out-name}.cis_long_table.tsv.gz`.  These files store LD data, LD index and covariate adjustment terms, and tabix index respectively. In general, these files are not intended to be human-readable, but can be queried and used for multiple-variant analysis using `apex meta` and `apexR`. <br />

## Command line arguments
A partial list of options is given below.  Please run `./apex store --help` to see a complete list of command line flags and options. 
 - **General options**
	  - `--window {BP}`, `-w {BP}` : Window size in base pairs for cis-xQTL analysis.  LD will be stored within twice the specified window size. 
 - **Computational resources** 
	 - `--threads {N}` : No. threads to be used (not to exceed no. available cores).
	 - `--low-mem` : Reduce memory usage by reading and processing genotypes in chunks.  
