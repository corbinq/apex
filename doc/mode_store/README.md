

# YAX: variance-covariance / LD storage guide
This page describes how to store variance-covariance matrices, which capture covariate-adjusted linkage disequilibrium (LD), in YAX. These files are required for multiple-variant analysis from summary statistics (without individual-level data), including meta-analysis.  Once installed, you can quickly get started by running ` ./yax store --help`. <br />

## Overview
We recommend running `yax store` separately on each chromosome, which can be accomplished by specifying `--region chr1`.  Each chromosomal vcov file is indexed by chromosomal position and byte offset, allowing fast access to variance-covariance data within a specific region or a specific variant.<br />

`yax store` should be run with the same set of input files, output prefix, and other options used to store association summary statistics with `yax cis`.  This ensures that both vcov and sumstat files include the same set of individuals, and the same set of covariates, and that `yax meta` and `yaxR` are able to easily link sumstat and vcov files.  <br />

Note that YAX vcov files capture covariated-adjusted LD (the covariance of genotype residuals), which can differ substantially from unadjusted LD in structured samples.  Results from `yax meta` and the `yaxR` R package are designed to be (nearly) numerically equivalent to results from individual-level data -- this requires covariate-adjusted (rather than raw) LD.  For detailed descriptions of input file formats, please see the [input file documentation page](/yax/doc/input_files/). <br />

##### Table of Contents  

 1. [Running YAX Store](#running-yax-store)
 2. [Command line arguments](#command-line-arguments)

 [*Return to YAX main page.*](/yax/)

## Running YAX Store
**Example command:** <br />
 `./yax store --vcf {vcf} --bed {trait-file} --cov {covariate-file} --prefix {out-name}` <br />

 **Output files.** The above command generates 3 output files, `{out-name}.vcov.bin`, `{out-name}.vcov.idx.gz`, `{out-name}.cis_long_table.tsv.gz`.  These files store LD data, LD index and covariate adjustment terms, and tabix index respectively. In general, these files are not intended to be human-readable, but can be queried and used for multiple-variant analysis using `yax meta` and `yaxR`. <br />

## Command line arguments
A partial list of options is given below.  Please run `./yax store --help` to see a complete list of command line flags and options. 
 - **General options**
	  - `--window {BP}`, `-w {BP}` : Window size in base pairs for cis-xQTL analysis.  LD will be stored within twice the specified window size. 
 - **Computational resources** 
	 - `--threads {N}` : No. threads to be used (not to exceed no. available cores).
	 - `--low-mem` : Reduce memory usage by reading and processing genotypes in chunks.  
