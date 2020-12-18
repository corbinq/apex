
# APEX: factor analysis guide
This page describes factor analysis using APEX. Once installed, you can quickly get started by running  `./apex factor --help`. <br />

## Overview
The command `apex factor`  provides factor analysis and principal components analysis to infer latent technical and biological variables from molecular trait data.  The inferred covariates can be included in cis and trans xQTL analysis as fixed-effect or random-effect covariates.   `apex factor`  accepts expression, covariate, and genotype input files; however,  expression is used in factor analysis.  Output files from `apex factor`  include inferred factor covariates together with observed covariates (if specified via `--cov {FILE}`) in a format suitable for association analysis using APEX.  For detailed descriptions of input file formats, please see the [input file documentation page](/apex/doc/input_files/). <br />

##### Table of Contents  
  1. [Running factor analysis or PCA in APEX](#factor-analysis-of-molecular-traits)
  2. [Command line options](#command-line-arguments) <br />

 [*Return to APEX main page.*](/apex/)

## Factor analysis of molecular traits
**Example command:** <br />
 `./apex factor --bed {expression-file} --cov {covariate-file} --prefix {output-prefix} --iter 3` <br />
 <br />
The above command runs factor analysis as described in our manuscript with 3 iterations.  To use PCA rather than factor analysis, specify  `--iter 0` (0 iterations).  We do not recommend attempting to iterate until convergence. 

## Command line arguments
A partial list of options is given below.  Please run `./apex factor --help` to see a complete list of command line flags and options. 
 - **General options**
	  - `--iter {N}` : Number of iterations (default: 3).
	  - `--prior-p {D}` : Prior weight (default 0, no prior).
	  - `--prior-tau {D}` : Prior variance (default 1). 
 - **Output options**
	  - `--prefix`, `-o` :  Output file prefix.
 -  **Scale and transform options**
	 - `--rankNormal` :  Apply rank normal transform to trait values.
 - **Computational resources** 
	 - `--threads {N}` : No. threads to be used (not to exceed no. available cores).
