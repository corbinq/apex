

# APEX: linear mixed model preprocessing guide
This page describes linear mixed model preprocessing using APEX. Once installed, you can quickly get started by running  `./apex lmm --help`. <br />

## Overview
cis-xQTL analysis in APEX uses either a) ordinary least squares (OLS) for unrelated samples or b) a linear mixed model (LMM) to account for familial relatedness or technical and biological variation.  For the former, we use a sparse genetic relatedness matrix (GRM), and for the latter, we use a low-rank random effect matrix. `./apex lmm` expedites LMM analysis for modes `cis` and `trans` by precomputing and storing a) LMM null models and trait residuals and b) spline terms for LMM genotypic variances.  Note that this preprocessing step is optional; these steps can also be performed internally using modes `cis` and `trans`; however, preprocessing can substantially reduce computation time, particularly for trans analysis. <br />

##### Table of Contents  

 1. [Precomputing LMM null models and trait residuals](#precomputing-lmm-null-models-and-trait-residuals)   
 2. [Precomputing LMM genotypic variances](#precomputing-lmm-genotypic-variances)
 3. [Command line options](#command-line-arguments) <br />

 [*Return to APEX main page.*](/apex/)

## Precomputing LMM null models and trait residuals 
**Example command:** <br />
 `./apex lmm --rankNormal --fit-null --save-resid --bcf {bcf} --bed {trait-file} --cov {covariate-file} --grm {grm} --prefix {out-name}` <br />
<br />
This command will rank-normalize traits, estimate LMM null models and residuals, and save variance component estimates and trait residuals for later use in modes `cis` and `trans`. 
<br />

## Precomputing LMM genotypic variances
**Example command:** <br />
 `./apex lmm --write-gvar --bcf {bcf} --bed {trait-file} --cov {covariate-file} --grm {grm} --prefix {out-name}` <br />
<br />
This command will calculate and store interpolation points for LMM genotype residual variances to expedite later LMM association analysis in modes `cis` and `trans`. <br />

## Command line arguments
A partial list of options is given below.  Please run `./apex lmm --help` to see a complete list of command line flags and options. 

 - **Analysis options**
	 - `--fit-null` :  Estimate LMM null models and save variance component estimates.
	 - `--save-resid` :  Save LMM trait residuals.
	 - `--write-gvar` :  Calculate and store interpolation terms for LMM genotypic variances.
 - **Output options**
	 - `--prefix`, `-o` :  Output file prefix.
	 - `--long`, `-l` :  Write cis-eQTL results in long-table format.
 - **Scale and transform options**
	 - `--rankNormal` :  Apply rank normal transform to trait values.
 - **Computational resources** 
	 - `--threads {N}` : No. threads to be used (not to exceed no. available cores).
	 - `--low-mem` : Reduce memory usage by reading and processing genotypes in chunks.  
 - **Filtering regions**
	 - `--region {chr:start-end} OR {chr}` : Only analysis variants and traits within specified region. 

