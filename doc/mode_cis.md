
# YAX: cis-xQTL analysis guide
This page describes cis-xQTL analysis using YAX. Once installed, you can quickly get started by running  ` ./yax cis --help`. <br />

## Overview
cis-xQTL analysis in YAX uses either a) ordinary least squares (OLS) for unrelated samples or b) a linear mixed model (LMM) to account for cryptic or familial relatedness fit by restricted maximum likelihood (REML).  For OLS, YAX requires 3 input files: molecular trait data, technical covariate data, and genotype data. For LMM, YAX additionally requires a kinship or genetic relatedness matrix (GRM). For detailed descriptions of input file formats, please see the [input file documentation page](https://github.com/corbinq/yax/blob/master/doc/input_files.md). <br />

##### Table of Contents  
[Vanilla cis-xQTL analysis (no related samples)](#ols-cis-xqtl-analysis-with-unrelated-samples)  
[LMM cis-xQTL analysis](#lmm-cis-xqtl-analysis)  
[Command line options](#command-line-arguments)  
<br />
[*Return to YAX main page.*](https://github.com/corbinq/yax)

## OLS cis-xQTL analysis with unrelated samples
**Example command:** <br />
 `./yax cis --vcf {vcf} --bed {expression-file} --cov {covariate-file} --prefix {output-prefix}` <br />
 <br />
**QTL software concordance.** When no GRM is specified, YAX single-variant output is equivalent to the R regression model `lm(traits[,j] ~ covariates + genotype[,k])` for each trait `j` and genotype `k`. YAX output is additionally equivalent to [FastQTL](http://fastqtl.sourceforge.net/) single-variant output.  Note that some tools, such as [QTLtools](https://qtltools.github.io/qtltools/), instead fit the model `lm(residuals[,j] ~ genotype[,k])` where `residuals[,j] = resid(lm(traits[,j] ~ covariates))`. YAX can mimic this model if the flag `--no-resid-geno` is specified.  This approach is slightly faster that standard OLS, but can cause [conservative p-values (loss of statistical power)](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22325). 
## LMM cis-xQTL analysis 
**Example command:** <br />
 `./yax cis --vcf {vcf} --bed {expression-file} --cov {covariate-file} --grm {grm-file} --prefix {output-prefix}` <br />
<br />
YAX uses a linear mixed model (LMM) to account for cryptic or familial relatedness in cis-eQTL analysis. To use this feature, specify a genetic relatedness matrix (GRM) file to YAX using  `--grm {grm-file}`.  Output files and options are otherwise similar to those from OLS cis-xQTL analysis (when `--grm` is not specified). <br />
 **LMM software concordance.** YAX's LMM estimates are consistent with the R packages [GMMAT](https://github.com/hanchenphd/GMMAT) and [GENESIS](http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) using AI-REML. 

## Command line arguments
A partial list of options is given below.  Please run `./yax cis --help` to see a complete list of command line flags and options. 
 - **General options**
	  - `--window {BP}`, `-w {BP}` : Window size in base pairs for cis-xQTL analysis.  Only variant-trait pairs within BP upstream or downstream of trait TSS will be analyzed (default: 1Mb, or `1000000`). 
 - **Output options**
	  - `--prefix`, `-o` :  Output file prefix.
	 - `--long`, `-l` :  Write cis-eQTL results in long-table format.
 -  **Scale and transform options**
	 - `--rankNormal` :  Apply rank normal transform to trait values.
	 - `--rankNormal-resid` :  Apply rank normal transform to residuals (can be used with rankNormal). [Not compatible with LMM].
	 - `--no-resid-geno` :  Do not residualize genotypes (not recommended). Output using this flag is concordant with QTLtools and some other tools. 
 - **Computational resources** 
	 - `--threads {N}` : No. threads to be used (not to exceed no. available cores).
	 - `--low-mem` : Reduce memory usage by reading and processing genotypes in chunks.  
 -  **Subsetting samples**
	 - `--exclude-iids {LIST}` : Comma-delimited list of sample IDs to exclude. 
	 - `--include-iids {LIST}` : Only include the specified comma-delimited sample IDs. 
 -  **Filtering regions and variants**
	 - `--region {chr:start-end}` : Only analysis variants and traits within specified region. 
	 - `--gene {LIST}` : Only analyze the specified comma-delimited molecular traits IDs. 
	 - `--exclude-snps {LIST}` : Comma-delimited list of SNPs to exclude. 
	 - `--include-snps {LIST}` : Only include the specified comma-delimited SNPs. 
