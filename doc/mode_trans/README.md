
# APEX: trans-xQTL analysis guide
This page describes trans-xQTL analysis using APEX. Once installed, you can quickly get started by running  ` ./apex trans --help`. <br />

## Overview
The command `apex trans` can be used to analyze genome-wide associations between molecular traits and all genetic variants.  This is in contrast to `apex cis`, which tests only analyzes genetic variants within a window of each moleculatr trait.  The underlying statistical methods are broadly similar between modes `cis` and `trans`; however, we introduce additional optimizations in mode `trans` to reduce computation time, memory, and storage. <br />
Similar to `apex cis` mode, trans-xQTL analysis in APEX (`apex trans`) uses either a) ordinary least squares (OLS) for unrelated samples or b) a linear mixed model (LMM) to account for cryptic or familial relatedness fit by restricted maximum likelihood (REML). For OLS, APEX requires 3 input files: molecular trait data, technical covariate data, and genotype data. For LMM, APEX additionally requires a kinship or genetic relatedness matrix (GRM). For detailed descriptions of input file formats, please see the [input file documentation page](/doc/input_files/). <br />

##### Table of Contents  
  1. [Vanilla trans-xQTL analysis (no related samples)](#ols-trans-xqtl-analysis-with-unrelated-samples)  
  2. [LMM trans-xQTL analysis](#lmm-trans-xqtl-analysis)  
  3. [Command line options](#command-line-arguments) <br />

 [*Return to APEX main page.*](https://github.com/corbinq/apex)

## OLS trans-xQTL analysis with unrelated samples
**Example command:** <br />
 `./apex trans --vcf {vcf} --bed {expression-file} --cov {covariate-file} --prefix {output-prefix}` <br />
 <br />
**QTL software concordance.** When no GRM is specified, APEX single-variant output is equivalent to the R regression model `lm(traits[,j] ~ covariates + genotype[,k])` for each trait `j` and genotype `k`. APEX output is additionally equivalent to [FastQTL](http://fastqtl.sourceforge.net/) single-variant output.  Note that some tools, such as [QTLtools](https://qtltools.github.io/qtltools/), instead fit the model `lm(residuals[,j] ~ genotype[,k])` where `residuals[,j] = resid(lm(traits[,j] ~ covariates))`. APEX can mimic this model if the flag `--no-resid-geno` is specified.  This approach is slightly faster that standard OLS, but can cause [conservative p-values (loss of statistical power)](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22325).  To see accepted input file formats, [please see here.](/doc/input_files/)
## LMM trans-xQTL analysis 
**Example command:** <br />
```
## Estimate null LMM models for all molecular traits and 
## store estimates for later use:
 ./apex trans --vcf {vcf} --bed {expression-file} --cov {covariate-file} --grm {grm-file} --fit-null --prefix {theta-prefix}
## Run trans-xQTL analysis, re-using variance component 
## estimates from the previous step:
 ./apex trans --vcf {vcf} --bed {expression-file} --cov {covariate-file} --grm {grm-file} --theta-file {theta-prefix}.theta.gz --prefix {output-prefix}
```
<br />
APEX uses a linear mixed model (LMM) to account for cryptic or familial relatedness in trans-eQTL analysis. To use this feature, specify a genetic relatedness matrix (GRM) file to APEX using  `--grm {grm-file}`. To see accepted input file formats, [please see here.](/doc/input_files/) <br />
Unlike `apex cis`, LMM analysis in `apex trans` is divided into two steps. First, we estimate variance component parameters for all molecular traits under the null hypothesis (no single-variant genetic effects), and store these estimates for later use. Second, we use these estimates to quickly calculate trans-xQTL association statistics. When jobs are parallelizes across chromosomes, this 2-step approach saves substantial computational resources, as the null model for each molecular trait need only be estimated once. <br />

 **LMM software concordance.** APEX's LMM estimates are consistent with the R packages [GMMAT](https://github.com/hanchenphd/GMMAT) and [GENESIS](http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) using AI-REML. 

## Command line arguments
A partial list of options is given below.  Please run `./apex trans --help` to see a complete list of command line flags and options. 
 - **General options**
	  - `--pvalue {P}` : Only report trans-xQTL associations with p-value <= {P}. 
 - **Output options**
	  - `--prefix`, `-o` :  Output file prefix.
	 - `--long`, `-l` :  Write trans-eQTL results in long-table format.
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
