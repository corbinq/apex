

# APEX: cis-xQTL analysis guide
This page describes cis-xQTL analysis using APEX. Once installed, you can quickly get started by running  `./apex cis --help`. <br />

## Overview
cis-xQTL analysis in APEX uses either a) ordinary least squares (OLS) or b) a linear mixed model (LMM) fit by restricted maximum likelihood (REML).  For OLS, APEX requires 3 input files: molecular trait data, technical covariate data, and genotype data.  LMM can be used to account for either cryptic familial relatedness using a kinship or genetic relatedness matrix (GRM), or to account for technical and biological variation using a low-rank matrix of random-effect covariates. For detailed descriptions of input file formats, please see the [input file documentation page](/apex/doc/input_files/). <br />

##### Table of Contents  

 1. [OLS cis-xQTL analysis](#ols-cis-xqtl-analysis)   
 2. [LMM cis-xQTL analysis with a GRM](#lmm-cis-xqtl-analysis)
 3. [Command line options](#command-line-arguments) <br />

 [*Return to APEX main page.*](/apex/)

## OLS cis-xQTL analysis
**Example command:** <br />
 `./apex cis --vcf {vcf} --bed {trait-file} --cov {covariate-file} --prefix {out-name} --long` <br />

 **Output files.** The above command generates 5 output files, `{out-name}.cis_sumstats.txt.gz`, `{out-name}.cis_sumstats.txt.gz.csi`, `{out-name}.cis_gene_table.txt.gz`, `{out-name}.cis_gene_table.txt.gz.csi`, `{out-name}.cis_long_table.txt.gz`.  The `cis_sumstats` output file contains association score statistics in a condensed format, which can be used for downstream analysis with the command `./apex meta`.  Human-readable output files are described below: <br />
 
`*.cis_long_table.tsv.gz` (flag `--long`) columns: 
 1. `#chrom` : Variant chromosome.
 2. `pos` : Variant chromosomal position (basepairs).
 3. `ref` : Variant reference allele (`A`, `C`, `T`, or `G`).
 4. `alt` : Variant alternate allele. 
 5. `gene` : Molecular trait identifier (as specified in `--bed {trait-file}`).
 6. `beta` : OLS regression slope for variant on trait. 
 7. `se` : Standard error of regression slope.
 8. `pval` : Single-variant association nominal p-value.  

`*.cis_gene_table.tsv.gz` columns:
 1. `#chrom` : Molecular trait chromosome.
 2. `start` : Molecular trait start position.
 3. `end` : Molecular trait end position.
 4. `gene` : Molecular trait identifier.
 5. `gene_pval` : Trait-level p-value calculated across all variants in the *cis* region using the [Cauchy combination test](https://arxiv.org/abs/1808.09011), comparable to beta-approximated permutation p-values. 
 6. `n_samples` : Number of samples included in analysis.
 7. `n_covar` : Number of covariates included in analysis, including intercept.
 8. `resid_sd` : Square root of regression mean squared error under the null model. 
 9. `n_cis_variants` : Number of variants in the *cis* region (which were used to calculate `gene_pval`). 

**QTL software concordance.** When no GRM is specified, APEX single-variant output is numerically equivalent to the R regression model `lm(traits[,j] ~ covariates + genotype[,k])` for each trait `j` and genotype `k`. APEX output is additionally equivalent to [FastQTL](http://fastqtl.sourceforge.net/) single-variant output.  Note that some tools, such as [QTLtools](https://qtltools.github.io/qtltools/), instead fit the model `lm(residuals[,j] ~ genotype[,k])` where `residuals[,j] = resid(lm(traits[,j] ~ covariates))`. APEX can mimic this model if the flag `--no-resid-geno` is specified.  This approach is slightly faster than standard OLS, but can cause [conservative p-values (loss of statistical power)](https://onlinelibrary.wiley.com/doi/abs/10.1002/gepi.22325). 
## LMM cis-xQTL analysis 
**Example command:** <br />
 `./apex cis --vcf {vcf} --bed {expression-file} --cov {covariate-file} --grm {grm-file} --prefix {out-name}` <br />
<br />
Here, APEX uses a linear mixed model (LMM) to account for cryptic or familial relatedness in cis-eQTL analysis of the form <img src="https://render.githubusercontent.com/render/math?math=y = X\beta %2B g %2B \varepsilon "> where <img src="https://render.githubusercontent.com/render/math?math=g\sim\mathcal{N}(0,\tau^{2}\GRM)"> and <img src="https://render.githubusercontent.com/render/math?math=\varepsilon\sim\mathcal{N}(0,\sigma^{2}I)">. To use this feature, specify a genetic relatedness matrix (GRM) file to APEX using  `--grm {grm-file}`.  Output files and options are otherwise similar to those from OLS cis-xQTL analysis (when `--grm` is not specified). [See here](/doc/input_files/) for accepted input file formats. <br />
 
**Example command:** <br />
 `./apex cis --vcf {vcf} --bed {trait-file} --cov {covariate-file} --prefix {out-name} --long` <br />

 **Output files.** Output files from LMM analysis are broadly similar to OLS. One additional output file, `{out-name}.theta.gz`, contains variance component parameter estimates from the LMM. The first 4 columns of this file list trait chromosomal position and identifier, and columns 5-7 list the residual variance component estimate <img src="https://render.githubusercontent.com/render/math?math=\sigma^2"> (independent error variance), heritable variance component estimate <img src="https://render.githubusercontent.com/render/math?math=\tau^2">, and their ratio <img src="https://render.githubusercontent.com/render/math?math=\phi=\tau^2/\sigma^2">. 
 6. Genetic variance component estimate (due to GRM). 
 7. Residual-genetic variance ratio.
 **LMM software concordance.** APEX's LMM estimates are consistent (nearly numerically equivalent) with the R packages [GMMAT](https://github.com/hanchenphd/GMMAT) and [GENESIS](http://www.bioconductor.org/packages/release/bioc/html/GENESIS.html) using AI-REML. 

## Command line arguments
A partial list of options is given below.  Please run `./apex cis --help` to see a complete list of command line flags and options. 
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
