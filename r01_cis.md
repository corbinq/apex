# YAX: cis-xQTL analysis guide
To quickly get started with cis-xQTL analysis in YAX, run  ` ./yax cis --help` to see a list of options. 
**Basic usage:**
 `./yax cis --vcf {vcf} --bed {expression-file} --cov {covariate-file} --prefix {output-prefix}` 
 **Advanced usage:**  
Use   `--grm {grm-file}` to specify a genetic relatedness matrix (GRM) file.  When specified, YAX uses a mixed model to account for cryptic or familial relatedness in cis-eQTL analysis. 
## Overview
You can use YAX for cis-xQTL analysis using either a) ordinary least squares (OLS) regression, or b) a linear mixed model (LMM) to account for cryptic or familial relatedness fit by restricted maximum likelihood (REML).  For OLS, YAX requires 3 input files: molecular trait data, technical covariate data, and genotype data. For LMM, YAX additionally requires a kinship or genetic relatedness matrix (GRM). Input file formats are generally consistent with those used in the [GTEx QTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl).  Further details on each required input file are given below.   
YAX uses the intersection of sample IDs present across trait, covariate, and genotype files for statistical analysis.  Sample IDs need not be listed in the same order across files. 
### Genotype data
YAX accepts genotype data in VCF and BCF format.  Genotype files should be indexed using [Tabix or BCFtools](http://samtools.github.io/bcftools/) (`.tbi` or `.csi` made using `bcftools index`).  Note that VCF/BCF and molecular trait files should be mapped to the same genome assembly (e.g., GRCh38); see [UCSC LiftOver](http://hgdownload.cse.ucsc.edu/downloads.html) to convert coordinates between assemblies. 
### Molecular trait data
Molecular trait data (`--bed {file}`) must be stored in bed file format, and should be compressed using BGZIP and index using Tabix. For example:

    #chr  start    end      gene_name        sample_1  sample_2
    1     65418    65419    ENSG00000186092  -0.0837   -0.3476
    1     827521   827522   ENSG00000225880   1.0369    1.3489
The 4 columns are required, where 1-3 specify the chromosomal  coordinates of each trait (for example, gene transcription start site [TSS] location) and 4 specifies the trait name or label (for example, Ensembl ID).   
### Covariate data
Covariate files (`--cov {file}`)  are stored similar to molecular traits, with 1 row per covariate and 1 column per sample. For example, 

    #ID   sample_1    sample_2    sample_3    sample_4
    PC1   0.0139      0.0145      0.0141      0.0135
    PC2  -0.0097     -0.0059     -0.0025     -0.0064
    PC3   0.0067      0.0096     -0.0036     -0.0041
Covariate files should be white-space delimited. Uncompressed or GZIP/BGZIP-compressed white-space delimited text files are supported.
### Genetic relatedness and kinship matrices
Genetic relatedness matrices (GRMs)  (`--grm {file}`) or kinship (`--kin {file}`) can be specified in a sparse matrix format as follows:

    #id1         id2            kinship
    sample_1     sample_1       0.50
    sample_1     sample_20      0.05
    sample_1     sample_22      0.15

If diagonal elements (self-kinship or inbreeding coefficients) are listed, then YAX will analyze only sample IDs present in the GRM file.  If diagonal elements are not listed, then YAX assumes that diagonal elements are equal to 0.5 for kinship matrices or 1 for GRMs, and assumes that sample not listed in the GRM file are unrelated.  Uncompressed or GZIP/BGZIP-compressed white-space delimited text files are supported.

