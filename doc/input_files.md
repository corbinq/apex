# YAX: Input file formats
This page summarizes individual-level data file formats accepted by the YAX commands `yax cis` (for cis-xQTL analysis), `yax trans` (for trans-xQTL analysis), and `yax store` (to store variance-covariance / LD information for data sharing or meta-analysis). File formats for xQTL data are generally consistent with those used in the [GTEx QTL pipeline](https://github.com/broadinstitute/gtex-pipeline/tree/master/qtl).  

YAX uses the intersection of sample IDs present across trait, covariate, and genotype files for statistical analysis.  Sample IDs need not be listed in the same order across files. Detailed descriptions of input data formats are provided below. 

##### Table of Contents  
[Genotype data](#genotype-data)  
[Molecular trait data](#molecular-trait-data)  
[Covariate data](#covariate-data)  
[Genetic relatedness and kinship matrices](#genetic-relatedness-and-kinship-matrices)  
<br />
[*Return to YAX main page.*](https://github.com/corbinq/yax)

## Genotype data
**Relevant flags:** `--bcf {FILE}`, `--vcf {FILE}`.  <br />
**File formats**. YAX accepts genotype data in VCF and BCF format.  Genotype files should be indexed using [Tabix or BCFtools](http://samtools.github.io/bcftools/) (`.tbi` or `.csi` made using `tabix` or `bcftools index`).  Note that VCF/BCF and molecular trait files should be mapped to the same genome assembly (e.g., GRCh38); see [UCSC LiftOver](http://hgdownload.cse.ucsc.edu/downloads.html) to convert coordinates between assemblies. <br />
**Genotype fields**.  By default, YAX  reads from the GT genotype field in VCF/BCF files (formatted as `0/0` or `0|0` for homozygous reference genotype).  Specify specify `--field DS` to instead use imputed genotype dosages (posterior mean alternate allele count) from the DS field (e.g., from [Minimac4](https://genome.sph.umich.edu/wiki/Minimac4)). <br />
**Missing genotypes**.  We highly recommend filtering and imputing genotype data prior to association analysis.   Genotype imputation software such as Beagle, IMPUTE, and Minimac can be used to accurately infer missing genotypes using a haplotype reference panel.  By default, YAX issues a warning when missing genotypes are detected, and sets missing genotypes to homozygous for the most common allele.    
## Molecular trait data
**Relevant flags:** `--bed {FILE}`. <br />
**File formats**. Molecular trait data must be stored in [BED file format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), and should be compressed using [BGZIP and indexed using Tabix](http://www.htslib.org/doc/tabix.html). For example, to BGZIP-compress and index a BED file called `ex.bed` one can run the command `bgzip ex.bed && tabix -p bed ex.bed.gz`. <br />
An example BED file with 2 individual samples and 2 genes is below:
```
#chr  start    end      gene_name        sample_1  sample_2
1     65418    65419    ENSG00000186092  -0.0837   -0.3476
1     827521   827522   ENSG00000225880   1.0369    1.3489
```
YAX requires the first 4 columns as shown above, where 1-3 (chr, start, end) specify the chromosomal  coordinates of each trait (for example, gene transcription start site [TSS] location) and 4 (gene_name) specifies the trait name or label (for example, Ensembl ID).   Any additional metadata columns (e.g., specifying strand orientation or other) will be ignored, and YAX will match the subsequent column names (`sample_1` and `sample_2` in the above example) to sample IDs present in other input files (e.g., genotype and covariate). <br />
**Missing data.** YAX does not currently support missing molecular trait values. We recommend filtering out traits with high proportions of missingness (e.g., >5%) prior to analysis, and imputing remaining missing values using [single-imputation](https://en.wikipedia.org/wiki/Imputation_%28statistics%29#Single_imputation) for other missing values.    
## Covariate data
**Relevant flags:** `--cov {FILE}`, `-c {FILE}`. <br />
**File format**. Covariate files are stored similar to molecular traits, with 1 row per covariate and 1 column per sample. YAX supports only numeric covariate data; any character-valued categorical variables must be converted to `0`-`1` dummy variables. For example, below is a covariate file with 4 samples and 3 covariates:
```
    #ID   sample_1    sample_2    sample_3    sample_4
    PC1   0.0139      0.0145      0.0141      0.0135
    PC2  -0.0097     -0.0059     -0.0025     -0.0064
    PC3   0.0067      0.0096     -0.0036     -0.0041
```
Covariate files should be white-space delimited.  Uncompressed or GZIP/BGZIP-compressed white-space delimited text files are supported.  Users are encouraged to verify that their covariate data matrix has full column rank.   

## Genetic relatedness and kinship matrices
**Relevant flags:** `--grm {FILE}`, `--kin {FILE}`. <br />
**File format**. Kinship or genetic relatedness matrices (GRMs) can be specified in a sparse matrix format as follows:
```
    #id1         id2            kinship
    sample_1     sample_1       0.50
    sample_1     sample_20      0.05
    sample_1     sample_22      0.15
```
If diagonal elements are listed in the kinship or GRM file (rows where `id1==id2`), then YAX analyzes the intersection of sample IDs present in the GRM/kinship,  genotype, trait, and covariate files.  Otherwise, YAX fixes diagonal elements of the kinship matrix to 0.5 (or 1 for GRMs), and assumes that any samples not listed in the GRM file (but present in other files) are unrelated.  Uncompressed or GZIP/BGZIP-compressed white-space delimited text files are currently supported.

