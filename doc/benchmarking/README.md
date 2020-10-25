
# YAX: cis-xQTL analysis guide
This page describes benchmarking experiments with YAX for various QTL analysis tasks.  This is intended to serve as a guide for time and memory usage requirements,  identify and delineate numerical  differences between software tools, and highlight unique features of YAX.  <br />

##### Table of Contents  

 1. [LMM Benchmarking](#lmm-benchmarking)   
 2. [Trait-Level Test Benchmarking](#trait-level-test-benchmarking)
 3. [Command line options](#command-line-arguments) <br />

 [*Return to YAX main page.*](/../../)


## LMM Benchmarking

We used empirical genotype data derived from the UK Biobank and simulated molecular phenotype data to benchmark time and memory for trans-xQTL analysis using a linear mixed model (LMM) using YAX LMM, FastGWA, BOLT-LMM, and GMMAT. 

### Input data
| Sample size |No. traits | No. SNPs | No. covariates |
|-------------|-----------------:|---------------:|---------------------:|
| 10,000      |           16,329 |        590,606 |                   10 |

### Time and memory usage
|                 |     CPU   hours    |     Time   speedup    |     Max   memory    |
|-----------------|-------------------:|----------------------:|--------------------:|
|     YAX, p < 5e-5    |             7.5    |               0.36    |        4.88   Gb    |
|     YAX         |            20.8    |              1.00*    |       4.88   Gb     |
|     FastGWA     |            52.1    |               2.50    |        0.14   Gb    |
|     BOLT-LMM    |         1,068.9    |              51.39    |        0.67   Gb    |
|     GMMAT       |       ~*5,692.5*     |            ~*273.68*    |             N/A     |


### Software commands
**YAX command:**
 - parallel over chromosomes, all traits at once
 - sparse GRM 
 - compressed BCF format genotype data
```
yax trans --bed $all_traits_bed --bcf $bcf --cov $covar_txt --grm $grm --region ${chr} --out trans_chr${chr}
```
**BOLT-LMM command:**
 - parallel over traits, all chroms at once 
 - 172,045 LD-pruned SNPs (no sparse GRM)
 - uncompressed BED/BIM/FAM format genotype data
```
bolt --lmm --LDscoresFile=$ldsc_f --bfile=$bfile --phenoFile=trait_${trait}.ped --phenoCol=${trait} --qCovarCol=PC{1:10} --covarFile=$covar_ped --modelSnps=${snp_file} --statsFile=bolt_${trait}
```
**fastGWA command:**
 - parallel over traits, all chroms at once 
 - Sparse GRM
 - uncompressed BED/BIM/FAM format genotype data
```
gcta64 --fastGWA-mlm --bfile $bfile --grm-sparse $sgrm --pheno trait_${trait}.ped --qcovar $covar_ped --threads 1 --out gcta_${trait}
```
**GMMAT command:**
 - parallel over traits, parallel over chroms
 - Sparse GRM
 - Compressed GDS format genotype data
```
eqn <- as.formula( Y ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10)
null_fit <- GMMAT::glmmkin(
	eqn, data = trait_data, id = "id", 
	kins = GRM, family = gaussian(), 
	method = "REML", method.optim = "AI", verbose = TRUE
)
GMMAT::glmm.score(obj = null_fit, infile = gds_file, outfile = trait_out_file, verbose = TRUE)
```

**Software concordance:**
LMM association tests from YAX and GMMAT are nearly numerically equivalent, as expected.  BOLT-LMM uses the pre-conditioned conjugate gradient method to avoid storing an explicit GRM, and a retrospective quasi-likelihood score test; these differences may explain differences with YAX and GMMAT.  FastGWA uses the GRAMMAR-Gamma approximation to calculate association tests, which may  explain  differences with YAX and GMMAT.  

## Trait-Level Test Benchmarking



