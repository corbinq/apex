
# YAX: cis-xQTL analysis guide
This page describes benchmarking experiments with YAX for various QTL analysis tasks.  This is intended to serve as a guide for time and memory usage requirements,  identify and delineate numerical  differences between software tools, and highlight unique features of YAX.  <br />

##### Table of Contents  

 1. [LMM Benchmarking](#lmm-benchmarking)   
 2. [Trait-Level Test Benchmarking](#trait-level-test-benchmarking)
 3. [Command line options](#command-line-arguments) <br />

 [*Return to YAX main page.*](https://github.com/corbinq/yax)


## LMM Benchmarking

We used empirical genotype data derived from the UK Biobank and simulated molecular phenotype data to benchmark time and memory for trans-xQTL analysis using a linear mixed model (LMM) using YAX LMM, FastGWA, BOLT-LMM, and GMMAT. 

**Input data:**
| Sample size |No. traits | No. SNPs | No. covariates |
|-------------|-----------------:|---------------:|---------------------:|
| 10.000      |           16,329 |        590,606 |                   10 |

**Time and memory usage:**
|                 |     CPU   hours    |     Time   speedup    |     Max   memory    |
|-----------------|-------------------:|----------------------:|--------------------:|
|     YAX, p < 5e-5    |             7.5    |               0.36    |        4.88   Gb    |
|     YAX         |            20.8    |              1.00*    |       4.88   Gb     |
|     FastGWA     |            52.1    |               2.50    |        0.14   Gb    |
|     BOLT-LMM    |         1,068.9    |              51.39    |        0.67   Gb    |
|     GMMAT       |       ~*5,692.5*     |            ~*273.68*    |             N/A     |

**Software concordance:**
LMM association tests from YAX and GMMAT are nearly numerically equivalent, as expected.  BOLT-LMM uses the pre-conditioned conjugate gradient method to avoid storing an explicit GRM, and a retrospective quasi-likelihood score test; these differences may explain differences with YAX and GMMAT.  FastGWA uses the GRAMMAR-Gamma approximation to calculate association tests, which may  explain  differences with YAX and GMMAT.  
## Trait-Level Test Benchmarking



sbatch --wrap="/usr/bin/time -v --output=chr20_yax.log ~/.local/bin/gqt store --bcf chr20.bcf --out test_chr20" -p xlin,shared --mem 8000 -t 2-00:00 -N 1 -n 1 -J chr20



