
# <img src="/doc/logo.svg" width="240" height="80"/>

APEX is a toolkit for analysis of molecular quantitative trait loci (xQTLs), such as mRNA expression or methylation levels. It performs cis and trans analysis, single- or multiple-variant analysis, and single-study or meta-analysis. 

To install APEX from source or download precompiled binaries, [**see installation guide below**](#installation).  

## Analysis modes and documentation
- [**xQTL factor analysis.**](/apex/doc/mode_factor/)  `./apex factor {options}`
	 - Rapid high-dimensional factor analysis to account for technical and biological variation in molecular trait data. 
- [**xQTL linear mixed model preprocessing.**](/apex/doc/mode_lmm/)  `./apex lmm {options}`
	 - Precompute null model terms for rapid LMM association analysis. 
- [**cis-xQTL analysis.**](/apex/doc/mode_cis/) `./apex cis {options}`
	 - Analyze xQTL associations within a genomic window of each molecular trait.  For example, analyze variant associations within 1 megabase (Mb) of each gene's transcription start site (TSS).  
 - [**trans-xQTL analysis.**](/apex/doc/mode_trans/)  `./apex trans {options}`
	 - Analyze xQTL genome-wide associations for each molecular trait, optionally excluding variants within the cis region. 
 - [**xQTL meta-analysis.**](/apex/doc/mode_meta/)  `./apex meta {options}`
	 -  Meta-analyze xQTL association summary statistics across one or multiple studies.   
 - [**xQTL summary data storage.**](/apex/doc/mode_store/)  `./apex store {options}`
	 -  Store genotypic variance-covariance (LD) data matrices for data-sharing or meta-analysis. These files enable accurate multiple-variant analysis from xQTL summary statistics, including joint and conditional regression analysis, Bayesian finemapping, aggregation tests, and penalized regression / variable selection.    
 
## Additional documentation pages
 See above for documentation pages on each analysis mode (cis, trans, meta, store).  Additional documentation pages are listed below. 
 - **[Input file formats for individual-level analysis](/apex/doc/input_files)** in modes `apex cis`, `apex trans`, and `apex store`. 
 - **APEX xQTL summary association data file formats** generated by `apex store` and used in `apex meta`. *Documentation coming soon.*
 - **Bayesian finemapping using APEX xQTL summary association data files** generated by `apex store`. *Documentation coming soon.*
 
## Release notes
 - apex-factor
	 - Currently, apex-factor assumes sample size is less than the number of molecular traits.  In the next release, we will provide optimizations for sample size > number of traits. 
	 - Currently, apex-factor calculates inferred covariates in a format suitable for OLS or GRM LMMs.  To model factor covariates as random effects, inferred covariates must be calculated internally in apex-cis and apex-trans.  In the next release, apex-factor will provide inferred factor covariates in a format suitable for low-rank LMMs in apex-cis and apex-trans. 
 - apex-lmm
	 - Currently, apex-lmm supports LMMs with GRMs, but not low-rank LMMs. Low-rank LMMs can be estimated directly in apex-cis and apex-trans.  In the next release, apex-lmm will also support low-rank LMMs. 
 - apex-cis and apex-trans
	 - Currently, pre-computed variance component estimates from apex-lmm can be used in apex-cis and apex-trans, but LMM trait residuals and genotypic variance interpolation points are precalculated internally in apex-cis and apex-trans.  In the next release, apex-cis and apex-trans will accept precomputed LMM trait residuals and genotypic variance interpolation points from apex-lmm.  
 - apex-store
	 - Currently, apex-store provides adjusted LD for OLS models.  In the next release, apex-store will provide options for adjusted LD in LMMs. 
 - apex-meta
	 - Currently, apex-meta supports only cis-xQTL meta-analysis. In the next release, apex-meta will also support trans-xQTL meta-analysis. 
 
## Installation
APEX is primarily written in C++. To compile APEX from source, run:
```
git clone https://github.com/corbinq/apex.git
cd apex 
bash get_dependencies.sh
make
```
Precompiled binaries are also available for 64-bit Linux systems as follows:
```
git clone https://github.com/corbinq/apex.git
cd apex/bin
gunzip apex_Linux_x86_64.gz
mv apex_Linux_x86_64 apex && chmod +x apex
./apex --help
```

## Software dependencies

 - [Eigen C++ matrix library](http://eigen.tuxfamily.org/)
 - [HTSlib C library](http://www.htslib.org/)
 - [github.com/jonathonl/shrinkwrap](https://github.com/jonathonl/shrinkwrap)
 - [github.com/Taywee/args](https://github.com/Taywee/args)
 - [Boost math C++ library](https://www.boost.org/)

## Contact, feedback, bug reports
APEX is a work in progress, and under active development. We appreciate feedback and bug reports, and welcome any questions or feature requests. **Contact:** <qcorbin@hsph.harvard.edu>.  

## Citation
If you use APEX, please cite `Quick, C; Guan, L; Lin, X (2020). URL: https://github.com/corbinq/apex` (preprint pending). Thanks!

## Acknowledgements
APEX is developed by Corbin Quick and Li Guan. Special thanks to 

 - Yaowu Liu
 - Xihao Li
 - Zilin Li
 - Rounak Dey
 - Laura Scott
 - Hyun Min Kang 
 - Xihong Lin 



