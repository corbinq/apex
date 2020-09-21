# YAX: *Y*et *a*nother *x*QTL toolkit
YAX is a comprehensive software toolkit for analysis of molecular quantitative trait loci (xQTLs) including eQTLs (mRNA expression) and  mQTL (methylation). Some unique features of YAX include

  **xQTL mapping**
 - Highly optimized linear mixed model (LMM) framework to account for cryptic or familial relatedness in cis and trans xQTL analysis. 
 - Variable selection and conditional analysis procedures to identify multiple association signals for a single molecular trait.
 
**xQTL meta-analysis**
 - Single and multiple-variant xQTL meta-analysis procedures. 
 - Highly condensed storage formats for xQTL summary statistics, including study-specific LD information.  These allow joint and conditional analysis, Bayesian finemapping, and aggregation tests across 1 or more xQTL studies while protecting genetic privacy.  
## Installation
YAX is primarily written in C++, and depends on the Eigen matrix library, HTSlib, and Boost. To compile YAX from source, run:

    git clone https://github.com/corbinq/yax.git
    cd yax 
    bash get_dependencies.sh
    make
Precompiled binaries are also available as follows:
## Analysis modes
**cis-xQTL analysis:** `./yax cis {options}`
**trans-xQTL analysis**  `./yax trans {options}`
**xQTL meta-analysis**  `./yax meta {options}`
**xQTL summary data storage**  `./yax store {options}`
## Software dependencies

 - [Eigen C++ matrix library](http://eigen.tuxfamily.org/)
 - [HTSlib C library](http://www.htslib.org/)
 - [github.com/Taywee/args](https://github.com/Taywee/args)
 - [Boost math C++ library](https://www.boost.org/)

## Acknowledgements
 - Li Guan
 - Laura Scott
 - Xihao Li
 - Zilin Li
 - Rounak Dey
 - Xihong Lin


