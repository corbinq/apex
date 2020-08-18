#!/bin/bash

# args library for processing command line arguments
git clone https://github.com/Taywee/args.git

# Eigen library for linear algebra
git clone https://gitlab.com/libeigen/eigen.git

# Spectra library for sparse, high-dimensional eigenvalue problems
git clone https://github.com/yixuan/spectra.git

# Brent's algorithm

# mkdir -p BRENT && cd BRENT
# wget --no-check-certificate https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.hpp
# wget --no-check-certificate https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.cpp
# cd -


# boost::math library 
module load boost
if [ $(whereis boost | cut -d: -f2 | wc -c) -lt 0 ]; 
then
	git clone https://github.com/boostorg/boost
	echo "WARNING: Could not find boost library on system.\n"
	echo "	 Preparing to install boost locally from source.\n"
	echo "NOTE: Installing boost will take several minutes.\n"
	sleep 2
	cd boost
	git submodule update --init
	bash bootstrap.sh
	./b2 headers
	cd ../
else
	echo "Found boost installed at $(whereis boost | cut -d: -f2)"
fi

# htslib library, which we use for BCF/VCF access and indexing
module load htslib
tbx_loc=$(echo '#include <htslib/tbx.h>' | cpp -H -o /dev/null 2>&1 | head -n1)

if [[ $tbx_loc == *"error"* ]]; then 
	echo "Could not find htslib on system."
	echo "Installing htslib locally."

	git clone https://github.com/samtools/htslib.git
	cd htslib
	autoreconf
	./configure --disable-lzma --disable-bz2 --disable-libcurl
	make install prefix=$PWD
	cd ../
else
	echo "Found htslib installed on system."
fi

# GDS format :
# git clone https://github.com/CoreArray/GDSFormat.git

# module load boost
:
