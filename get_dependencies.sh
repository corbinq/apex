#!/bin/bash
set -euxo pipefail
IFS=$'\n\t'

# args library for processing command line arguments
git clone https://github.com/Taywee/args.git src/args

# Eigen library for linear algebra
git clone https://gitlab.com/libeigen/eigen.git src/eigen

# Spectra library for sparse, high-dimensional eigenvalue problems
#git clone https://github.com/yixuan/spectra.git src/spectra
git clone --single-branch --branch 0.9.x https://github.com/yixuan/spectra.git src/spectra

git clone https://github.com/jonathonl/shrinkwrap.git src/shrinkwrap

# Brent's algorithm

# mkdir -p BRENT && cd BRENT
# wget --no-check-certificate https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.hpp
# wget --no-check-certificate https://people.sc.fsu.edu/~jburkardt/cpp_src/brent/brent.cpp
# cd -


# try loading modules  
if [ `which module` ]; then
	module load boost
#	module load htslib
fi

# boost::math library
if [ ! `which boost` ];
then
	git clone https://github.com/boostorg/boost src/boost
	echo "WARNING: Could not find boost library on system.\n"
	echo "	 Preparing to install boost locally from source.\n"
	echo "NOTE: Installing boost will take several minutes.\n"
	sleep 2
	cd src/boost
	git submodule update --init
	bash bootstrap.sh
	./b2 headers
	cd ../../
else
	echo "Found boost installed at $(which boost)"
fi

# htslib library, which we use for BCF/VCF access and indexing
tbx_loc=$(echo '#include <htslib/tbx.h>' | cpp -H -o /dev/null 2>&1 | head -n1)

# if [[ $tbx_loc == *"error"* ]]; then 
# 	echo "Could not find htslib on system."
	echo "Installing htslib locally."

	git clone https://github.com/samtools/htslib.git src/htslib
	cd src/htslib
	git submodule update --init --recursive
	autoreconf
	./configure --disable-lzma --disable-bz2 --disable-libcurl
	make install prefix=$PWD
	cd ../../
# else
#	echo "Found htslib installed on system."
# fi

# GDS format :
# git clone https://github.com/CoreArray/GDSFormat.git

# module load boost


