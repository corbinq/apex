/*  
    Copyright (C) 2020 
    Author: Corbin Quick <qcorbin@hsph.harvard.edu>

    This file is a part of APEX.

    APEX is distributed "AS IS" in the hope that it will be 
    useful, but WITHOUT ANY WARRANTY; without even the implied 
    warranty of MERCHANTABILITY, NON-INFRINGEMENT, or FITNESS 
    FOR A PARTICULAR PURPOSE.

    The above copyright notice and disclaimer of warranty must 
    be included in all copies or substantial portions of APEX.
*/

#ifndef APEX_CLI_HPP
#define APEX_CLI_HPP

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <cstdlib>
#include <csignal>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "args.hxx"
#include "setOptions.hpp"
#include "htsWrappers.hpp"
#include "setRegions.hpp"
#include "readBED.hpp"
#include "readTable.hpp"
#include "genotypeData.hpp"
#include "transMapping.hpp"
#include "cisMapping.hpp"
#include "metaAnalysis.hpp"
#include "fitUtils.hpp"
#include "mapID.hpp"
#include "miscUtils.hpp"
#include "mathStats.hpp"
#include "processVCOV.hpp"

// ------------------------------------
//  Analysis modes
// ------------------------------------

int cis(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int trans(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int factor(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int lmm(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int meta(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);
int store(const std::string &progname, std::vector<std::string>::const_iterator beginargs, std::vector<std::string>::const_iterator endargs);

using mode_fun = std::function<int(const std::string &, std::vector<std::string>::const_iterator, std::vector<std::string>::const_iterator)>;

#endif

