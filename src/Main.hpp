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


#ifndef APEX_HPP
#define APEX_HPP

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

#include "args/args.hxx"

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


#endif

