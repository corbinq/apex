cmake_minimum_required(VERSION 3.14)
project(apex VERSION 0.0.1)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)

if(CGET_PREFIX)
  include_directories(${CGET_PREFIX}/include)
  include_directories(${CGET_PREFIX}/include/eigen3)
endif()

find_library(LIB_ZSTD zstd HINT ${CGET_PREFIX}/lib)
find_library(LIB_LZMA lzma HINT ${CGET_PREFIX}/lib)
find_library(LIB_Z z HINT ${CGET_PREFIX}/lib)
find_library(LIB_HTS hts HINT ${CGET_PREFIX}/lib)

find_package(BZip2)
find_package(OpenMP)

#set(BOOST_INCLUDEDIR ${CGET_PREFIX}/include)
#set(BOOST_LIBRARYDIR ${CGET_PREFIX}/lib)
#find_package(Boost REQUIRED)

#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(SOURCE_FILES 
  src/Main.cpp
  src/Main.hpp
  src/apexLMM.cpp
  src/apexLMM.hpp
  src/cisMapping.cpp
  src/cisMapping.hpp
  src/dataParser.hpp
  src/eigenData.cpp
  src/eigenData.hpp
  src/fitModels.cpp
  src/fitModels.hpp
  src/fitUtils.cpp
  src/fitUtils.hpp
  src/genotypeData.cpp
  src/genotypeData.hpp
  src/htsWrappers.cpp
  src/htsWrappers.hpp
  src/mapID.cpp
  src/mapID.hpp
  src/mathStats.cpp
  src/mathStats.hpp
  src/metaAnalysis.cpp
  src/metaAnalysis.hpp
  src/miscUtils.cpp
  src/miscUtils.hpp
  src/processVCOV.cpp
  src/processVCOV.hpp
  src/readBED.cpp
  src/readBED.hpp
  src/readTable.cpp
  src/readTable.hpp
  src/scanSignals.cpp
  src/setOptions.cpp
  src/setOptions.hpp
  src/setRegions.cpp
  src/setRegions.hpp
  src/transMapping.cpp
  src/transMapping.hpp
  src/xzFormat.hpp
  src/BRENT/brent_fmin.cpp
  src/BRENT/brent_fmin.hpp)
add_executable(apex ${SOURCE_FILES})

# Do not use LTO in debug mode, symbols get messed up & breakpoints will not work correctly
target_compile_options(apex PUBLIC "$<$<NOT:$<PLATFORM_ID:Darwin>>:-static>" "$<$<CONFIG:RELEASE>:-flto>" -Wall -Wpedantic -fno-math-errno)

# -static not supported on MacOS, can't statically link to C runtime
target_link_options(apex PUBLIC "$<$<NOT:$<PLATFORM_ID:Darwin>>:-static>" "$<$<CONFIG:RELEASE>:-flto>" -Wno-unused-function)
target_include_directories(apex PUBLIC ${CGET_PREFIX}/include src/BRENT ${BZIP2_INCLUDE_DIRS})
target_link_libraries(apex PUBLIC ${LIB_ZSTD} ${LIB_LZMA} ${LIB_Z} ${LIB_HTS} ${Boost_LIBRARIES} ${BZIP2_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_compile_options(apex PUBLIC "${OpenMP_CXX_FLAGS}")
  target_link_libraries(apex PUBLIC OpenMP::OpenMP_CXX)
endif()
