cmake_minimum_required(VERSION 3.14)
project(apex VERSION 0.0.1)
include(CTest)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 14)

if(NOT APPLE)
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  set(BUILD_SHARED_LIBS OFF)
  set(CMAKE_EXE_LINKER_FLAGS "-static")
endif()

if(CGET_PREFIX)
  include_directories(${CGET_PREFIX}/include)
  include_directories(${CGET_PREFIX}/include/eigen3)
  include_directories(${CGET_PREFIX}/include/Rmath)
endif()

find_library(LIB_ZSTD zstd HINT ${CGET_PREFIX}/lib)
find_library(LIB_LZMA lzma HINT ${CGET_PREFIX}/lib)
find_library(LIB_Z z HINT ${CGET_PREFIX}/lib)
find_library(LIB_HTS hts HINT ${CGET_PREFIX}/lib)
find_library(LIB_RMATH Rmath HINT ${CGET_PREFIX}/lib)

find_package(BZip2)
find_package(OpenMP)

#set(BOOST_INCLUDEDIR ${CGET_PREFIX}/include)
#set(BOOST_LIBRARYDIR ${CGET_PREFIX}/lib)
#find_package(Boost REQUIRED)

#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

set(APEX_HEADERS
    src/apexLMM.hpp
    src/cisMapping.hpp
    src/dataParser.hpp
    src/eigenData.hpp
    src/fitUtils.hpp
    src/fitModels.hpp
    src/htsWrappers.hpp
    src/genotypeData.hpp
    src/rmathWrappers.hpp
    src/mathStats.hpp
    src/metaAnalysis.hpp
    src/miscUtils.hpp
    src/processVCOV.hpp
    src/mapID.hpp
    src/readBED.hpp
    src/readTable.hpp
    src/setOptions.hpp
    src/setRegions.hpp
    src/transMapping.hpp
    src/xzFormat.hpp
    src/BRENT/brent_fmin.hpp
)

set(APEX_SOURCE
  src/apexLMM.cpp
  src/cisMapping.cpp
  src/eigenData.cpp
  src/fitModels.cpp
  src/fitUtils.cpp
  src/genotypeData.cpp
  src/htsWrappers.cpp
  src/rmathWrappers.cpp
  src/mapID.cpp
  src/mathStats.cpp
  src/metaAnalysis.cpp
  src/miscUtils.cpp
  src/processVCOV.cpp
  src/readBED.cpp
  src/readTable.cpp
  src/scanSignals.cpp
  src/setOptions.cpp
  src/setRegions.cpp
  src/transMapping.cpp
  src/BRENT/brent_fmin.cpp
)

add_library(apex ${APEX_HEADERS} ${APEX_SOURCE})

# Do not use LTO in debug mode, symbols get messed up & breakpoints will not work correctly
target_compile_options(apex PUBLIC "$<$<NOT:$<PLATFORM_ID:Darwin>>:-static>" "$<$<CONFIG:RELEASE>:-flto>" -Wall -Wpedantic -fno-math-errno)

# -static not supported on MacOS, can't statically link to C runtime
target_link_options(apex PUBLIC "$<$<NOT:$<PLATFORM_ID:Darwin>>:-static>" "$<$<CONFIG:RELEASE>:-flto>" -Wno-unused-function)
target_include_directories(apex PUBLIC ${CGET_PREFIX}/include src src/BRENT ${BZIP2_INCLUDE_DIRS})
target_link_libraries(apex PUBLIC ${LIB_ZSTD} ${LIB_LZMA} ${LIB_Z} ${LIB_HTS} ${LIB_RMATH} ${Boost_LIBRARIES} ${BZIP2_LIBRARIES})

if(OpenMP_CXX_FOUND)
  target_compile_options(apex PUBLIC "${OpenMP_CXX_FLAGS}")
  target_link_libraries(apex PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(apex-bin
  src/Main.cpp
  src/Main.hpp
)
target_link_libraries(apex-bin PUBLIC apex)
set_target_properties(apex-bin PROPERTIES OUTPUT_NAME apex)

if(BUILD_TESTING)
  add_subdirectory(tests)
endif()