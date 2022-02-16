enable_testing()

find_library(LIB_GTEST gtest HINT ${CGET_PREFIX}/lib)
find_library(LIB_GTEST_MAIN gtest_main HINT ${CGET_PREFIX}/lib)

get_target_property(APEX_INCLUDES apex INCLUDE_DIRECTORIES)

add_executable(
    tests
    reader_util.h
    meta_stepwise_reader.h
    meta_stepwise_reader.cpp
    test_meta.cpp
    test_rmath.cpp
)

target_include_directories(tests PUBLIC ${APEX_INCLUDES})

target_link_libraries(
    tests
    apex
    ${LIB_GTEST}
    ${LIB_GTEST_MAIN}
)