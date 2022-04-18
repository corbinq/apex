enable_testing()

find_library(LIB_GTEST gtest HINT ${CGET_PREFIX}/lib)
find_library(LIB_GTEST_MAIN gtest_main HINT ${CGET_PREFIX}/lib)
find_library(LIB_BOOST_IOSTREAMS boost_iostreams HINT ${CGET_PREFIX}/lib)

get_target_property(APEX_INCLUDES apex INCLUDE_DIRECTORIES)

add_executable(
    tests
    main.cpp
    reader_util.h
    meta_svar_reader.h
    meta_svar_reader.cpp
    meta_stepwise_reader.h
    meta_stepwise_reader.cpp
    cis_gene_table_reader.h
    cis_gene_table_reader.cpp
    cis_long_table_reader.h
    cis_long_table_reader.cpp
    test_meta.cpp
    test_rmath.cpp
    test_cis.cpp
    test_store.cpp
)

target_include_directories(tests PUBLIC ${APEX_INCLUDES})

target_link_libraries(
    tests
    apex
    ${LIB_GTEST}
    ${LIB_BOOST_IOSTREAMS}
)

add_dependencies(tests apex-bin)