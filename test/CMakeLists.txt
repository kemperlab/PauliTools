add_executable ( pauli_tests Pauli.cc )
target_link_libraries( pauli_tests PUBLIC pauli GTest::gtest_main )
#target_include_directories( pauli_tests PUBLIC  )

add_executable( dla_tests DLA.cc )
target_link_libraries( dla_tests PUBLIC pauli dla models GTest::gtest_main )
#target_include_directories

#${PROJECT_SOURCE_DIR}/test 
#${PROJECT_SOURCE_DIR}/src )


enable_testing()


include(GoogleTest)
gtest_discover_tests(pauli_tests)
gtest_discover_tests(dla_tests)
