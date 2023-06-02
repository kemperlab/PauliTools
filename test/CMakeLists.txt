
add_executable( run_all_tests 
    run.cc 
    Pauli.cc
    DLA.cc
    )

target_link_libraries( run_all_tests PUBLIC 
    dla
    pauli 
    models
    GTest::gtest_main )


add_executable ( pauli_tests Pauli.cc )
target_link_libraries( pauli_tests PUBLIC pauli GTest::gtest_main )

add_executable( dla_tests DLA.cc )
target_link_libraries( dla_tests PUBLIC dla pauli models GTest::gtest_main )

add_executable ( pauli_vector_tests PauliVector.cc )
target_link_libraries( pauli_vector_tests PUBLIC paulivector pauli GTest::gtest_main )

enable_testing()


include(GoogleTest)
gtest_discover_tests(pauli_tests)
gtest_discover_tests(dla_tests)
