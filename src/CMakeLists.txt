add_library( pauli Pauli.cc involutions.cc symplectic.cc )
target_include_directories( pauli PUBLIC ${PROJECT_SOURCE_DIR}/include )

add_library( dla dla_generator.cc )
target_include_directories( dla PUBLIC ${PROJECT_SOURCE_DIR}/include )

add_library( models models.cc )
target_include_directories( models PUBLIC ${PROJECT_SOURCE_DIR}/include )

add_library( paulivector PauliVector.cc )
target_include_directories( paulivector PUBLIC ${PROJECT_SOURCE_DIR}/include )

# openmp
if( TARGET OpenMP::OpenMP_CXX )
  target_link_libraries( pauli PUBLIC OpenMP::OpenMP_CXX )
  target_link_libraries( dla PUBLIC OpenMP::OpenMP_CXX )
endif()
