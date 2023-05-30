cmake_minimum_required( VERSION 3.26.3 FATAL_ERROR )

# project
project( CARTAN VERSION 0.1.0 LANGUAGES CXX )

# openmp
find_package( OpenMP )


execute_process(
    COMMAND ping www.github.com -c 2
    ERROR_QUIET
    RESULT_VARIABLE NO_CONNECTION
    )

#if(NOT NO_CONNECTION EQUAL 0)
#set(FETCHCONTEXT_FULLY_DISCONNECTED ON)
#message(WARNING "Fetch offline mode: requires already populated _deps")
#else()
#set(FETCHCONTENT_FULLY_DISCONNECTED OFF)
#endif()

# fetch content

include( FetchContent )

# GTest
FetchContent_Declare( gtest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG main
)

FetchContent_MakeAvailable( gtest )



# src
add_subdirectory( src )

# test
add_subdirectory( test )

# examples
add_subdirectory( examples )

# iterative cartan
# add_subdirectory( iterative_cartan )



SET(GCC_COVERAGE_COMPILE_FLAGS "-O3 -fopenmp")
#SET(GCC_COVERAGE_LINK_FLAGS    "-lgcov")
SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

