# PauliTools

PauliTools is a set of libraries useful for working with Pauli strings. 
The code is written in C++ and takes full advantage of the symplectic binary representation
of Pauli strings.  Current capabilities include generating the Hamiltonian algebra (AKA
Dynamical Lie Algebra) of a set of Pauli strings.

This the algorithm can be found in:

Fixed Depth Hamiltonian Simulation via Cartan Decomposition
Efekan Kökcü, Thomas Steckmann, J. K. Freericks, Eugene F. Dumitrescu, Alexander F. Kemper
Phys. Rev. Lett. 129, 070501 (2022)

https://kemperlab.github.io/papers/Kokcu-Cartan/



## How to run?

The PauliTools package uses the CMake build system (CMake version ≥ 3.16).
The recommended way of building PauliTools is as follows:

1. Install

        git clone https://github.com/kemperlab/PauliTools.git

2. CMake

        cd PauliTools
        mkdir build
        cd build
        cmake ..
        make -j8

3. Run tests

        ./test/dla_tests
        ./test/pauli_tests
        
        
## Developers - NC State University
- [Alexander Kemper](https://kemperlab.github.io/) - akemper@ncsu.edu
