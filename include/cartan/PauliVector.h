
#pragma once

#include "cartan/Pauli.h"
#include <map>
#include <complex>
#include <utility>
#include <iostream>

typedef std::complex<double> cdouble;
inline cdouble II(0,1);

const double PAULI_VECTOR_IS_ZERO_TOL = 1e-10;

/** Struct
 * \struct PauliSingleton
 * \brief Special case: a PauliVector with a single string
 */
struct PauliSingleton
{
    Pauli pauli;
    cdouble coef;
};

  /**
   * \class PauliVector
   * \brief Class to handle elements of the Pauli algebra/vector space
   */
class PauliVector
{
    public:

    // Constructors
    PauliVector();
    PauliVector(const double, const Pauli&);
    PauliVector(const Pauli&, const double);
    PauliVector(const cdouble, const Pauli&);
    PauliVector(const Pauli&, const cdouble);
    PauliVector(const PauliSingleton&);

    // Destructor
    ~PauliVector() {}

    // Addition
    PauliVector& operator+= (const PauliVector&);

    // Commutator
    PauliVector ad(const PauliVector&);

    // Take out null elements
    PauliVector& squeeze();

    // Utility function for adding to the data structure
    void add_or_insert(const Pauli& newPauli, const cdouble coef);

    std::map<Pauli, cdouble> data;
};




std::ostream& operator<< (std::ostream& os, const PauliVector& pv);

PauliVector operator+(PauliVector p1, const PauliVector& p2);

// Syntactic sugar for making PauliVectors
PauliVector operator*(double a, Pauli p);
PauliVector operator*( Pauli p, double a);
PauliVector operator*( cdouble a, Pauli p);
PauliVector operator*( Pauli p, cdouble a);

// Commutators etc
PauliVector commute(const PauliVector&, const PauliVector&);
