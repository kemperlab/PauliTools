
#include "cartan/PauliVector.h"
#include "cartan/symplectic.h"

PauliVector::PauliVector() {};

PauliVector::PauliVector(const PauliSingleton& p)
{
    data[p.pauli.code] = p.coef;
}

PauliVector::PauliVector(const cdouble a, const Pauli& p)
{
    data[p.code] = a;
}

PauliVector::PauliVector(const Pauli& p, const cdouble a)
{
    data[p.code] = a;
}

PauliVector::PauliVector(const double a, const Pauli& p)
{
    data[p.code] = a;
}

PauliVector::PauliVector(const Pauli& p, const double a)
{
    data[p.code] = a;
}

void PauliVector::add_or_insert(const Pauli& newPauli, const cdouble coef)
{
    if(auto search = data.find(newPauli.code); search != data.end())
        data[newPauli.code] += coef;
    else
        data[newPauli.code] = coef;
}

PauliVector& PauliVector::operator+= (const PauliVector& pv2)
{
    for (auto const& [key, val] : pv2.data)
    {
        this->add_or_insert(key, val);
    }
    return *this;
}

PauliVector operator+(PauliVector p1, const PauliVector& p2)
{
    p1 += p2;
    return p1;
}


std::ostream& operator<< (std::ostream& os, const PauliVector& pv)
{
    for (auto const& [key, val] : pv.data)
    {
        os << key << " " << val << std::endl;
    }
    return os;
}

PauliVector operator*( cdouble a, Pauli p)
{
    return PauliVector( p, a);
}

PauliVector operator*( Pauli p, cdouble a)
{
    return PauliVector( p, a);
}

PauliVector operator*( double a, Pauli p)
{
    return PauliVector( p, a);
}

PauliVector operator*( Pauli p, double a)
{
    return PauliVector( p, a);
}


PauliVector& PauliVector::squeeze()
{
    const auto count = std::erase_if(data, 
            [](const auto& item) {   // Lambda function
                auto const& [pauli, coef] = item;
                return abs(coef) < PAULI_VECTOR_IS_ZERO_TOL;
            });
    return *this;
}

// Commutators
PauliVector PauliVector::ad(const PauliVector& pv2)
{
    return commute(*this, pv2);
}

PauliVector commute(const PauliVector& pv1, const PauliVector& pv2)
{
    PauliVector p3;
    for (auto const& [Pauli1, coef1] : pv1.data)
    for (auto const& [Pauli2, coef2] : pv2.data)
    {
        auto comm = symplectic_binary_commute_with_prefactor(
                Pauli1.code, Pauli2.code);
        /*
        std::cout << Pauli1 << " " << coef1 << std::endl;
        std::cout << Pauli2 << " " << coef2 << std::endl;
        std::cout << Pauli(comm.first) << " " << comm.second * coef1 * coef2 << std::endl;
        std::cout << "\n" << std::endl;
        */
        p3.add_or_insert(comm.first, comm.second * coef1 * coef2);
    }
    //std::cout << p3 << std::endl;
    //std::cout << &p3 << std::endl;
    return p3.squeeze();
}
