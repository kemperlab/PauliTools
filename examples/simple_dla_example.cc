// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022

#include "cartan/Pauli.h"
#include "cartan/algebra_types.h"
#include "cartan/dla_generator.h"
#include "cartan/models.h"
#include "cartan/involutions.h"

int main()
{
    setNq(6);
    std::vector<std::string> ham = heisenberg(false);

    // Fill the Hamiltonian algebra objects
    algebra algebra_g, algebra_h;

    for(const auto& hamterm : ham)
    {
        algebra_g.insert(Pauli(hamterm));
    }

    // Produce an initial guess for algebra_h
    for(size_t i=0; i < ham.size(); i += 6)
    {
        algebra_h.insert(Pauli(ham[i]));
        algebra_h.insert(Pauli(ham[i+1]));
    }


    auto&& [algebra_k_1, algebra_m_1 ] = get_algebra_by_commuting(
            &algebra_g, &algebra_h, parityY);


    algebra_g.clear();
    for(const auto& hamterm : ham)
    {
        algebra_g.insert(Pauli(hamterm));
    }

    get_algebra_by_commuting(&algebra_g);

    algebra_g.clear();
    for(const auto& hamterm : ham)
    {
        algebra_g.insert(Pauli(hamterm));
    }
    get_algebra_by_commuting(&algebra_g, nullptr, parityY);

    return 0;

}
