/// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022

#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>
#include <complex>
#include <sstream>
    using namespace std;

#include "cartan/Pauli.h"
#include "cartan/PauliVector.h"
#include "cartan/symplectic.h"

typedef complex<double> cdouble;

pair<cdouble, string> reference_pauli_product(string p1, string p2)
{
    char buffer[Nq+1];
    buffer[Nq] = '\0';
    cdouble prefactor = 1;
    for(int i=0; i < Nq; i++)
    {
        if(p1[i] == 'I' or p1[i] == '-')
        {
            buffer[i] = p2[i];
            continue;
        }

        if(p2[i] == 'I' or p2[i] == '-')
        {
            buffer[i] = p1[i];
            continue;
        }

        if(p1[i] == p2[i])
        {
            buffer[i] = 'I';
            continue;
        }

        if(p1[i] == 'X' and p2[i] == 'Y')
        {
            buffer[i] = 'Z';
            prefactor *= II;
        }
        else if(p1[i] == 'Y' and p2[i] == 'X')
        {
            buffer[i] = 'Z';
            prefactor *= -II;
        }
        else if(p1[i] == 'Y' and p2[i] == 'Z')
        {
            buffer[i] = 'X';
            prefactor *= II;
        }
        else if(p1[i] == 'Z' and p2[i] == 'Y')
        {
            buffer[i] = 'X';
            prefactor *= -II;
        }
        else if(p1[i] == 'Z' and p2[i] == 'X')
        {
            buffer[i] = 'Y';
            prefactor *= II;
        }
        else if(p1[i] == 'X' and p2[i] == 'Z')
        {
            buffer[i] = 'Y';
            prefactor *= -II;
        }
    }
    return pair(prefactor, std::string(buffer));
}

vector<string> build_all_paulis_N4()
{
    int ix=0;
    char p_arr[] = {'I','X','Y','Z'};
    char buffer[] = {'I','I','I','I','\0'};

    vector<string> paulis;
    for(int i=0; i < 4; i++)
    {
        buffer[0] = p_arr[i];
        for(int j=0; j < 4; j++)
        {
            buffer[1] = p_arr[j];
            for(int k=0; k < 4; k++)
            {
                buffer[2] = p_arr[k];
                for(int l=0; l < 4; l++)
                {
                    buffer[3] = p_arr[l];
                    string s = std::string(buffer);
                    paulis.push_back(s);
                    //if(ix++ > 10)
                        //return paulis;
                }
            }
        }
    }
    return paulis;
}


void check_pv_commutator(vector<string> allpaulis)
{
    // Build full vector
    cdouble coef(1.0, sqrt(2.0)-1.0);
    vector<cdouble> coefs1;
    PauliVector pv1;
    for(auto const& pauli: allpaulis)
    {
        //cout << pauli << endl;
        pv1 += coef * Pauli(pauli);
        coefs1.push_back(coef);
        coef *= M_PI;
        coef -= complex(1.0*int(coef.real()),1.0*int(coef.imag()));
    }

    // Second vector, not the same as the first
    PauliVector pv2;
    vector<cdouble> coefs2;
    for(auto const& pauli: allpaulis)
    {
        //cout << pauli << endl;
        pv2 += coef * Pauli(pauli);
        coefs2.push_back(coef);
        coef *= M_PI;
        coef -= complex(1.0*int(coef.real()),1.0*int(coef.imag()));
    }
    PauliVector pv_comm = commute(pv1, pv1);
    EXPECT_EQ( pv_comm.data.size(), 0);


    // Build the reference vector
    PauliVector pv_reference;
    for(const auto& [Pauli1, coef1] : pv1.data)
    for(const auto& [Pauli2, coef2] : pv2.data)
    {
        std::ostringstream stream1;
        stream1 << Pauli1;
        string s1 = stream1.str();

        std::ostringstream stream2;
        stream2 << Pauli2;
        string s2 = stream2.str();


        // Count the number of not-same Paulis that are also
        // not the identity
        int count = 0 ;
        for(int i=0; i < Nq; i++)
        {
            if(s1[i] == '-' or s2[i] == '-')
                continue;

            if(not(s1[i] == s2[i]))
                count++;
        }
        //cout << "Count: " << count << endl;
        if(count % 2 == 0) // They commuted, so do nothing
            continue;



        pair<cdouble, string> newterm = 
            reference_pauli_product(s1, s2);

        //cout << s1 << " " << s2 << endl;
        //cout << newterm.second << endl;

        pv_reference += 2. * coef1 * coef2 * newterm.first * Pauli(newterm.second);
    }

    pv_comm = commute(pv1, pv2);

    // Now compare
    for(auto const& [pauli, coef] : pv_reference.data)
    {
        cdouble coef2 = pv_comm.data[pauli];
        EXPECT_EQ(coef, coef2);
        //std::cout << pauli << " " << coef << " " << endl;
                  
    }
}

TEST( pauli_vector_test, Commutator)
{
    setNq(4);
    auto allpaulis = build_all_paulis_N4();
    check_pv_commutator(allpaulis);
}
