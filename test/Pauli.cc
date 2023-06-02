// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022

#include <gtest/gtest.h>
#include <algorithm>
#include <string>
#include <vector>
#include <complex>
    using namespace std;

#include "cartan/Pauli.h"
#include "cartan/symplectic.h"

typedef complex<double> cdouble;
static cdouble II(0,1);

/*
inline int Nq = 3;
inline long a_mask;
inline long b_mask;
inline bool masks_assigned = false;
*/




vector<string> build_all_paulis_N4()
{
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
                }
            }
        }
    }
    return paulis;
}


void check_involutions_3qubit()
{
    setNq(3);

    char p_arr[] = {'I','X','Y','Z'};
    char buffer[] = {'I','I','I','\0'};

    for(int i=0; i < 4; i++)
    {
        buffer[0] = p_arr[i];
        for(int j=0; j < 4; j++)
        {
            buffer[1] = p_arr[j];
            for(int k=0; k < 4; k++)
            {
                buffer[2] = p_arr[k];
                string s = std::string(buffer);
                //std::cout << buffer << " " << s << endl;

                std::string::difference_type nZ = std::count(s.begin(), s.end(), 'Z'); 
                std::string::difference_type nY = std::count(s.begin(), s.end(), 'Y'); 
                std::string::difference_type nX = std::count(s.begin(), s.end(), 'X'); 

                int nX_pauli =  Pauli(s).countX();
                int nY_pauli =  Pauli(s).countY();
                int nZ_pauli =  Pauli(s).countZ();

                EXPECT_EQ( nX, nX_pauli );
                EXPECT_EQ( nY, nY_pauli );
                EXPECT_EQ( nZ, nZ_pauli );

                nX_pauli =  countX(Pauli(s).code);
                nY_pauli =  countY(Pauli(s).code);
                nZ_pauli =  countZ(Pauli(s).code);

                EXPECT_EQ( nX, nX_pauli );
                EXPECT_EQ( nY, nY_pauli );
                EXPECT_EQ( nZ, nZ_pauli );
           
            }
        }
    }
}

pair<cdouble, string> reference_pauli_product(string p1, string p2)
{
    char buffer[Nq+1];
    buffer[Nq] = '\0';
    cdouble prefactor = 1;
    for(int i=0; i < Nq; i++)
    {
        if(p1[i] == 'I')
        {
            buffer[i] = p2[i];
            continue;
        }

        if(p2[i] == 'I')
        {
            buffer[i] = p1[i];
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

void check_symplectic_comm(vector<string>& paulis)
{
    setNq(4);

    for(const auto& p1: paulis)
    for(const auto& p2: paulis)
    {
        Pauli P1(p1);
        Pauli P2(p2);
        /*
        cout << "Checking comm: " << P1 << " " << P2 << endl;
        cout << "Checking comm: " << p1 << " " << p2 << endl;
        cout << "P1: ";
        print_int_as_binary(P1.code, 2*Nq);
        cout << endl;
        cout << "P2: ";
        print_int_as_binary(P2.code, 2*Nq);
        cout << endl;
        */

        // Count the number of not-same Paulis that are also
        // not the identity
        int count = 0 ;
        for(int i=0; i < Nq; i++)
        {
            if(p1[i] == 'I' or p2[i] == 'I')
                continue;

            if(not(p1[i] == p2[i]))
                count++;
        }

        pauli_int comm = symplectic_binary_comm(P1.code, P2.code);

        // If there's an even number then they commute
        if(count % 2 == 0)
        {
            EXPECT_EQ(comm , 0);
            continue;
        }

        // They didn't commute, so compute the product
        pair<cdouble, string> product = reference_pauli_product(p1, p2);
        EXPECT_EQ( comm, Pauli(product.second).code);
    }
}


void check_symplectic_product(vector<string>& paulis)
{
    setNq(4);

    for(const auto& p1: paulis)
    for(const auto& p2: paulis)
    {
        Pauli P1(p1);
        Pauli P2(p2);
        /*
        cout << "Checking comm: " << P1 << " " << P2 << endl;
        cout << "Checking comm: " << p1 << " " << p2 << endl;
        cout << "P1: ";
        print_int_as_binary(P1.code, 2*Nq);
        cout << endl;
        cout << "P2: ";
        print_int_as_binary(P2.code, 2*Nq);
        cout << endl;
        */

        // Count the number of not-same Paulis that are also
        // not the identity
        int count = 0 ;
        for(int i=0; i < Nq; i++)
        {
            if(p1[i] == 'I' or p2[i] == 'I')
                continue;

            if(not(p1[i] == p2[i]))
                count++;
        }

        pair<pauli_int, cdouble> comm = symplectic_binary_product(P1.code, P2.code);

        // Compute the reference product
        pair<cdouble, string> product = reference_pauli_product(p1, p2);
        EXPECT_EQ( comm.first, Pauli(product.second).code);

        // Now check the prefactor
        EXPECT_EQ( product.first, comm.second );
    }
}


TEST( pauli_test, Involutions ){
    check_involutions_3qubit();
}

TEST( pauli_test, SymplecticCommutator){
    vector<string> paulis = build_all_paulis_N4();
    check_symplectic_comm(paulis);
}


TEST( pauli_test, SymplecticProduct){
    vector<string> paulis = build_all_paulis_N4();
    check_symplectic_product(paulis);
}




/*
 *
TEST(TestSuiteName, TestName) {
  ... test body ...
}
 * 
 */
