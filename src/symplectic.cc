#include <iostream>
#include <utility>
#include <complex>
    using namespace std;


#include "cartan/Pauli.h"
#include "cartan/PauliVector.h"

pauli_int symplectic_binary_comm(const pauli_int c1, const pauli_int c2)
{
    pauli_int a1 = (c1 & a_mask);
    pauli_int b1 = (c1 & b_mask) >> Nq;
    pauli_int a2 = (c2 & a_mask);
    pauli_int b2 = (c2 & b_mask) >> Nq;

    
    /*
    cout << "a1:a2 = " ;
    print_int_as_binary(a1,Nq) ; cout << ":";
    print_int_as_binary(a2,Nq) ; cout << "  ";
    cout << "b1:b2 = " ;
    print_int_as_binary(b1,Nq) ; cout << ":";
    print_int_as_binary(b2,Nq) ; cout << endl;
    
    cout << "a1 & b2 = ";
    print_int_as_binary(a1 & b2,Nq);
    cout << endl;

    cout << "a2 & b1 = ";
    print_int_as_binary(a2 & b1,Nq);
    cout << endl;
    */


    if(__builtin_parity(a1 & b2) ^ __builtin_parity(a2 & b1))
        return c1 ^ c2;
    else
        return 0;
}




pair<pauli_int,cdouble> symplectic_binary_product(const pauli_int c1, const pauli_int c2)
{
    /*
     * Work out the product with the sign included
     * A B = C
     * The terms that have a negative sign are:
     *     Y X = -iZ 
     *     Z Y = -iX 
     *     X Z = -iY 
     * the X component is a, the Z component is b. Printed with the
     * convention above, that is written b|a. With that,
     * the three cases above are
     *     1|1  0|1  =  - 1|0
     *     1|0  1|1  =  - 0|1
     *     0|1  1|0  =  - 1|1
     * in logic gates,
     *   ( b1 & a1) & (~b2 & a2)
     *   ( b1 &~a1) & ( b2 & a2)
     *   (~b1 & a1) & ( b2 & ~a2)
     *
     * These result in either a -1 (per qubit). The total number of
     * negative signs is the parity of this operation.
     *
     * The total number of Is is the number of not-equal Paulis
     */



    const int N = Nq;


    pauli_int a1 = (c1 & a_mask);
    pauli_int b1 = (c1 & b_mask) >> N;
    pauli_int a2 = (c2 & a_mask) ;
    pauli_int b2 = (c2 & b_mask) >> N;

//#define verbose
    
    
#ifdef verbose
    cout << "Going into product." << endl;
    cout << "op1:\t";
    print_int_as_binary(a1, Nq);
    cout << "|";
    print_int_as_binary(b1, Nq);
    cout << endl;
    cout << "op2:\t";
    print_int_as_binary(a2, Nq);
    cout << "|";
    print_int_as_binary(b2, Nq);
    cout << endl;

    cout << "(a1 & ( b1) & a2 & ~b2) : ";
    print_int_as_binary((a1 & ( b1) & a2 & ~b2), Nq);
    cout << endl;
    cout << "(~a1 & b1 & (a2) & b2) : ";
    print_int_as_binary((~a1 & b1 & (a2) & b2), Nq);
    cout << endl;
    cout << "((a1) & ~b1 & ~a2 & (b2)) : ";
    print_int_as_binary(((a1) & b1 & ~a2 & (b2)), Nq);
    cout << endl;
    cout << "a1 & a2 & (~b1) & (~b2): "; // X equality
    print_int_as_binary(a1 & a2 & (~b1) & (~b2), Nq);
    cout << endl;
    cout << "a1 & a2 & b1 & b2: ";      // Y equality
    print_int_as_binary(a1 & a2 & b1 & b2, Nq);
    cout << endl;
    cout << "(~a1) & (~a2) & (b1 & b2) : "; // Z equality
    print_int_as_binary((~a1) & (~a2) & (b1 & b2), Nq);
    cout << endl;
#endif   
    
    

    pauli_int product = c1 ^ c2;


    pauli_int a3 = (product & a_mask);
    pauli_int b3 = (product & b_mask) >> N;

    
    //int p1 = __builtin_popcount(a1 | b1);
    //int p2 = __builtin_popcount(a2 | b2);
    //int p3 = __builtin_popcount(a3 | b3);

    // Count overlapping Paulis
    int p4 = __builtin_popcount((a1 | b1) & (a2 | b2));

    pauli_int sign_counter = (a1 & (b1) & a2 & ~b2) |    // Y X = -iZ
                       (~a1 & b1 & (a2) & b2) |     // Z Y = -iX
                       ((a1) & ~b1 & ~a2 & (b2)) ;  // X Z = -iY
    pauli_int equal_paulis = 
                       a1 & a2 & (~b1) & (~b2) | // X equality
                       (~a1) & (~a2) & b1 & b2 | // Y equality
                       a1 & a2 & b1 & b2;        // Z equality

#ifdef verbose
    cout << "op3:\t";
    print_int_as_binary(a3, Nq);
    cout << "|";
    print_int_as_binary(b3, Nq);
    cout << endl;
    
    cout << "Sign counter: ";
    print_int_as_binary(sign_counter,2*Nq);
    cout << " " << __builtin_parity(sign_counter) << endl;
                       
    
    cout << "Paulis in p1 : " << p1 << endl;
    cout << "Paulis in p2 : " << p2 << endl;
    cout << "Paulis in p3 : " << p3 << endl;
    cout << "Overlapping Paulis : " << p4 << endl;
    cout << "Negative signs: " << __builtin_popcount(sign_counter) << endl;
    cout << "Sign parity: " << __builtin_parity(sign_counter) << endl;
    cout << "Equal Paulis: " << __builtin_popcount(equal_paulis) << endl;
#endif   
    

    int icount = (p4 -  __builtin_popcount(equal_paulis)) % 4;
    while(icount < 0)
        icount += 4;
    while(icount >= 4)
        icount -= 4;

    const cdouble sign_arr[4] = {1, II, -1, -II};
    cdouble sign = sign_arr[icount];


                      
    if(__builtin_parity(sign_counter))
        sign *= -1;

    return pair(product, sign);
}


pair<pauli_int, cdouble> symplectic_binary_commute_with_prefactor(
        const pauli_int c1, const pauli_int c2)
{    
    pauli_int a1 = (c1 & a_mask);
    pauli_int b1 = (c1 & b_mask) >> Nq;
    pauli_int a2 = (c2 & a_mask);
    pauli_int b2 = (c2 & b_mask) >> Nq;

    if(not (__builtin_parity(a1 & b2) ^ __builtin_parity(a2 & b1))) // They commute!
    {
        return pair(0, 0);
    }

    // They don't commute, so return 2x product
    pair<pauli_int, cdouble> product = symplectic_binary_product(c1, c2);
    product.second *= 2;
    return product;
}
