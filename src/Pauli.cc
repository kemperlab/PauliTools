// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022




#include<cstdint>
#include<cstring>
#include<string>
#include<iostream>
#include<cassert>
    using namespace std;



#include "cartan/Pauli.h"

    /*
extern int Nq;
extern long a_mask;
extern long b_mask;
extern bool masks_assigned;
*/





Pauli::Pauli()
{
    this->pauli_setup(0);
}

Pauli::Pauli(pauli_int init)
{
    this->pauli_setup(init);
}

Pauli::Pauli(string init)
{
    this->code = string_to_binary_rep(init);
}

void Pauli::pauli_setup(pauli_int init)
{
    assert(masks_assigned);
    this->code = init;
}

bool Pauli::operator<( const Pauli& otherPauli) const
{ 
    return this->code < otherPauli.code;
}

inline bool operator<( const Pauli& a1, const Pauli& a2 )
{
    return a1.code < a2.code;
}

void set_masks()
{    
    const int N = Nq;
    a_mask = 0;
    b_mask = 0;

    pauli_int x = 1;

    for(int i=0; i < N; i++)
    {
        a_mask ^= (x<<i);
        b_mask ^= (x<<(i+N));
    }

     cout << "a binary mask: " ;
     print_int_as_binary(a_mask,2*N);
     cout << endl;
     cout << "b binary mask: " ;
     print_int_as_binary(b_mask,2*N);
     cout << endl;

    masks_assigned = true;
}

ostream& operator<< (ostream& os, const Pauli& p)
{
    dump_pauli(p.code, cout);
    return os;
}

void dump_pauli(const pauli_int encoded_pauli, ostream& out)
{
    pauli_int x = 1;
    for(int i=0; i < Nq; i++)
    {
        bool a = encoded_pauli & (x << i);
        bool b = encoded_pauli & (x << (i+Nq));
 
        if(a and b)
            out << 'Y';
        else if(a)
            out << 'X';
        else if(b)
            out << 'Z';
        else
            out << '-';
    }
}


void setNq(int newNq)
{
    Nq = newNq;
    set_masks();
}


bool isValidPauli(char s)
{
    if(s=='I' or s=='-' or s=='X' or s=='Y' or s=='Z')
        return true;
    return false;
}


pauli_int char_to_binary_rep(const char* paulistring)
{
    bool data[2*Nq];
    bool* a = data;
    bool* b = &(data[Nq]);
    memset(data, 0, 2*Nq*sizeof(bool));

    for(int j=0; j < Nq; j++)
    {
        char pauli = paulistring[j];
        assert(isValidPauli(pauli));
        if(pauli=='X' or pauli=='Y')
            a[j] = 1;
        if(pauli=='Z' or pauli=='Y')
            b[j] = 1;
    }


    // Bool array to int
    const int M = 2 * Nq;

    pauli_int drep = 0;
    pauli_int current_idx = 1;
    for(int i=0; i < M; i++)
    {
        drep += data[i] * current_idx;
        current_idx = current_idx << 1;
    }
    return drep;

}


pauli_int string_to_binary_rep(const string paulistring)
{
    return char_to_binary_rep(paulistring.c_str());
}


bool Pauli::parityY() const
{
    return __builtin_parity( (code & a_mask) & ((code & b_mask) >> Nq) );
}

int Pauli::countY() const
{
    return __builtin_popcount( (code & a_mask) & ((code & b_mask) >> Nq) );
}


bool Pauli::parityX() const
{
    return __builtin_parity( (code & a_mask) & ~((code & b_mask) >> Nq) );
}

int Pauli::countX() const
{
    return __builtin_popcount( (code & a_mask) & ~((code & b_mask) >> Nq) );
}

bool Pauli::parityZ() const
{
    return __builtin_parity(  ~(code & a_mask) & ((code & b_mask) >> Nq) );
}

int Pauli::countZ() const
{
    return __builtin_popcount(  ~(code & a_mask) & ((code & b_mask) >> Nq) );
}

void print_int_as_binary(pauli_int a, int maxbits)
{
    pauli_int x=1;
    for(int i=maxbits-1; i >= 0; i--)
        cout << bool(a & (x << i));
}

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



// These are using the integer directly
bool parityY(pauli_int c)
{
    return __builtin_parity(  (c & a_mask) & ((c & b_mask) >> Nq) );
}

bool parityZ(pauli_int c)
{
    return __builtin_parity(  ~(c & a_mask) & ((c & b_mask) >> Nq) );
}

bool parityX(pauli_int c)
{
    return __builtin_parity(  (c & a_mask) & ~((c & b_mask) >> Nq) );
}

int countY(pauli_int c)
{
    return __builtin_popcount( (c & a_mask) & ((c & b_mask) >> Nq) );
}

int countZ(pauli_int c)
{
    return __builtin_popcount(  ~(c & a_mask) & ((c & b_mask) >> Nq) );
}

int countX(pauli_int c)
{
    return __builtin_popcount( (c & a_mask) & ~((c & b_mask) >> Nq) );
}


