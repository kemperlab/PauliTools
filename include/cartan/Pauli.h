// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022


#pragma once
#include<string>
#include<iostream>

typedef uint64_t pauli_int;

class Pauli
{
    public:

    // Constructors
    Pauli();
    Pauli(pauli_int init);
    Pauli(std::string init);

    // Destructor
    ~Pauli() {}

    // Modifier function
    void pauli_setup(pauli_int init);

    // Comparator
    bool operator<(const Pauli& otherPauli) const;

    // Involutions
    int countX() const;
    int countY() const;
    int countZ() const;
    bool parityX() const;
    bool parityY() const;
    bool parityZ() const;

    pauli_int code;
};


// Global vars necessary for Pauli ops
inline int Nq = -1;
inline long a_mask = -1;
inline long b_mask = -1;
inline bool masks_assigned = false;


 // Pile of helpers
void setNq(int Nq);
bool isValidPauli(char s);
pauli_int char_to_binary_rep(const char* paulistring);
pauli_int string_to_binary_rep(const std::string paulistring);
std::ostream& operator<< (std::ostream& os, const Pauli& p);
void dump_pauli(const pauli_int encoded_pauli, std::ostream& out);
void print_int_as_binary(pauli_int a, int maxbits);
void set_masks();
pauli_int symplectic_binary_comm(const pauli_int c1, const pauli_int c2);


int countY(pauli_int c);
int countZ(pauli_int c);
int countX(pauli_int c);



