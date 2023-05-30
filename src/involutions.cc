// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022


#include "cartan/involutions.h"

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

