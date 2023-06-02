// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022


#include "cartan/models.h"

extern int Nq;

vector<string> heisenberg(bool pbc)
{
    const int N = Nq;

    vector<string> ham;
    char* base_string = new char[N+1];
    base_string[N] = '\0';
    for(int i=0; i < N; i++)
        base_string[i] = 'I';

    for(int i=0; i < N; i++)
    {
        int j = (i+1) % N;
        if(j==0 and not pbc) continue;

        base_string[i] = 'X';
        base_string[j] = 'X';
        ham.push_back(string(base_string));

        base_string[i] = 'Y';
        base_string[j] = 'Y';
        ham.push_back(string(base_string));

        base_string[i] = 'Z';
        base_string[j] = 'Z';
        ham.push_back(string(base_string));

        base_string[i] = 'I';
        base_string[j] = 'I';
        

    }

    return ham;
}

vector<string> tfxy(bool pbc)
{
    const int N = Nq;

    vector<string> ham;
    char* base_string = new char[N+1];
    base_string[N] = '\0';
    for(int i=0; i < N; i++)
        base_string[i] = 'I';

    for(int i=0; i < N; i++)
    {
        base_string[i] = 'Z';
        ham.push_back(string(base_string));

        base_string[i] = 'I';
    }

    for(int i=0; i < N; i++)
    {
        int j = (i+1) % N;
        if(j==0 and not pbc) continue;

        base_string[i] = 'X';
        base_string[j] = 'X';
        ham.push_back(string(base_string));

        base_string[i] = 'Y';
        base_string[j] = 'Y';
        ham.push_back(string(base_string));

        base_string[i] = 'I';
        base_string[j] = 'I';
        

    }

    return ham;
}

