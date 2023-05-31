// (C)  Copyright Alexander (Lex) Kemper akemper@ncsu.edu 2022

#include "cartan/Pauli.h"
#include "cartan/PauliVector.h"
//#include "cartan/algebra_types.h"
//#include "cartan/models.h"
//#include "cartan/involutions.h"

#include<utility>
    using namespace std;

typedef pair<Pauli, cdouble> singlePauli;

int main()
{
    setNq(4);


    
    PauliVector pv1 = 0.5 * Pauli("IIII");
    pv1 += PauliVector( Pauli("XXII"), 0.33);
    pv1 += PauliVector( Pauli("IXYZ"), -0.25*II);

    cout << "Vector 1:\n" << pv1 << endl;

    PauliVector pv2 = 1.0 * Pauli("YIII");
    pv2 += 0.1 * Pauli("IIIX");
    cout << "Vector 2:\n" << pv2 << endl;

    PauliVector pv3 = commute(pv1, pv2);
    cout << "[pv1, pv2]:\n" << pv3 << endl;

    PauliVector pv4 = commute(pv2, pv1);
    cout << "[pv2, pv1]:\n" << pv4 << endl;

    return 0;

}
