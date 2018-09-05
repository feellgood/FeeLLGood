#include "fem.h"

using namespace std;

void Fem::affichage(void)
{
cout << "\n\t regions\t\t" << REG << endl;
cout << "\t noeuds\t\t\t" << NOD << endl;
cout << "\t faces\t\t\t" << FAC << endl;
cout << "\t tetraedres\t\t" << TET << endl;
cout << "\t surface\t\t"  << surf << endl;
cout << "\t volume\t\t\t" << vol << endl;
}
