#include "fem.h"

void Fem::affichage(void)
{
std::cout << "This is feeLLGood SHA1= " + std::string(SHAnumber) << std::endl;
std::cout << "diam bounding box ="<< diam << std::endl;
std::cout << "\t noeuds\t\t\t" << node.size() << std::endl;
std::cout << "\t faces\t\t\t" << fac.size() << std::endl;
std::cout << "\t tetraedres\t\t" << tet.size() << std::endl;
std::cout << "\t Total surface\t\t"  << surf << std::endl;
std::cout << "\t Total volume\t\t\t" << vol << std::endl;
}
