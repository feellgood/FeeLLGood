#include "fem.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <string>

void Fem::infos(void) const
{
std::cout << "This is feeLLGood SHA1= " + std::string(SHAnumber) << std::endl;
msh.infos();
}




