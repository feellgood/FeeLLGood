#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <stdexcept>

#include <boost/format.hpp>
#include <boost/tokenizer.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "rev_CMK.h"


int main(int argc, char **argv)
{
if (argc !=2) {
   std::cout << "rev_CMK expect a mesh filename input to apply reverse Cuthill McKee algorithm. The original file is unchanged, a new mesh file is created with original input filename with an extra extension .r_cmk" << std::endl;
   exit(1);
   }

std::string fileName = argv[1];

std::vector <int> old2newlabel;
std::vector <int> new2oldlabel;

std::ifstream fin(fileName);
if (!fin) {
   std::cerr << "Cannot open file " << fileName << std::endl;
   exit(1);
   }

reverse_cmk(fin, old2newlabel, new2oldlabel);
update_labelling(fileName, old2newlabel, new2oldlabel);

return 0;
}
