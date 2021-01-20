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
   std::cout << "usage : cmk file" << std::endl;
   exit(1);
   }

std::string str=argv[1];

std::vector <int> old2newlabel;
std::vector <int> new2oldlabel;
reverse_cmk(str, old2newlabel, new2oldlabel);
update_labelling(str, old2newlabel, new2oldlabel);

return 0;
}
