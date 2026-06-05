#include <iostream>
#include <fstream>
#include <cstring>

#include "tags.h"

void on_fail_msg_error(const std::ifstream &f_in, const std::string& strWhat)
    {
    if (f_in.fail())
        {
        std::cerr << strWhat << ": " << strerror(errno) << "\n";
        exit(1);
        }
    }
