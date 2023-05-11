#include <iostream>
#include <fstream>
#include <cstring>

#include "tags.h"

void on_fail_msg_error(std::ifstream &f_in, const std::string strWhat)
    {
    if (f_in.fail())
        {
        std::cerr << strWhat << ": " << strerror(errno) << std::endl;
        exit(1);
        }
    }

// carefull ! if strWhat tag contains a space it does not work ...
bool tags::lookFor(const bool _b, std::ifstream &f_in, const std::string strWhat)
    {
    std::string symb = "";
    while ((f_in.peek() != EOF) && (symb != strWhat))
        {
        f_in >> symb;
        }

    if (_b) on_fail_msg_error(f_in, "could not find tag " + strWhat);

    return !(f_in.fail());
    }

