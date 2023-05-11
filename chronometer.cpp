#include "chronometer.h"

std::string date(void)
    {
    std::stringstream realWorldTime;
    std::time_t _t = std::time(nullptr);
    realWorldTime << std::put_time(std::localtime(&_t), "%FT%H:%M:%S%z");

    return realWorldTime.str();
    }
