#ifndef CHRONOMETER_H
#define CHRONOMETER_H

#include<chrono>
#include <string>
#include <iomanip>
#include <sstream>

std::string date(void);

/** \class chronometer
convenient class to measure durations in ms or µs. Constructor is initialized with a starting value.
Then each successive calls to millis() or micros() are measuring the elapsed time since last call.
The chronometer can be reset calling reset()
*/

class chronometer
    {
    public:
        /** default constructor, 0 digit after comma */
        chronometer():nb(0) { reset(); }
 
        /** constructor, nb is the number of digit after comma */
        chronometer(const int _nb):nb(_nb) { reset(); }
    
        /** returns number of elapsed seconds since last call in double precision */
        double fp_elapsed()
            {
            start = end;
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double,std::ratio<1>> d = end-start;
            return d.count();
            }

        /** stringify number of seconds */
        std::string convertSeconds(double d)
            {
            std::stringstream ss;
            ss << d << " s";
            if (d >= 3600)
                ss << " (" << d / 3600 << " h)";
            else if (d >= 60)
                ss << " (" << d / 60 << " min)";
            ss << '\n';
            return ss.str();
            }

        /** elapsed time since last call in milliseconds */
        std::string millis() { return measure<std::milli>(" ms"); }
        
        /** elapsed time since last call in microseconds */
        std::string micros() { return measure<std::micro>(" µs"); }
    
        /** reset the chronometer */
        void reset(void)
            { end = start = std::chrono::high_resolution_clock::now(); };
    private:
        /** number of digit after comma */
        const int nb;
        
        /** time value */
        std::chrono::time_point<std::chrono::high_resolution_clock> start;
        
        /** time value */
        std::chrono::time_point<std::chrono::high_resolution_clock> end;
    
        /**
        template function that measure elapsed time since last call
        */
        template<typename T>
        std::string measure(const std::string str)
            {
            start = end;
            end = std::chrono::high_resolution_clock::now();
            
            std::chrono::duration<double,T> d = end-start;
            
            std::stringstream ss;
            ss << std::fixed << std::setprecision(nb) << d.count() << str;
            return ss.str();
            }
    };

#endif
