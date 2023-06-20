#ifndef vecND_h
#define vecND_h

#include <vector>
#include <algorithm>
#include <execution>

/**
* \brief class vecND for small dense vectors, using parallelization and vectorization for operators
template class for double precision vectors of dimension _DIM: template parameter
*/


//operators are twice slower with par_unseq than unseq 
#define strategy std::execution::unseq

template <int _DIM>
class vecND
    {
    public:
        
        /**
        constructor, all coefficients are set to zero
        */
        vecND(void)
            {
            std::fill_n(std::execution::unseq, std::begin(_x), _DIM, 0);
            }

        /**
        coefficient setter
        */
        void set(std::initializer_list<double> const & _v)
            {
            //slightly faster than std::copy
            std::copy_n(std::execution::unseq, _v.begin(), _DIM, std::begin(_x));
            }
            
        /** getter by int : carefull, no test on i ! */
        double operator()(const int i) const { return _x[i]; }
        
        /** setter by int : carefull, no test on i ! */
        double & operator()(const int i) { return _x[i]; }
        
        /** setter by int : carefull, no test on i ! */
        void operator()(const int i, const double val) { _x[i] = val; }

    /**
     * algebric += components by components
     */
    inline vecND &operator+=(const vecND &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x),
                        std::begin(a._x), std::begin(_x),
                        [](double &v_loc, const double &val) -> double { return v_loc + val; } );
        return *this;
        }

    /**
     * algebric -= components by components
     */
    inline vecND &operator-=(const vecND &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x),
                        std::begin(a._x), std::begin(_x),
                        [](double &v_loc, const double &val) -> double { return v_loc - val; } );
        return *this;
        }

    /** algebric *= with a double */
    inline vecND &operator*=(const double &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x), std::begin(_x),
                        [a](double &v_loc) -> double { return v_loc*=a; } );
        return *this;
        }
    
    protected:
        //a std_vector<double> is slower : depending on policy from 25% to 100% of added time on full_test.py
        
        /**
        the vector components
        */
        double _x[_DIM];
    };

#endif
