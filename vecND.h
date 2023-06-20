#ifndef vecND_h
#define vecND_h

#include <vector>
#include <algorithm>
#include <execution>

/**
* \brief class vecND for small dense vectors, using parallelization and vectorization for operators
template class for double precision vectors of dimension _DIM: template parameter
*/

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
        
        /** constructor copy */
/*
        vecND(const vecND &_v)
            {
            _x.assign(_v._x.begin(),_v._x.end());
             //memcpy(_x, v._x, sizeof _x);
            }
*/
        /** operator= */
/*
        vecND &operator=(vecND const & _v)
            {
            _x.assign(_v._x.begin(),_v._x.end());
            return *this;
            }
*/
    
    protected:
        //a std_vector<double> is slower : depending on policy from 25% to 100% of added time on full_test.py
        
        /**
        the vector components
        */
        double _x[_DIM];
    };

#endif
