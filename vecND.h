#ifndef vecND_h
#define vecND_h

#include <vector>
#include <algorithm>
#include <execution>

/** \file vecND.h
* \brief template class vecND for small dense vectors, using parallelization and vectorization for operators
template class for double precision vectors of dimension _DIM: template parameter
*/


//operators are twice slower with par_unseq than unseq 
/**
execution policy for various operators
*/
#define strategy std::execution::unseq

/** \class vecND
\brief vector of dimension _DIM
template class for a vector of dimension _DIM with basic operators implemented using transform algorithm 
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
     * algebra += components by components
     */
    inline vecND &operator+=(const vecND &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x),
                        std::begin(a._x), std::begin(_x),
                        [](double &v_loc, const double &val) -> double { return v_loc + val; } );
        return *this;
        }

    /**
     * algebra -= components by components
     */
    inline vecND &operator-=(const vecND &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x),
                        std::begin(a._x), std::begin(_x),
                        [](double &v_loc, const double &val) -> double { return v_loc - val; } );
        return *this;
        }

    /** algebra *= with a double */
    inline vecND &operator*=(const double &a)
        {
        std::transform(strategy, std::begin(_x), std::end(_x), std::begin(_x),
                        [a](double &v_loc) -> double { return v_loc*=a; } );
        return *this;
        }

    /**
     * algebra /= , if a is zero do nothing and send a message on cerr
     */
    inline vecND &operator/=(const double &a)
        {
        if (a != 0)
            {
            *this *= 1.0/a;
            }
        else
            std::cerr << "division by zero in vecND::operator/=" << std::endl;
        return *this;
        }
    
    /**
    * algebra : returns square of the norm \f$ \mathcal{R}^n \f$
    */
    inline double norm2(void) const
        {
        return std::transform_reduce(strategy, std::begin(_x), std::end(_x), 0.0, std::plus{},
                                     [](double val){return val*val;} );
        }

    /**
     * returns norm \f$ \mathcal{R}^n \f$
     */
    inline double norm(void) const
        {
        return sqrt(norm2());
        }

/**
algebra : euclidian \f$ \mathcal{R}^n \f$ scalar product with vector a
 */
    inline double pScal(const vecND &a) const
        {
        return std::transform_reduce(strategy, std::begin(_x), std::end(_x), std::begin(a._x), 0.0, std::plus{},
                                     [](const double v1,const double v2){return v1*v2;} );
        }

    protected:
        //a std_vector<double> is slower : depending on policy from 25% to 100% of added time on full_test.py
        
        /**
        the vector components
        */
        double _x[_DIM];
    };



#endif
