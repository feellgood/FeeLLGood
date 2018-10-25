#ifndef matBlocDiag_h
#define matBlocDiag_h

/** \file matBlocDiag.h
 * \brief class matBlocDiag and their algebra
 * header containing matBlocDiag class, some functions to operate algebric operations on bloc diagonal(4x4) matrix
 */


#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>

/**
 \namespace matBlocDiag
 to grab altogether some dedicated functions and some enum for class matBlocDiag
 */
namespace matBlocDiag
{
    const int DIM = 4;/**< space intrinsic dimension, usefull for many tables */
    
    /** \enum index
     specify by integer indices the coordinate
     */
    enum index {IDX_A = 0,IDX_B = 1,IDX_C = 2,IDX_D = 3};

/**
 * sign(x) = 1.0 if x>0, otherwise -1.0
*/
inline double sign(double x) {if (x>0) return 1.0; else return -1.0; }

/**
\return \f$ x^2 \f$
*/
inline double sq(double x) {return x*x;}



/** \class BlocElem
 \brief BlocElem is a class to perform algebric operations on \f$ \mathbb{R}^4 \f$ points
 all usefull operators are defined : += ; -= ; + - <br>
 \f$ \mathbb{R}^4 \f$ scalar product and norm are defined, return double
 The purpose is to provide some convenient operators to perform algebra of 4x4 diagonal matrix, including direct product, hadamard product
  */
class BlocElem
{
public:
    /**
     * default constructor : initialization by zero values
     */
    inline BlocElem() { _x[0] = _x[1] = _x[2] = _x[3] = 0.0; }
    
    /**
     * constructor by values
     */
    inline BlocElem(double a,double b,double c,double d) { _x[0] = a; _x[1] = b; _x[2] = c; _x[3] = d;}
    
    
    /**
     * constructor for a unit vector built by coordinate index, usefull to buid basis <br>
     */
    inline BlocElem(enum index idx)
        {_x[0] = _x[1] = _x[2] = _x[3] = 0.0; _x[idx] = 1.0;}
    
    
     inline BlocElem(const BlocElem & p) { _x[0]=p.a(); _x[1]=p.b(); _x[2]=p.b(); _x[3]=p.d(); }
    /**< constructor copy */
    
    
    /**
     * getter for a
     */
    inline double a(void) const {return _x[0];}
    
    /**
     * getter for b
     */
    inline double b(void) const {return _x[1];}
    
    /**
     * getter for c
     */
    inline double c(void) const {return _x[2];}
    
    /**
     * getter for d
     */
    inline double d(void) const {return _x[3];}
    
    
    /**
     * setter for a
     */
    inline void a(double x) {_x[0] = x;}
    
    /**
     * setter for b
     */
    inline void b(double y) {_x[1] = y;}
    
    /**
     * setter for c;
     */
    inline void c(double z) {_x[2] = z;}

    /**
     * setter for d;
     */
    inline void d(double t) {_x[3] = t;}
    
    /**
     * getter by single index
     */
    inline double operator() (index i) { return _x[i]; }
    
    /**
     * getter by pair of index
     */
    inline double operator() (index i,index j) { if(i != j) {return 0.0;} else { return _x[i];} }
    
    
    inline BlocElem& operator=(BlocElem const& p) {_x[0] = p.a(); _x[1] = p.b(); _x[2] = p.c(); _x[3] = p.d(); return *this;} /**< operator= */
    
    /**
     * algebric += components by components
     */
    inline BlocElem& operator+=(const BlocElem& p) { _x[0] += p.a(); _x[1] += p.b(); _x[2] += p.c(); _x[3] += p.d(); return *this; }
    
    /**
     * algebric -= components by components
     */
    inline BlocElem& operator-=(const BlocElem& p) { _x[0] -= p.a(); _x[1] -= p.b(); _x[2] -= p.c(); _x[3] -= p.d(); return *this; }
    

/** algebric *= with a double */
    inline BlocElem& operator*=(const double& t) { _x[0] *= t;_x[1] *= t; _x[2] *= t; _x[3] *= t; return *this; }    

/**
     * algebric /= , if t is zero do nothing and send a message on cerr
     */
    inline BlocElem& operator/=(const double& t)
    { if (t!=0) { _x[0] /= t; _x[1] /= t; _x[2] /= t; _x[3] /= t;}
    else std::cerr << "division by zero in BlocElem::operator/=" << std::endl;
    return *this; }
    
	/**
     * return algebric norm \f$ \mathcal{R}^4 \f$
     */
	inline double norm(void) {return sqrt(_x[0]*_x[0] + _x[1]*_x[1] + _x[2]*_x[2] + _x[3]*_x[3]);}

    /**
     * algebric normalization : divide each components by the norm \f$ \mathcal{R}^4 \f$ in place
     */
    inline void normalize(void)
    { double r = norm();
    if (r > 0.0) { _x[0] /= r; _x[1] /= r; _x[2] /= r; _x[3] /= r; }
    return; }
    
    /**
     * printing function, called by operator<<
     */
    inline void affiche(std::ostream &flux) const { flux << _x[0] << "\t" << _x[1] << "\t" << _x[2] << "\t" << _x[3]; }

    
    
private:
    double _x[DIM];/**< rectangular \f$ x_i \f$ coordinates */
};

/** operator<< for BlocElem, coordinates are tab separated */
inline std::ostream &operator<<(std::ostream &flux, BlocElem const& p) { p.affiche(flux);return flux; }

/**
 * algebra : +
 */
inline BlocElem operator+(BlocElem const& a,BlocElem const& b) { BlocElem r = BlocElem(a); r += b; return r; }

/**
 * algebra : -
 */
inline BlocElem operator-(BlocElem const& a,BlocElem const& b) { BlocElem r = BlocElem(a); r -= b; return r; }


/**
 * algebra : left scalar product
 */
inline BlocElem operator*(double const& a,BlocElem const& p) { return BlocElem(a*p.a(),a*p.b(),a*p.c(),a*p.d()); }

/**
 * algebra : right scalar product
 */
inline BlocElem operator*(BlocElem const& p,double const& a) { return BlocElem(a*p.a(),a*p.b(),a*p.c(),a*p.d()); }

/**
 * algebra : division components by components by a scalar to the right
 */
inline BlocElem operator/(BlocElem const & p,double const& b) { return BlocElem(p.a()/b, p.b()/b, p.c()/b, p.d()/b ); }

/**
algebra : R^4 scalar product
 */
inline double pScal(BlocElem const& x,BlocElem const& y) { return( x.a()*y.a() + x.b()*y.b() + x.c()*y.c() + x.d()*y.d() ); }


/**
algebra : R^4 euclidian distance
 */
inline double dist(BlocElem const& x,BlocElem const& y) { return sqrt( sq(x.a()-y.a()) + sq(x.b()-y.b()) + sq(x.c()-y.c()) + sq(x.d()-y.d())); }


/**
algebra : returns square of the R^4 norm
 */
inline double norme2(BlocElem const & p) {return( p.a()*p.a() + p.b()*p.b() + p.c()*p.c() + p.d()*p.d() );}

/**
 * algebra : returns R^4 norm
 */
inline double norme(BlocElem const & p) {return sqrt(norme2(p));}


/**
returns a random vector with each componant inferior to 1
*/
inline BlocElem rand(void) 
    {return BlocElem(std::rand() / (RAND_MAX+1.), std::rand() / (RAND_MAX+1.), std::rand() / (RAND_MAX+1.), std::rand() / (RAND_MAX+1.));}



/**
 * function template to write into a file the vector<typename>, as a column of text : use overloaded operator<<
 */
template <class T>
int writeToFile(std::string fileName,std::vector<T> v)
    { using namespace matBlocDiag; // pour le bon operator<<
    std::ofstream f_out( fileName ,std::ios::out);
    
    std::for_each(v.begin(),v.end(),[&f_out] (T &p) {f_out << p << std::endl;} );
    f_out.close();
    return 0;
    }

    struct matBloc
    {
        static const int I=2;
        static const int J=3;
        BlocElem D[I][J];
    };
    
} // fin namespace

#endif /* matBlocDiag_h */
