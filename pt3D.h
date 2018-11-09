#ifndef pt3D_h
#define pt3D_h

/** \file pt3D.h
 * \brief class pt2D and pt3D and their algebra
 * header containing pt2D and pt3D class, some functions to operate algebric operations;  << >> and
 * a template to write text files from std::vector<>
 */


#include <iostream>
#include <cstdlib>
#include <vector>
#include <math.h>
#include <fstream>
#include <algorithm>

/**
 \namespace Pt
 to grab altogether some dedicated functions and some enum for class pt2D and pt3D
 */
namespace Pt
{
    const int DIM = 3;/**< space dimension, usefull for many tables */
    
    enum index /**< convenient enum mainly to avoid direct indices values to specify coordinates in calculations */
        {IDX_X = 0,IDX_Y = 1,IDX_Z = 2};

/**
sign(x) = 1.0 if x>0, otherwise -1.0
*/
inline double sign(double x) {if (x>0) return 1.0; else return -1.0; }

/**
\return \f$ x^2 \f$
*/
inline double sq(double x) {return x*x;}

/** \class pt2D
\brief convenient class to perform geometrical operations on \f$ \mathbb{R}^2 \f$ points
 */
class pt2D
{
    
public:
    /**
     * default constructor : initialization by zero values
     */
    inline pt2D() {_x = _y = 0.0;}
    
    /**
     * constructor by values
     */
    inline pt2D(double a,double b) {_x=a; _y=b;}

    /**
     * returns a unit vector built by coordinate index, usefull to buid basis
     * example : pt2D p = pt2D(IDX_Y); //p will contain (x,y) = (0.0,1.0)
     */
    inline pt2D(enum index idx) {_x = _y = 0.0;if(idx==0) {_x=1.0;} else {_y=1.0;}}
    
    /**
     * constructor copy
    */
     inline pt2D(const pt2D & p) { _x=p.x(); _y=p.y(); }
    
    /**
     * getter for x coordinate : example : double x0 = pt.x();
     */
    inline double x(void) const {return _x;}
    
    /**
     * getter for y coordinate  : example : double y0 = pt.y();
     */
    inline double y(void) const {return _y;}
    
    /**
    * setter for x coordinate : example : pt.x(0.707);
     */
     inline void x(double a) {_x = a;}
    
    /**
     * setter for y coordinate : example : pt.y(3.14);
     */
    inline void y(double b) {_y = b;}
    
    /**
     * setter/getter by index : example : p(IDX_X) = 3.14; is equivalent to p.x(3.14);
     */
    inline double& operator() (const unsigned i) { if(i==0) {return _x;} else return _y; }
    
    /**
     * algebric += components by components
     */
    inline pt2D& operator+=(const pt2D& a) { _x += a.x(); _y += a.y(); return *this; }
    
    /**
     * algebric -= components by components
     */
    inline pt2D& operator-=(const pt2D& a) { _x -= a.x(); _y -= a.y(); return *this; }
    
    /**
     * algebric /= , if a is zero do nothing and send a message on cerr
     */
    inline pt2D& operator/=(const double& a)
    { if (a!=0) { _x /= a;_y /= a; } else std::cerr << "division by zero in pt2D::operator/=" << std::endl;
    return *this; }
    
    /**
     * algebric 2D normalization : divide each components x and y by the norm R^2 in place
     */
    inline void normalize(void)
    { double r = sqrt(_x*_x + _y*_y);
        if (r > 0.0) { _x /= r; _y /= r; }
    return; }
    
    /**
     * printing function, called by operator<<
     */
    inline void affiche(std::ostream &flux) const { flux << _x << "\t" << _y ;}
    
private:
    double _x;/**< rectangular X coordinate */
    double _y;/**< rectangular Y coordinate */
};

/** operator<< for pt2D */
inline std::ostream &operator<<(std::ostream &flux, pt2D const& p) { p.affiche(flux);return flux; }

    /** operator>> for pt2D */
inline std::istream &operator>>(std::istream &flux,pt2D & p) { double a,b; flux >> a >> b; p.x(a); p.y(b); return flux;}

/**
 * algebra : +
 */
inline pt2D operator+(pt2D const& a,pt2D const& b) { pt2D r = pt2D(a); r += b; return r; }

/**
 * algebra : -
 */
inline pt2D operator-(pt2D const& a,pt2D const& b) { pt2D r = pt2D(a); r -= b; return r; }

/**
* algebra : left scalar product
*/
 inline pt2D operator*(double const& a,pt2D const& b) { return pt2D(a*b.x(),a*b.y()); }

/**
 * algebra : division by a scalar on the right
 */
inline pt2D operator/(pt2D const & a,double const& b) { return pt2D(a.x()/b,a.y()/b); }

/**
 * algebra : euclidian scalar product
 */
inline double pScal(pt2D const& a,pt2D const& b) { return( a.x()*b.x() + a.y()*b.y() ); }

/** 
 * algebra : returns euclidian norm <br>
 The norm of \f$(x,y)\f$ is \f$\sqrt{x^2+y^2}\f$.
 */
inline double norme(pt2D const & p) {return sqrt(p.x()*p.x() + p.y()*p.y() );}



/** \class pt3D
 \brief pt3D is a class to perform algebric operations on \f$ \mathbb{R}^3 \f$ points
 all usefull operators are defined : += ; -= ; + - * <br>
 \f$ \mathbb{R}^3 \f$ scalar product and norm are defined, return double
 */
class pt3D
{
public:
    /**
     * default constructor : initialization by zero values
     */
    inline pt3D() { _x = _y = _z = 0.0; }
    
    /**
     * constructor by values
     */
    inline pt3D(double a,double b,double c) : _x(a),_y(b),_z(c) {}
    
    
    /**
     * returns a unit vector built by coordinate index, usefull to buid basis <br>
     * example : pt3D p = pt3D(IDX_Y); <br>
     p will contain \f$ (x,y,z) = (0.0,1.0,0.0) \f$
     */
    inline pt3D(enum index idx)
        {_x = _y = _z = 0.0;if(idx==IDX_X) {_x=1.0;} else {if(idx==IDX_Y) {_y=1.0;} else _z=1.0;}}
    
    
     inline pt3D(const pt3D & p) { _x=p.x(); _y=p.y(); _z=p.z(); }
    /**< constructor copy */
    
    
    /**
     * getter for x coordinate : example : double x0 = pt.x();
     */
    inline double x(void) const {return _x;}
    
    /**
     * getter for y coordinate : example : double y0 = pt.y();
     */
    inline double y(void) const {return _y;}
    
    /**
     * getter for z coordinate : example : double z0 = pt.z();
     */
    inline double z(void) const {return _z;}
    
    /**
     * \return \f$ \rho \f$ in cylindrical coordinates \f$ (\rho,\theta,z) \f$
     */
    inline double rho(void) const {return sqrt(_x*_x + _y*_y);}
    
    /**
     \return \f$ \theta \f$ in cylindrical cordinates \f$ (\rho,\theta,z) \f$
     */
    inline double theta(void) const {return atan2(_y,_x);}
    
    /**
     * setter for x coordinate : example : pt.x(0.707);
     */
    inline void x(double a) {_x = a;}
    
    /**
     * setter for y coordinate : example : pt.y(0.707);
     */
    inline void y(double b) {_y = b;}
    
    /**
     * setter for z coordinate : example : pt.z(0.707);
     */
    inline void z(double c) {_z = c;}
    
    /**
     * setter/getter by index
     @param i is an unsigned (0|1|2) use of IDX_X IDX_Y IDX_Z is recommended
     example : p(IDX_X) = 3.14; is equivalent to p.x(3.14);
     */
    inline double& operator() (unsigned i) { if(i==IDX_X) {return _x;} else {if(i==IDX_Y) {return _y;} else return _z;} }
    
    
    inline pt3D& operator=(pt3D const& p) {_x = p.x(); _y = p.y(); _z = p.z(); return *this;} /**< operator= */
    
    /**
     * algebric += components by components
     */
    inline pt3D& operator+=(const pt3D& a) { _x += a.x(); _y += a.y(); _z += a.z(); return *this; }
    
    /**
     * algebric -= components by components
     */
    inline pt3D& operator-=(const pt3D& a) { _x -= a.x(); _y -= a.y(); _z -= a.z(); return *this; }
    

/** algebric *= with a double */
    inline pt3D& operator*=(const double& a) { _x *= a;_y *= a;_z *= a; return *this; }    

/**
     * algebric /= , if a is zero do nothing and send a message on cerr
     */
    inline pt3D& operator/=(const double& a)
    { if (a!=0) { _x /= a;_y /= a;_z /= a; }
    else std::cerr << "division by zero in pt3D::operator/=" << std::endl;
    return *this; }
    
	/**
     * return algebric norm \f$ \mathcal{R}^3 \f$
     */
	inline double norm(void) {return sqrt(_x*_x + _y*_y + _z*_z);}

    /**
     * algebric 3D normalization : divide each components x,y and z by the norm \f$ \mathcal{R}^3 \f$ in place
     */
    inline void normalize(void)
    { double r = sqrt(_x*_x + _y*_y + _z*_z);
    if (r > 0.0) { _x /= r; _y /= r; _z /= r; }
    return; }
    
    /**
     * printing function, called by operator<<
     */
    inline void affiche(std::ostream &flux) const { flux << _x << "\t" << _y << "\t" << _z ; }

    /**
     scaling factor mainly for writing files in user units
     */
    inline void rescale(double scaling) {_x *= scaling; _y *= scaling; _z *= scaling; }
   
    /** \return max length coordinate */
    inline double maxLength(void) {return std::max(_z,std::max(_y,_x));}   
   
private:
    double _x;/**< rectangular X coordinate */
    double _y;/**< rectangular Y coordinate */
    double _z;/**< rectangular Z coordinate */
};

/** operator<< for pt3D, coordinates are tab separated */
inline std::ostream &operator<<(std::ostream &flux, pt3D const& p) { p.affiche(flux);return flux; }

/** operator>> for pt3D */
inline std::istream &operator>>(std::istream &flux,pt3D & p) {flux >> p(IDX_X) >> p(IDX_Y) >> p(IDX_Z); return flux;}

/**
 * algebra : +
 */
inline pt3D operator+(pt3D const& a,pt3D const& b) { pt3D r = pt3D(a); r += b; return r; }

/**
 * algebra : -
 */
inline pt3D operator-(pt3D const& a,pt3D const& b) { pt3D r = pt3D(a); r -= b; return r; }

/**
 * algebra :  vector product
 */
inline pt3D operator*(pt3D const& a,pt3D const& b)
{ return pt3D(a.y()*b.z() - a.z()*b.y() , a.z()*b.x() - a.x()*b.z() , a.x()*b.y() - a.y()*b.x()); }

/**
 * algebra : left scalar product
 */
inline pt3D operator*(double const& a,pt3D const& b) { return pt3D(a*b.x(),a*b.y(),a*b.z()); }

/**
 * algebra : right scalar product
 */
inline pt3D operator*(pt3D const& a,double const& b) { return pt3D(b*a.x(),b*a.y(),b*a.z()); }

/**
 * algebra : division components by components by a scalar to the right
 */
inline pt3D operator/(pt3D const & a,double const& b) { return pt3D(a.x()/b,a.y()/b,a.z()/b); }

/**
algebra : R^3 scalar product
 */
inline double pScal(pt3D const& a,pt3D const& b) { return( a.x()*b.x() + a.y()*b.y() + a.z()*b.z() ); }


/**
algebra : R^3 euclidian distance
 */
inline double dist(pt3D const& a,pt3D const& b) { return sqrt( sq(a.x()-b.x()) + sq(a.y()-b.y()) + sq(a.z()-b.z()) ); }


/**
algebra : returns square of the R^3 norm
 */
inline double norme2(pt3D const & p) {return( p.x()*p.x() + p.y()*p.y() + p.z()*p.z() );}

/**
 * algebra : returns R^3 norm
 */
inline double norme(pt3D const & p) {return sqrt(norme2(p));}

/** 
\return \f$ \det M \f$
*/
inline double det(const double M[DIM][DIM])
{return( M[0][0]*(M[1][1]*M[2][2] - M[1][2]*M[2][1]) 
    - M[0][1]*(M[1][0]*M[2][2] - M[1][2]*M[2][0]) 
    + M[0][2]*(M[1][0]*M[2][1] - M[1][1]*M[2][0]) );}

/** compute in place \f$ M^{-1} \f$, detM must be non zero  */
inline void inverse(double M[DIM][DIM],double detM) 
{
    double m00 = M[0][0];
    double m01 = M[0][1];
    double m02 = M[0][2]; 
    double m10 = M[1][0];
    double m11 = M[1][1];
    double m12 = M[1][2];
    double m20 = M[2][0];
    double m21 = M[2][1];
    double m22 = M[2][2];

    M[0][0] = (m11*m22 - m12*m21)/detM;
    M[0][1] = (m02*m21 - m01*m22)/detM;
    M[0][2] = (m01*m12 - m02*m11)/detM;
    M[1][0] = (m12*m20 - m10*m22)/detM;
    M[1][1] = (m00*m22 - m02*m20)/detM;
    M[1][2] = (m02*m10 - m00*m12)/detM;
    M[2][0] = (m10*m21 - m11*m20)/detM;
    M[2][1] = (m01*m20 - m00*m21)/detM;
    M[2][2] = (m00*m11 - m01*m10)/detM;
}    

/**
returns a random vector with each componant inferior to 1
*/
inline pt3D rand(void) {return pt3D(std::rand() / (RAND_MAX+1.), std::rand() / (RAND_MAX+1.), std::rand() / (RAND_MAX+1.));}

/** 
 geometry : returns iso-barycenter
*/
inline pt3D barycentre(pt3D a,pt3D b,pt3D c) { return( (a+b+c)/3.0 );}


/**
 * function template to write into a file the vector<typename>, as a column of text : use overloaded operator<<
 */
template <class T>
int writeToFile(std::string fileName,std::vector<T> v)
    { using namespace Pt; // pour le bon operator<<
    std::ofstream f_out( fileName ,std::ios::out);
    
    std::for_each(v.begin(),v.end(),[&f_out] (T &p) {f_out << p << std::endl;} );
    f_out.close();
    return 0;
    }

} // fin namespace Pt

#endif /* pt3D_h */
