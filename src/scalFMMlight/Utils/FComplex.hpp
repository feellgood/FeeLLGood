// ===================================================================================
// olivier.coulaud@inria.fr, berenger.bramas@inria.fr
// This software is a computer program whose purpose is to compute the FMM.
//
// This software is governed by the CeCILL-C and LGPL licenses and
// abiding by the rules of distribution of free software.  
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public and CeCILL-C Licenses for more details.
// "http://www.cecill.info". 
// "http://www.gnu.org/licenses".
// ===================================================================================
#ifndef FCOMPLEXE_HPP
#define FCOMPLEXE_HPP

/**
* @author Berenger Bramas (berenger.bramas@inria.fr)
* @class FComplex<FReal>
* This class is a basic implementation of Complex numbers.
* Please read the license
* Do not modify the attributes of this class, it can be passed to blas fonction and has to be 2 x FReal size only.
*/

/** FComplex class */
template <class FReal>
class FComplex {
    /** Real & Imaginary container */
    FReal complex[2];

public:
    /** Default Constructor (set real and imaginary part to 0) */
    FComplex()
        {
        complex[0] = 0;
        complex[1] = 0;
        }

    /** Constructor with values
      * @param inReal the real
      * @param inImag the imaginary
      */
    explicit FComplex(const FReal inReal, const FReal inImag)
        {
        complex[0] = inReal;
        complex[1] = inImag;
        }

    /** Copy constructor */
    FComplex(const FComplex<FReal>& other)
        {
        complex[0] = other.complex[0];
        complex[1] = other.complex[1];
        }

    /** Copy operator */
    FComplex<FReal>& operator=(const FComplex<FReal>& other)
        {
        this->complex[0] = other.complex[0];
        this->complex[1] = other.complex[1];
        return *this;
        }

    /** Get imaginary part */
    FReal getImag() const{ return this->complex[1]; }

    /** Get real part */
    FReal getReal() const{ return this->complex[0]; }

    /** Set Imaginary part */
    void setImag(const FReal inImag) { this->complex[1] = inImag; }

    /** Set Real part */
    void setReal(const FReal inReal) { this->complex[0] = inReal; }

    /** set to zero real and imaginary */
    void setToZero()
        {
        this->complex[0] = FReal(0.0);
        this->complex[1] = FReal(0.0);
        }

    /** Set Real and imaginary */
    void setRealImag(const FReal inReal, const FReal inImag)
        {
        this->complex[0] = inReal;
        this->complex[1] = inImag;
        }

    
    /** Operator += with a RHS complex */
    FComplex<FReal>& operator+=(const FComplex<FReal>& other)
        {
        this->complex[0] += other.complex[0];
        this->complex[1] += other.complex[1];
        return *this;
        }

    /** Operator + with a RHS complex (unused, just for dev/debug) */
    FComplex<FReal> operator+(const FComplex<FReal>& other)
        { return FComplex<FReal>(this->complex[0] + other.complex[0], this->complex[1] + other.complex[1]); }

    /** Increment real part by RHS FReal */
   void incReal(const FReal inIncReal) { this->complex[0] += inIncReal; }

    /** Increment imaginary part by RHS FReal */
    void incImag(const FReal inIncImag) { this->complex[1] += inIncImag; }

    /** operator *= with a RHS complex */
    FComplex<FReal>& operator*=(const FComplex<FReal>& other)
        {
        const FReal tempReal = this->complex[0];
        this->complex[0] = (tempReal * other.complex[0]) - (this->complex[1] * other.complex[1]);
        this->complex[1] = (tempReal * other.complex[1]) + (this->complex[1] * other.complex[0]);
        return *this;
        }
};

#endif //FCOMPLEXE_HPP
