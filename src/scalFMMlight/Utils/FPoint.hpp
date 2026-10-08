// See LICENCE file at project root
//
#ifndef FPOINT_HPP
#define FPOINT_HPP

#include <array>
#include <iterator>
#include <ostream>

/** 3-dimensional cartesian coordinates
 *
 * \author Berenger Bramas <berenger.bramas@inria.fr>, Quentin Khan <quentin.khan@inria.fr>
 *
 * Fixed size array that represents coordinates in space. This class adds a few convenience 
 * operations such as addition, scalar multiplication and division and formated stream output.
 *
 * @param _Real The floating number type
 **/

/** Space dimension */
const std::size_t _Dim=3;

/** template class FPoint, dimension set to three, cartesian coordinates x,y,z in a std::array */
template<typename _Real>
class FPoint : public std::array<_Real, _Dim> {
public:
    /** Floating number type */
    using FReal = _Real;

    /** space dimension = 3 */
    constexpr static const std::size_t Dim = _Dim;

private:

    /// Type used in SFINAE to authorize arithmetic types only in template parameters
    template<class T>
    using must_be_arithmetic = typename std::enable_if<std::is_arithmetic<T>::value>::type*;
    /// Type used in SFINAE to authorize floating point types only in template parameters
    template<class T>
    using must_be_floating = typename std::enable_if<std::is_floating_point<T>::value>::type*;
    /// Type used in SFINAE to authorize integral types only in template parameters
    template<class T>
    using must_be_integral = typename std::enable_if<std::is_integral<T>::value>::type*;

public:

    /** Default constructor */
    FPoint() = default;
    /** Copy constructor */
    FPoint(const FPoint&) = default;

    /** Copy constructor from other point type */
    template<typename A, must_be_arithmetic<A> = nullptr>
    FPoint(const FPoint<A>& other)
        {
        this->data()[0] = other.data()[0];
        this->data()[1] = other.data()[1];
        this->data()[2] = other.data()[2];
        }

    /** Constructor from array */
    FPoint(const FReal array_[Dim])
        {
        this->data()[0] = array_[0];
        this->data()[1] = array_[1];
        this->data()[2] = array_[2];
        }

    /** Constructor from 3 FReals */
    template<typename FReal>
    FPoint(const FReal& X, const FReal& Y, const FReal& Z)
        {
        this->data()[0] = X;
        this->data()[1] = Y;
        this->data()[2] = Z;
        }

    /** Additive constructor, same as FPoint(other + add_value) */
    FPoint(const FPoint& other, const FReal add_value)
        {
        this->data()[0] = other.data()[0] + add_value;
        this->data()[1] = other.data()[1] + add_value;
        this->data()[2] = other.data()[2] + add_value;
        }

    /** Assignment operator
     * \param other A FPoint object.
     */
    template<typename T, must_be_arithmetic<T> = nullptr>
    FPoint<FReal>& operator=(const FPoint<T>& other) {
        this->copy(other);
        return *this;
    }

    /** Sets the point value */
    template<typename FReal>
    void setPosition(const FReal& X, const FReal& Y, const FReal& Z)
        {
        this->data()[0] = X;
        this->data()[1] = Y;
        this->data()[2] = Z;
        }

    /** \brief Get x
     * \return this->data()[0]
     */
    FReal getX() const { return this->data()[0]; }

    /** \brief Get y
     * \return this->data()[1]
     */
    FReal getY() const { return this->data()[1]; }

    /** \brief Get z
     * \return this->data()[2]
     */
    FReal getZ() const { return this->data()[2]; }

    /** \brief Set x
     * \param inX the new x
     */
    void setX(const FReal inX) { this->data()[0] = inX; }

    /** \brief Set y
     * \param inY the new y
     */
    void setY(const FReal inY) { this->data()[1] = inY; }

    /** \brief Set z
     * \param inZ the new z
     */
    void setZ(const FReal inZ) { this->data()[2] = inZ; }

    /** \brief Add to the x-dimension the inX value
     * \param inX the increment in x
     */
    void incX(const FReal inX) { this->data()[0] += inX; }

    /** \brief Add to the y-dimension the inY value
     * \param inY the increment in y
     */
    void incY(const FReal inY) { this->data()[1] += inY; }

    /** \brief Add to z-dimension the inZ value
     * \param inZ the increment in z
     */
    void incZ(const FReal inZ) { this->data()[2] += inZ; }

    /** \brief Get a pointer on the coordinate of FPoint<FReal>
     * \return the data value array
     */
    FReal * getDataValue() { return this->data(); }

    /** \brief Get a pointer on the coordinate of FPoint<FReal>
     * \return the data value array
     */
    const FReal * getDataValue() const { return this->data(); }

    /** \brief Compute the distance to the origin
     * \return the norm of the FPoint
     */
    FReal norm() const { return sqrt(norm2()); }

    /** \brief Compute the distance to the origin
     * \return the square norm of the FPoint
     */
    FReal norm2() const
        { return sq(this->data()[0]) + sq(this->data()[1]) + sq(this->data()[2]); }

    /** Addition assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator +=(const FPoint<T>& other)
        {
        auto other_it = other.begin();
        auto this_it  = this->begin();
        *this_it += *other_it; this_it++; other_it++;
        *this_it += *other_it; this_it++; other_it++;
        *this_it += *other_it;
        return *this;
        }

    /** Addition operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator+(FPoint<FReal> lhs, const FPoint<T>& rhs)
        {
        lhs += rhs;
        return lhs;
        }

    /** Scalar assignment addition */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator+=(const T& val)
        {
        this->data()[0] += val;
        this->data()[1] += val;
        this->data()[2] += val;
        return *this;
        }

    /** Scalar addition */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator+(FPoint<FReal> lhs, const T& val)
        {
        lhs += val;
        return lhs;
        }

    /** Subtraction assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator -=(const FPoint<T>& other)
        {
        auto other_it = other.begin();
        auto this_it  = this->begin();
        *this_it -= *other_it; this_it++; other_it++;
        *this_it -= *other_it; this_it++; other_it++;
        *this_it -= *other_it;
        return *this;
        }

    /** Subtraction operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator-(FPoint<FReal> lhs, const FPoint<T>& rhs)
        {
        lhs -= rhs;
        return lhs;
        }

    /** Scalar subtraction assignment */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint& operator-=(const T& val)
        {
        this->data()[0] -= val;
        this->data()[1] -= val;
        this->data()[2] -= val;
        return *this;
        }

    /** Scalar subtraction */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator-(FPoint<FReal> lhs, const T& rhs)
        {
        lhs -= rhs;
        return lhs;
        }

    /** Right scalar multiplication assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint<FReal>& operator *=(const T& val)
        {
        this->data()[0] *= val;
        this->data()[1] *= val;
        this->data()[2] *= val;
        return *this;
        }

    /** Right scalar multiplication operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator*(FPoint<FReal> lhs, const T& val)
        {
        lhs *= val;
        return lhs;
        }

    /** Left scalar multiplication operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator*(const T& val, FPoint<FReal> rhs)
        {
        rhs *= val;
        return rhs;
        }

    /** Data to data division operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator/(FPoint<FReal> lhs, const FPoint<FReal>& rhs)
        {
        lhs /= rhs;
        return lhs;
        }

    /** Right scalar division assignment operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    FPoint<FReal>& operator /=(const T& val)
        {
        this->data()[0] /= val;
        this->data()[1] /= val;
        this->data()[2] /= val;
        return *this;
        }

    /** Right scalar division operator */
    template<class T, must_be_arithmetic<T> = nullptr>
    friend FPoint<FReal> operator/(FPoint<FReal> lhs, const T& val)
        {
        lhs /= val;
        return lhs;
        }

    /** Formated output stream operator */
    friend std::ostream& operator<<(std::ostream& os, const FPoint<FReal>& pos)
        {
        os << "[" << pos->data()[0] << ", " << pos->data()[1] << ", " << pos->data()[2] << "]";
        return os;
        }
};

#endif
