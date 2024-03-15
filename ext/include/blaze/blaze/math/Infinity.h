//=================================================================================================
/*!
//  \file blaze/math/Infinity.h
//  \brief Numerical infinity for built-in data types.
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZE_MATH_INFINITY_H_
#define _BLAZE_MATH_INFINITY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Platform.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/Limits.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Negative infinity for built-in data types.
// \ingroup math
//
// The NegativeInfinity class is a wrapper class around the functionality of the blaze::Limits
// class to provide the possibility to assign negative infinity values to built-in data types.
// As negative infinity value, the largest possible negative value of the corresponding data
// type is used. In order to assign the negative infinity value, the NegativeInfinity class
// can be implicitly converted to all signed integral and floating point data types:
//
// <ul>
//    <li>integers</li>
//    <ul>
//       <li>signed char, char, wchar_t</li>
//       <li>short</li>
//       <li>int</li>
//       <li>long</li>
//       <li>ptrdiff_t (for certain 64-bit compilers)</li>
//    </ul>
//    <li>floating points</li>
//    <ul>
//       <li>float</li>
//       <li>double</li>
//       <li>long double</li>
//    </ul>
// </ul>
//
// \note The NegativeInfinity class is a helper class for the Infinity class. It cannot be
// instantiated on its own, but can only be used by the Infinity class.
*/
template< typename I >  // Positive infinity type
class NegativeInfinity
{
 public:
   //**Type definitions****************************************************************************
   using PositiveType = I;  //!< The positive infinity type.
   //**********************************************************************************************

 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr NegativeInfinity();
   NegativeInfinity( const NegativeInfinity& ) = default;
   //@}
   //**********************************************************************************************

 public:
   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~NegativeInfinity() = default;
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   constexpr operator signed char() const;
   constexpr operator char()        const;
   constexpr operator wchar_t()     const;
   constexpr operator short()       const;
   constexpr operator int()         const;
   constexpr operator long()        const;
#if BLAZE_WIN32_PLATFORM || BLAZE_WIN64_PLATFORM
   constexpr operator ptrdiff_t()   const;
#endif
   constexpr operator float()       const;
   constexpr operator double()      const;
   constexpr operator long double() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename T >
   constexpr bool equal( const T& rhs ) const;
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   NegativeInfinity& operator=( const NegativeInfinity& ) = delete;
   void* operator&() const = delete;
   //@}
   //**********************************************************************************************

 private:
   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   friend class Infinity;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor of the NegativeInfinity class.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::NegativeInfinity()
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to the signed char built-in type.
//
// The conversion operator returns the smallest possible signed char value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator signed char() const
{
   return Limits<signed char>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the char built-in type.
//
// The conversion operator returns the smallest possible char value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator char() const
{
   return Limits<char>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the wchar_t built-in type.
//
// The conversion operator returns the smallest possible wchar_t value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator wchar_t() const
{
   return Limits<wchar_t>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the short built-in type.
//
// The conversion operator returns the smallest possible short value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator short() const
{
   return Limits<short>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the int built-in type.
//
// The conversion operator returns the smallest possible int value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator int() const
{
   return Limits<int>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the long built-in type.
//
// The conversion operator returns the smallest possible long value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator long() const
{
   return Limits<long>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_WIN32_PLATFORM || BLAZE_WIN64_PLATFORM
/*!\brief Conversion operator to the ptrdiff_t built-in type.
//
// The conversion operator returns the smallest possible ptrdiff_t value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator ptrdiff_t() const
{
   return Limits<ptrdiff_t>::ninf();
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the float built-in type.
//
// The conversion operator returns the smallest possible float value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator float() const
{
   return Limits<float>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the double built-in type.
//
// The conversion operator returns the smallest possible double value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator double() const
{
   return Limits<double>::ninf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the long double built-in type.
//
// The conversion operator returns the smallest possible long double value.
*/
template< typename I >  // Positive infinity type
constexpr NegativeInfinity<I>::operator long double() const
{
   return Limits<long double>::ninf();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Equality comparison to a built-in data type.
//
// This function compares built-in data types with their largest possible value. The function
// only works for built-in data types. The attempt to compare user-defined class types will
// result in a compile time error.
*/
template< typename I >  // Positive infinity type
template< typename T >  // Built-in data type
constexpr bool NegativeInfinity<I>::equal( const T& rhs ) const
{
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( T );
   return Limits<T>::ninf() == rhs;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name NegativeInfinity operators */
//@{
template< typename I1, typename I2 >
constexpr bool operator==( const NegativeInfinity<I1>& lhs, const NegativeInfinity<I2>& rhs );

template< typename I, typename T >
constexpr bool operator==( const NegativeInfinity<I>& lhs, const T& rhs );

template< typename I, typename T >
constexpr bool operator==( const T& lhs, const NegativeInfinity<I>& rhs );

template< typename I1, typename I2 >
constexpr bool operator!=( const NegativeInfinity<I1>& lhs, const NegativeInfinity<I2>& rhs );

template< typename I, typename T >
constexpr bool operator!=( const NegativeInfinity<I>& lhs, const T& rhs );

template< typename I, typename T >
constexpr bool operator!=( const T& lhs, const NegativeInfinity<I>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two NegativeInfinity objects.
// \ingroup math
//
// \return \a true.
*/
template< typename I1    // Left-hand side positive infinity type
        , typename I2 >  // Right-hand side positive infinity type
constexpr bool operator==( const NegativeInfinity<I1>& /*lhs*/, const NegativeInfinity<I2>& /*rhs*/ )
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an NegativeInfinity object and a built-in data type.
// \ingroup math
//
// \param lhs The left-hand side NegativeInfinity object.
// \param rhs The right-hand side built-in data value.
// \return \a true if the built-in data value is negative infinity, \a false if not.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename I    // Positive infinity type
        , typename T >  // Built-in data type
constexpr bool operator==( const NegativeInfinity<I>& lhs, const T& rhs )
{
   return lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a built-in data type and an NegativeInfinity object.
// \ingroup math
//
// \param lhs The left-hand side built-in data value.
// \param rhs The right-hand side NegativeInfinity object.
// \return \a true if the built-in data value is negative infinity, \a false if not.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename I    // Positive infinity type
        , typename T >  // Built-in data type
constexpr bool operator==( const T& lhs, const NegativeInfinity<I>& rhs )
{
   return rhs.equal( lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two NegativeInfinity objects.
// \ingroup math
//
// \return \a false.
*/
template< typename I1    // Left-hand side positive infinity type
        , typename I2 >  // Right-hand side positive infinity type
constexpr bool operator!=( const NegativeInfinity<I1>& /*lhs*/, const NegativeInfinity<I2>& /*rhs*/ )
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between an NegativeInfinity object and a built-in data type.
// \ingroup math
//
// \param lhs The left-hand side NegativeInfinity object.
// \param rhs The right-hand side built-in data value.
// \return \a true if the built-in data value is not negative infinity, \a false if it is.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename I    // Positive infinity type
        , typename T >  // Built-in data type
constexpr bool operator!=( const NegativeInfinity<I>& lhs, const T& rhs )
{
   return !lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a built-in data type and an NegativeInfinity object.
// \ingroup math
//
// \param lhs The left-hand side built-in data value.
// \param rhs The right-hand side NegativeInfinity object.
// \return \a true if the built-in data value is not negative infinity, \a false if it is.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename I    // Positive infinity type
        , typename T >  // Built-in data type
constexpr bool operator!=( const T& lhs, const NegativeInfinity<I>& rhs )
{
   return !rhs.equal( lhs );
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Positive infinity for built-in data types.
// \ingroup math
//
// The Infinity class is a wrapper class around the functionality of the blaze::Limits class
// to provide the possiblity to assign a positive infinity value to built-in data types.
// As positive infinity value, the largest possible positive value of the corresponding
// data type is used. In order to assign the positive infinity value, the Infinity class
// can be implicitly converted to the following 13 built-in integral and floating point
// data types:
//
// <ul>
//    <li>integers</li>
//    <ul>
//       <li>unsigned char, signed char, char, wchar_t</li>
//       <li>unsigned short, short</li>
//       <li>unsigned int, int</li>
//       <li>unsigned long, long</li>
//    </ul>
//    <li>floating points</li>
//    <ul>
//       <li>float</li>
//       <li>double</li>
//       <li>long double</li>
//    </ul>
// </ul>
//
// In order to be able to assign infinity values, the global Infinity instance blaze::inf
// is provided, which can be used wherever a built-in data type is required.

   \code
   int i    =  inf;  // Assigns a positive infinity value
   double d = -inf;  // Assigns a negative infinity value
   ...
   \endcode
*/
class Infinity
{
 public:
   //**Type definitions****************************************************************************
   using NegativeType = NegativeInfinity<Infinity>;  //!< The negative infinity type.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr Infinity();
   Infinity( const Infinity& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Infinity() = default;
   //@}
   //**********************************************************************************************

   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   constexpr operator unsigned char()  const;
   constexpr operator signed char()    const;
   constexpr operator char()           const;
   constexpr operator wchar_t()        const;
   constexpr operator unsigned short() const;
   constexpr operator short()          const;
   constexpr operator unsigned int()   const;
   constexpr operator int()            const;
   constexpr operator unsigned long()  const;
   constexpr operator long()           const;
#if BLAZE_WIN32_PLATFORM || BLAZE_WIN64_PLATFORM
   constexpr operator size_t()         const;
   constexpr operator ptrdiff_t()      const;
#endif
   constexpr operator float()          const;
   constexpr operator double()         const;
   constexpr operator long double()    const;
   //@}
   //**********************************************************************************************

   //**Arithmetic operators************************************************************************
   /*!\name Arithmetic operators */
   //@{
   constexpr const Infinity&    operator+() const;
   constexpr const NegativeType operator-() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   template< typename T >
   constexpr bool equal( const T& rhs ) const;
   //@}
   //**********************************************************************************************

   //**Forbidden operations************************************************************************
   /*!\name Forbidden operations */
   //@{
   Infinity& operator=( const Infinity& ) = delete;
   void* operator&() const = delete;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor of the Infinity class.
*/
constexpr Infinity::Infinity()
{}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator to the unsigned char built-in type.
//
// The conversion operator returns the largest possible unsigned char value.
*/
constexpr Infinity::operator unsigned char() const
{
   return Limits<unsigned char>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the char built-in type.
//
// The conversion operator returns the largest possible char value.
*/
constexpr Infinity::operator char() const
{
   return Limits<char>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the signed char built-in type.
//
// The conversion operator returns the largest possible signed char value.
*/
constexpr Infinity::operator signed char() const
{
   return Limits<signed char>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the wchar_t built-in type.
//
// The conversion operator returns the largest possible wchar_t value.
*/
constexpr Infinity::operator wchar_t() const
{
   return Limits<wchar_t>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the unsigned short built-in type.
//
// The conversion operator returns the largest possible unsigned short value.
*/
constexpr Infinity::operator unsigned short() const
{
   return Limits<unsigned short>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the short built-in type.
//
// The conversion operator returns the largest possible short value.
*/
constexpr Infinity::operator short() const
{
   return Limits<short>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the unsigned int built-in type.
//
// The conversion operator returns the largest possible unsigned int value.
*/
constexpr Infinity::operator unsigned int() const
{
   return Limits<unsigned int>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the int built-in type.
//
// The conversion operator returns the largest possible int value.
*/
constexpr Infinity::operator int() const
{
   return Limits<int>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the unsigned long built-in type.
//
// The conversion operator returns the largest possible unsigned long value.
*/
constexpr Infinity::operator unsigned long() const
{
   return Limits<unsigned long>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the long built-in type.
//
// The conversion operator returns the largest possible long value.
*/
constexpr Infinity::operator long() const
{
   return Limits<long>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_WIN32_PLATFORM || BLAZE_WIN64_PLATFORM
/*!\brief Conversion operator to the size_t built-in type.
//
// The conversion operator returns the largest possible size_t value.
*/
constexpr Infinity::operator size_t() const
{
   return Limits<size_t>::inf();
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_WIN32_PLATFORM || BLAZE_WIN64_PLATFORM
/*!\brief Conversion operator to the ptrdiff_t built-in type.
//
// The conversion operator returns the largest possible ptrdiff_t value.
*/
constexpr Infinity::operator ptrdiff_t() const
{
   return Limits<ptrdiff_t>::inf();
}
#endif
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the float built-in type.
//
// The conversion operator returns the largest possible float value.
*/
constexpr Infinity::operator float() const
{
   return Limits<float>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the double built-in type.
//
// The conversion operator returns the largest possible double value.
*/
constexpr Infinity::operator double() const
{
   return Limits<double>::inf();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion operator to the long double built-in type.
//
// The conversion operator returns the largest possible long double value.
*/
constexpr Infinity::operator long double() const
{
   return Limits<long double>::inf();
}
//*************************************************************************************************




//=================================================================================================
//
//  ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the positive infinity value for all built-in data types.
//
// \return The positive infinity value.
*/
constexpr const Infinity& Infinity::operator+() const
{
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the negative infinity value for all built-in data types.
//
// \return The negative infinity value.
*/
constexpr const Infinity::NegativeType Infinity::operator-() const
{
   return NegativeType();
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Equality comparison to a built-in data type.
//
// This function compares built-in data types with their largest possible value. The function
// only works for built-in data types. The attempt to compare user-defined class types will
// result in a compile time error.
*/
template< typename T >
constexpr bool Infinity::equal( const T& rhs ) const
{
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( T );
   return Limits<T>::inf() == rhs;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Infinity operators */
//@{
constexpr bool operator==( const Infinity& lhs, const Infinity& rhs );

template< typename I >
constexpr bool operator==( const Infinity& lhs, const NegativeInfinity<I>& rhs );

template< typename I >
constexpr bool operator==( const NegativeInfinity<I>& lhs, const Infinity& rhs );

template< typename T >
constexpr bool operator==( const Infinity& lhs, const T& rhs );

template< typename T >
constexpr bool operator==( const T& lhs, const Infinity& rhs );

constexpr bool operator!=( const Infinity& lhs, const Infinity& rhs );

template< typename I >
constexpr bool operator!=( const Infinity& lhs, const NegativeInfinity<I>& rhs );

template< typename I >
constexpr bool operator!=( const NegativeInfinity<I>& lhs, const Infinity& rhs );

template< typename T >
constexpr bool operator!=( const Infinity& lhs, const T& rhs );

template< typename T >
constexpr bool operator!=( const T& lhs, const Infinity& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two Infinity objects.
// \ingroup math
//
// \return \a true.
*/
constexpr bool operator==( const Infinity& /*lhs*/, const Infinity& /*rhs*/ )
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an Infinity object and a NegativeInfinity object.
// \ingroup math
//
// \return \a false.
*/
template< typename I >  // Positive infinity type
constexpr bool operator==( const Infinity& /*lhs*/, const NegativeInfinity<I>& /*rhs*/ )
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a NegativeInfinity object and an Infinity object.
// \ingroup math
//
// \return \a false.
*/
template< typename I >  // Positive infinity type
constexpr bool operator==( const NegativeInfinity<I>& /*lhs*/, const Infinity& /*rhs*/ )
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an Infinity object and a built-in data type.
// \ingroup math
//
// \param lhs The left-hand side Infinity object.
// \param rhs The right-hand side built-in data value.
// \return \a true if the built-in data value is infinity, \a false if not.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename T >
constexpr bool operator==( const Infinity& lhs, const T& rhs )
{
   return lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a built-in data type and an Infinity object.
// \ingroup math
//
// \param lhs The left-hand side built-in data value.
// \param rhs The right-hand side Infinity object.
// \return \a true if the built-in data value is infinity, \a false if not.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename T >
constexpr bool operator==( const T& lhs, const Infinity& rhs )
{
   return rhs.equal( lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two Infinity objects.
// \ingroup math
//
// \return \a false.
*/
constexpr bool operator!=( const Infinity& /*lhs*/, const Infinity& /*rhs*/ )
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between an Infinity object and a NegativeInfinity object.
// \ingroup math
//
// \return \a true.
*/
template< typename I >  // Positive infinity type
constexpr bool operator!=( const Infinity& /*lhs*/, const NegativeInfinity<I>& /*rhs*/ )
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a NegativeInfinity object and an Infinity object.
// \ingroup math
//
// \return \a true.
*/
template< typename I >  // Positive infinity type
constexpr bool operator!=( const NegativeInfinity<I>& /*lhs*/, const Infinity& /*rhs*/ )
{
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between an Infinity object and a built-in data type.
// \ingroup math
//
// \param lhs The left-hand side Infinity object.
// \param rhs The right-hand side built-in data value.
// \return \a true if the built-in data value is not infinity, \a false if it is.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename T >
constexpr bool operator!=( const Infinity& lhs, const T& rhs )
{
   return !lhs.equal( rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a built-in data type and an Infinity object.
// \ingroup math
//
// \param lhs The left-hand side built-in data value.
// \param rhs The right-hand side Infinity object.
// \return \a true if the built-in data value is not infinity, \a false if it is.
//
// This operator works only for built-in data types. The attempt to compare user-defined class
// types will result in a compile time error.
*/
template< typename T >
constexpr bool operator!=( const T& lhs, const Infinity& rhs )
{
   return !rhs.equal( lhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL INFINITY VALUE
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global Infinity instance.
// \ingroup math
//
// The blaze::inf instance can be used wherever a built-in data type is expected. It is implicitly
// converted to the corresponding built-in data type and represents its largest possible data
// value.
*/
constexpr Infinity inf;
//*************************************************************************************************

} // namespace blaze

#endif
