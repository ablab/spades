//=================================================================================================
/*!
//  \file blaze/math/shims/IsDefault.h
//  \brief Header file for the isDefault shim
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

#ifndef _BLAZE_MATH_SHIMS_ISDEFAULT_H_
#define _BLAZE_MATH_SHIMS_ISDEFAULT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <blaze/math/Accuracy.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsBuiltin.h>


namespace blaze {

//=================================================================================================
//
//  ISDEFAULT SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the given value/object is in default state.
// \ingroup math_shims
//
// \param v The value/object to be tested for its default state.
// \return \a true in case the given value/object is in its default state, \a false otherwise.
//
// The \a isDefault shim represents an abstract interface for testing a value/object whether
// it is in its default state or not. In case the value/object is in its default state, the
// function returns \a true, otherwise it returns \a false. For integral built-in data types,
// the function returns \a true in case the current value is zero:

   \code
   int i1 = 0;  // isDefault( i1 ) returns true
   int i2 = 1;  // isDefault( i2 ) returns false
   \endcode

// For floating point built-in data types, the function by default uses relaxed semantics and
// returns \a true in case the current value is close to zero within a certain accuracy:

   \code
   double d1 = 0.0;   // isDefault( d1 ) returns true
   double d2 = 1E-9;  // isDefault( d2 ) returns true since d2 is below 1E-8
   double d3 = 1.0;   // isDefault( d3 ) returns false
   \endcode

// Optionally, it is possible to switch between relaxed semantics (blaze::relaxed) and strict
// semantics (blaze::strict). In case of strict semantics, for floating point built-in data types
// the function returns \a true in case the current value is exactly zero:

   \code
                      // isDefault<strict>( ... ) | isDefault<relaxed>( ... )
   double d1 = 0.0;   //    true                  |    true
   double d2 = 1E-9;  //    false (not 0.0)       |    true (below 1E-8)
   double d3 = 1.0;   //    false                 |    false
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type >    // Type of the given value/object
BLAZE_ALWAYS_INLINE bool isDefault( const Type& v ) noexcept( IsBuiltin_v<Type> )
{
   return v == Type();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given single precision floating point value is zero.
// \ingroup math_shims
//
// \param v The single precision floating point value to be tested for zero.
// \return \a true in case the given value is zero, \a false otherwise.
//
// This overload of the \a isDefault shim tests whether the given single precision floating point
// value is exactly zero or within an epsilon range to zero. In case the value is zero or close
// to zero the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool isDefault( float v ) noexcept
{
   if( RF == relaxed )
      return std::fabs( v ) <= accuracy;
   else
      return v == 0.0F;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given double precision floating point value is zero.
// \ingroup math_shims
//
// \param v The double precision floating point value to be tested for zero.
// \return \a true in case the given value is zero, \a false otherwise.
//
// This overload of the \a isDefault shim tests whether the given double precision floating point
// value is exactly zero or within an epsilon range to zero. In case the value is zero or close
// to zero the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool isDefault( double v ) noexcept
{
   if( RF == relaxed )
      return std::fabs( v ) <= accuracy;
   else
      return v == 0.0;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given extended precision floating point value is zero.
// \ingroup math_shims
//
// \param v The extended precision floating point value to be tested for zero.
// \return \a true in case the given value is zero, \a false otherwise.
//
// This overload of the \a isDefault shim tests whether the given extended precision floating
// point value is exactly zero or within an epsilon range to zero. In case the value is zero or
// close to zero the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF >  // Relaxation flag
BLAZE_ALWAYS_INLINE bool isDefault( long double v ) noexcept
{
   if( RF == relaxed )
      return std::fabs( v ) <= accuracy;
   else
      return v == 0.0L;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given complex number is zero.
// \ingroup math_shims
//
// \param v The complex number to be tested for zero.
// \return \a true in case the given value is zero, \a false otherwise.
//
// This overload of the \a isDefault shim tests whether both the real and the imaginary part of
// the given complex number are exactly zero or within an epsilon range to zero. In case the both
// parts are zero or close to zero the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename T >       // Value type of the complex number
BLAZE_ALWAYS_INLINE bool isDefault( const complex<T>& v ) noexcept( IsBuiltin_v<T> )
{
   return isDefault<RF>( real( v ) ) && isDefault<RF>( imag( v ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given value/object is in default state.
// \ingroup math_shims
//
// \param v The value/object to be tested for its default state.
// \return \a true in case the given value/object is in its default state, \a false otherwise.
//
// The \a isDefault shim represents an abstract interface for testing a value/object whether
// it is in its default state or not. In case the value/object is in its default state, the
// function returns \a true, otherwise it returns \a false. For integral built-in data types,
// the function returns \a true in case the current value is zero:

   \code
   int i1 = 0;  // isDefault( i1 ) returns true
   int i2 = 1;  // isDefault( i2 ) returns false
   \endcode

// For floating point built-in data types, the function by default uses relaxed semantics and
// returns \a true in case the current value is close to zero within a certain accuracy:

   \code
   double d1 = 0.0;   // isDefault( d1 ) returns true
   double d2 = 1E-9;  // isDefault( d2 ) returns true since d2 is below 1E-8
   double d3 = 1.0;   // isDefault( d3 ) returns false
   \endcode

// Optionally, it is possible to switch between relaxed semantics (blaze::relaxed) and strict
// semantics (blaze::strict). In case of strict semantics, for floating point built-in data types
// the function returns \a true in case the current value is exactly zero:

   \code
                      // isDefault<strict>( ... ) | isDefault<relaxed>( ... )
   double d1 = 0.0;   //    true                  |    true
   double d2 = 1E-9;  //    false (not 0.0)       |    true (below 1E-8)
   double d3 = 1.0;   //    false                 |    false
   \endcode
*/
template< typename Type >  // Type of the given value/object
BLAZE_ALWAYS_INLINE bool isDefault( const Type& v ) noexcept( IsBuiltin_v<Type> )
{
   return isDefault<relaxed>( v );
}
//*************************************************************************************************

} // namespace blaze

#endif
