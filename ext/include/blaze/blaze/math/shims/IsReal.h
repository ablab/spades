//=================================================================================================
/*!
//  \file blaze/math/shims/IsReal.h
//  \brief Header file for the isReal shim
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

#ifndef _BLAZE_MATH_SHIMS_ISREAL_H_
#define _BLAZE_MATH_SHIMS_ISREAL_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/system/Inline.h>
#include <blaze/util/Complex.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/typetraits/IsBuiltin.h>


namespace blaze {

//=================================================================================================
//
//  ISREAL SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the given value/object represents a real number.
// \ingroup math_shims
//
// \param v The value to be tested.
// \return \a true in case the given value represents a real number, \a false otherwise.
//
// The \a isReal shim provides an abstract interface for testing a value/object of any type
// whether it represents the a real number. In case the value/object is of built-in type, the
// function returns \a true:

   \code
   int    i = 1;    // isReal( i ) returns true
   double d = 1.0;  // isReal( d ) returns true
   \endcode

// In case the value/object is of complex type, the function returns \a true if the imaginary
// part is close to 0 (relaxed semantics):

   \code
   complex<double> c1( 1.0, 0.0 );  // isReal( c1 ) returns true
   complex<double> c2( 0.0, 1.0 );  // isReal( c2 ) returns false
   \endcode

// Optionally, it is possible to switch between relaxed semantics (blaze::relaxed) and strict
// semantics (blaze::strict). In case of strict semantics, for complex types the function returns
// \a true in case the imaginary part is exactly zero:

   \code
                                     // isReal<strict>( ... ) | isReal<relaxed>( ... )
   complex<double> c1( 1.0, 0.0  );  //    true               |    true
   complex<double> c2( 1.0, 1E-9 );  //    false              |    true (below 1E-8)
   complex<double> c3( 0.0, 1.0  );  //    false              |    false
   \endcode

// For all other types the function returns \a false.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type >    // Type of the given value
BLAZE_ALWAYS_INLINE bool isReal( const Type& v ) noexcept
{
   MAYBE_UNUSED( v );

   return IsBuiltin_v<Type>;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Overload of the \a isReal function for complex data types.
// \ingroup math_shims
//
// \param v The complex number to be tested.
// \return \a true in case the imaginary part is equal to 0, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type >    // Type of the given value
BLAZE_ALWAYS_INLINE bool isReal( const complex<Type>& v ) noexcept( IsBuiltin_v<Type> )
{
   return IsBuiltin_v<Type> && isZero<RF>( v.imag() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given value/object represents a real number.
// \ingroup math_shims
//
// \param v The value to be tested.
// \return \a true in case the given value represents a real number, \a false otherwise.
//
// The \a isReal shim provides an abstract interface for testing a value/object of any type
// whether it represents the a real number. In case the value/object is of built-in type, the
// function returns \a true. In case the value/object is of complex type, the function returns
// \a true if the imaginary part is equal to 0. Otherwise it returns \a false.

   \code
   int    i = 1;                      // isReal( i ) returns true
   double d = 1.0;                    // isReal( d ) returns true

   complex<double> c1( 1.0, 0.0 );    // isReal( c1 ) returns true
   complex<double> c2( 0.0, 1.0 );    // isReal( c2 ) returns false

   blaze::DynamicVector<int> vec;     // isReal( vec ) returns false
   blaze::DynamicMatrix<double> mat;  // isReal( mat ) returns false
   \endcode
*/
template< typename Type >
BLAZE_ALWAYS_INLINE bool isReal( const Type& v ) noexcept
{
   return isReal<relaxed>( v );
}
//*************************************************************************************************

} // namespace blaze

#endif
