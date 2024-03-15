//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsScalar.h
//  \brief Header file for the IsScalar type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSCALAR_H_
#define _BLAZE_MATH_TYPETRAITS_ISSCALAR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/HasCompositeType.h>
#include <blaze/math/typetraits/IsProxy.h>
#include <blaze/math/typetraits/IsSIMDPack.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/IsPointer.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for scalar types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is scalar type (i.e. neither
// a vector, nor a matrix, proxy, SIMD pack, reference or pointer type). In case the type is a
// scalar type, the \a value member constant is set to \a true, the nested type definition \a Type
// is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to \a false,
// \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::IsScalar< int >::value                         // Evaluates to 1
   blaze::IsScalar< const double >::Type                 // Results in TrueType
   blaze::IsScalar< volatile complex<float> >            // Is derived from TrueType
   blaze::IsScalar< StaticVector<double,3UL> >::value    // Evaluates to 0
   blaze::IsScalar< const DynamicMatrix<double> >::Type  // Results in FalseType
   blaze::IsScalar< volatile CompressedMatrix<int> >     // Is derived from FalseType
   \endcode
*/

template< typename T >
struct IsScalar
   : public BoolConstant< !HasCompositeType_v<T> &&
                          !IsProxy_v<T> &&
                          !IsSIMDPack_v<T> &&
                          !IsReference_v<T> &&
                          !IsPointer_v<T> >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsScalar type trait.
// \ingroup math_type_traits
//
// The IsScalar_v variable template provides a convenient shortcut to access the nested \a value
// of the IsScalar class template. For instance, given the type \a T the following two statements
// are identical:

   \code
   constexpr bool value1 = blaze::IsScalar<T>::value;
   constexpr bool value2 = blaze::IsScalar_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsScalar_v = IsScalar<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
