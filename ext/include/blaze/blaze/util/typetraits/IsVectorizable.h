//=================================================================================================
/*!
//  \file blaze/util/typetraits/IsVectorizable.h
//  \brief Header file for the IsVectorizable type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ISVECTORIZABLE_H_
#define _BLAZE_UTIL_TYPETRAITS_ISVECTORIZABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Vectorization.h>
#include <blaze/util/Complex.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsVectorizable type trait.
// \ingroup type_traits
*/
template< typename T >
struct IsVectorizableHelper
{
 public:
   //**********************************************************************************************
   static constexpr bool value = ( sizeof(T) <= 8UL ) &&
                                 ( ( bool( BLAZE_SSE_MODE      ) && IsFloat_v<T>   ) ||
                                   ( bool( BLAZE_SSE2_MODE     ) && IsNumeric_v<T> ) ||
                                   ( bool( BLAZE_AVX512BW_MODE ) && IsNumeric_v<T> ) ||
                                   ( bool( BLAZE_AVX512F_MODE || BLAZE_MIC_MODE )
                                     && IsNumeric_v<T> && sizeof(T) >= 4UL ) );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsVectorizableHelper class template for 'complex'.
// \ingroup type_traits
*/
template< typename T >
struct IsVectorizableHelper< complex<T> >
{
 public:
   //**********************************************************************************************
   static constexpr bool value = IsVectorizableHelper< RemoveCV_t<T> >::value;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsVectorizableHelper class template for 'void'.
// \ingroup type_traits
*/
template<>
struct IsVectorizableHelper<void>
{
 public:
   //**********************************************************************************************
   static constexpr bool value = false;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for vectorizable types.
// \ingroup type_traits
//
// Depending on the available instruction set (SSE, SSE2, SSE3, SSE4, AVX, AVX2, MIC, ...),
// this type trait tests whether or not the given template parameter is a vectorizable type,
// i.e. a type for which intrinsic vector operations and optimizations can be used. Currently,
// all built-in data types except \c bool and the according complex numbers are considered to
// be vectorizable types. In case the type is vectorizable, the \a value member constant is
// set to \a true, the nested type definition \a Type is \a TrueType, and the class derives
// from \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType, and the
// class derives from \a FalseType.

   \code
   blaze::IsVectorizable< int >::value         // Evaluates to 'true'
   blaze::IsVectorizable< const float >::Type  // Results in TrueType
   blaze::IsVectorizable< volatile double >    // Is derived from TrueType
   blaze::IsVectorizable< void >::value        // Evaluates to 'false'
   blaze::IsVectorizable< const bool >::Type   // Results in FalseType
   blaze::IsVectorizable< volatile MyClass >   // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsVectorizable
   : public BoolConstant< IsVectorizableHelper< RemoveCV_t<T> >::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsVectorizable type trait.
// \ingroup type_traits
//
// The IsVectorizable_v variable template provides a convenient shortcut to access the nested
// \a value of the IsVectorizable class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsVectorizable<T>::value;
   constexpr bool value2 = blaze::IsVectorizable_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsVectorizable_v = IsVectorizable<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
