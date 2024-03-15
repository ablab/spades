//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsInitializer.h
//  \brief Header file for the IsInitializer type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISINITIALIZER_H_
#define _BLAZE_MATH_TYPETRAITS_ISINITIALIZER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for custom data types.
// \ingroup math_type_traits
//
// This type trait tests whether the given data type represents an initializer list, i.e. is an
// initializer vector or an initializer matrix. In case the data type represents an initializer
// list, the \a value member constant is set to \a true, the nested type definition \a Type is
// \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set to \a false,
// \a Type is \a FalseType, and the class derives from \a FalseType. Examples:

   \code
   using blaze::InitializerVector;
   using blaze::DynamicVector;
   using blaze::InitializerMatrix;
   using blaze::DynamicMatrix;
   using blaze::columnVector;
   using blaze::rowMajor;

   blaze::IsInitializer< InitializerVector<int,columnVector> >::value     // Evaluates to 1
   blaze::IsInitializer< const InitializerVector<int,rowVector> >::Type   // Results in TrueType
   blaze::IsInitializer< volatile InitializerMatrix<double,rowMajor> >    // Is derived from TrueType
   blaze::IsInitializer< int >::value                                     // Evaluates to 0
   blaze::IsInitializer< const DynamicVector<float,columnVector> >::Type  // Results in FalseType
   blaze::IsInitializer< volatile DynamicMatrix<int,rowMajor> >           // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsInitializer
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsInitializer type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsInitializer< const T >
   : public IsInitializer<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsInitializer type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsInitializer< volatile T >
   : public IsInitializer<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsInitializer type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsInitializer< const volatile T >
   : public IsInitializer<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsInitializer type trait.
// \ingroup math_type_traits
//
// The IsInitializer_v variable template provides a convenient shortcut to access the nested
// \a value of the IsInitializer class template. For instance, given the type \a T the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::IsInitializer<T>::value;
   constexpr bool value2 = blaze::IsInitializer_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsInitializer_v = IsInitializer<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
