//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsUniform.h
//  \brief Header file for the IsUniform type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISUNIFORM_H_
#define _BLAZE_MATH_TYPETRAITS_ISUNIFORM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T > struct IsUniform;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsUniform type trait.
// \ingroup math_traits
*/
template< typename T
        , typename = void >
struct IsUniformHelper
   : public IsZero<T>
{};

template< typename T >  // Type of the operand
struct IsUniformHelper< T, EnableIf_t< IsExpression_v<T> && !IsSame_v<T,typename T::ResultType> > >
   : public IsUniform< typename T::ResultType >::Type
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for uniform vectors and matrices.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a uniform vector or
// matrix type (i.e. a data type that is guaranteed to be a uniform vector or matrix at compile
// time). In case the type is a uniform vector or matrix type, the \a value member constant is
// set to \a true, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class
// derives from \a FalseType.

   \code
   blaze::IsUniform< UniformVector<int> >::value          // Evaluates to 1
   blaze::IsUniform< const UniformMatrix<float> >::Type   // Results in TrueType
   blaze::IsUniform< volatile UniformVector<double> >     // Is derived from TrueType
   blaze::IsUniform< DynamicVector<int> >::value          // Evaluates to 0
   blaze::IsUniform< const DynamicMatrix<float> >::Type   // Results in FalseType
   blaze::IsUniform< volatile CompressedVector<double> >  // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsUniform
   : public IsUniformHelper<T>
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsUniform type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsUniform< const T >
   : public IsUniform<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsUniform type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsUniform< volatile T >
   : public IsUniform<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsUniform type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsUniform< const volatile T >
   : public IsUniform<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsUniform type trait.
// \ingroup math_type_traits
//
// The IsUniform_v variable template provides a convenient shortcut to access the nested
// \a value of the IsUniform class template. For instance, given the type \a T the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::IsUniform<T>::value;
   constexpr bool value2 = blaze::IsUniform_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsUniform_v = IsUniform<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
