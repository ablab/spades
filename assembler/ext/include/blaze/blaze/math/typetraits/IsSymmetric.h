//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSymmetric.h
//  \brief Header file for the IsSymmetric type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSYMMETRIC_H_
#define _BLAZE_MATH_TYPETRAITS_ISSYMMETRIC_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsExpression.h>
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
template< typename T > struct IsSymmetric;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsSymmetric type trait.
// \ingroup math_traits
*/
template< typename T
        , typename = void >
struct IsSymmetricHelper
   : public FalseType
{};

template< typename T >  // Type of the operand
struct IsSymmetricHelper< T, EnableIf_t< IsExpression_v<T> && !IsSame_v<T,typename T::ResultType> > >
   : public IsSymmetric< typename T::ResultType >::Type
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for symmetric matrices.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a symmetric matrix type
// (i.e. a matrix type that is guaranteed to be symmetric at compile time). In case the type is
// a symmetric matrix type, the \a value member constant is set to \a true, the nested type
// definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value
// is set to \a false, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   using blaze::rowMajor;

   using StaticMatrixType     = blaze::StaticMatrix<double,rowMajor>;
   using DynamicMatrixType    = blaze::DynamicMatrix<float,rowMajor>;
   using CompressedMatrixType = blaze::CompressedMatrix<int,rowMajor>;

   using SymmetricStaticType     = blaze::SymmetricMatrix<StaticMatrixType>;
   using SymmetricDynamicType    = blaze::SymmetricMatrix<DynamicMatrixType>;
   using SymmetricCompressedType = blaze::SymmetricMatrix<CompressedMatrixType>;

   blaze::IsSymmetric< SymmetricStaticType >::value        // Evaluates to 1
   blaze::IsSymmetric< const SymmetricDynamicType >::Type  // Results in TrueType
   blaze::IsSymmetric< volatile SymmetricCompressedType >  // Is derived from TrueType
   blaze::IsSymmetric< StaticMatrixType >::value           // Evaluates to 0
   blaze::IsSymmetric< const DynamicMatrixType >::Type     // Results in FalseType
   blaze::IsSymmetric< volatile CompressedMatrixType >     // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsSymmetric
   : public IsSymmetricHelper<T>
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSymmetric type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSymmetric< const T >
   : public IsSymmetric<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSymmetric type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSymmetric< volatile T >
   : public IsSymmetric<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSymmetric type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSymmetric< const volatile T >
   : public IsSymmetric<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSymmetric type trait.
// \ingroup math_type_traits
//
// The IsSymmetric_v variable template provides a convenient shortcut to access the nested
// \a value of the IsSymmetric class template. For instance, given the type \a T the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::IsSymmetric<T>::value;
   constexpr bool value2 = blaze::IsSymmetric_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsSymmetric_v = IsSymmetric<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
