//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSquare.h
//  \brief Header file for the IsSquare type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSQUARE_H_
#define _BLAZE_MATH_TYPETRAITS_ISSQUARE_H_


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
/*!\brief Compile time check for square matrices.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a square matrix type
// (i.e. a matrix type that is guaranteed to be square at compile time). In case the type is
// a square matrix type, the \a value member constant is set to \a true, the nested type
// definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise
// \a value is set to \a false, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   using blaze::rowMajor;

   // Type definitions of square matrix types
   using Mat2x2 = blaze::StaticMatrix<double,2UL,2UL,rowMajor>;
   using Mat3x3 = blaze::StaticMatrix<double,3UL,3UL,rowMajor>;
   using Mat4x4 = blaze::StaticMatrix<double,4UL,4UL,rowMajor>;

   // Type definitions of non-square matrix types
   using Mat2x3            = blaze::StaticMatrix<double,2UL,3UL,rowMajor>;
   using DynamicMatrixType = blaze::DynamicMatrix<double,rowMajor>;
   using HybridMatrixType  = blaze::HybridMatrix<double,3UL,3UL,rowMajor>;

   blaze::IsSquare< Mat2x2 >::value              // Evaluates to 1
   blaze::IsSquare< const Mat3x3 >::Type         // Results in TrueType
   blaze::IsSquare< volatile Mat4x4 >            // Is derived from TrueType
   blaze::IsSquare< DynamicMatrixType >::value   // Evaluates to 0
   blaze::IsSquare< const Mat2x3 >::Type         // Results in FalseType
   blaze::IsSquare< volatile HybridMatrixType >  // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsSquare
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSquare type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSquare< const T >
   : public IsSquare<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSquare type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSquare< volatile T >
   : public IsSquare<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSquare type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSquare< const volatile T >
   : public IsSquare<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSquare type trait.
// \ingroup math_type_traits
//
// The IsSquare_v variable template provides a convenient shortcut to access the nested \a value
// of the IsSquare class template. For instance, given the type \a T the following two statements
// are identical:

   \code
   constexpr bool value1 = blaze::IsSquare<T>::value;
   constexpr bool value2 = blaze::IsSquare_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsSquare_v = IsSquare<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
