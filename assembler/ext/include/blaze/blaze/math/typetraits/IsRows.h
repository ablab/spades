//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsRows.h
//  \brief Header file for the IsRows type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISROWS_H_
#define _BLAZE_MATH_TYPETRAITS_ISROWS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/views/Forward.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for row selections.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a row selection (i.e. a
// view on rows of a dense or sparse matrix). In case the type is a row selection, the \a value
// member constant is set to \a true, the nested type definition \a Type is \a TrueType, and the
// class derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType,
// and the class derives from \a FalseType.

   \code
   using blaze::aligned;

   using MatrixType1 = blaze::StaticMatrix<int,10UL,16UL>;
   using MatrixType2 = blaze::DynamicMatrix<double>;
   using MatrixType3 = blaze::CompressedMatrix<float>;

   MatrixType1 A;
   MatrixType2 B( 100UL, 200UL );
   MatrixType3 C( 200UL, 250UL );

   using RowsType1 = decltype( blaze::rows<2UL,4UL>( A ) );
   using RowsType2 = decltype( blaze::rows( B, 8UL, 24UL ) );
   using RowsType3 = decltype( blaze::rows( C, 5UL, 13UL ) );

   blaze::IsRows< RowsType1 >::value         // Evaluates to 1
   blaze::IsRows< const RowsType2 >::Type    // Results in TrueType
   blaze::IsRows< volatile RowsType3 >       // Is derived from TrueType
   blaze::IsRows< MatrixType1 >::value       // Evaluates to 0
   blaze::IsRows< const MatrixType2 >::Type  // Results in FalseType
   blaze::IsRows< volatile MatrixType3 >     // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsRows
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRows type trait for 'Rows'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, typename... CRAs >
struct IsRows< Rows<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRows type trait for 'const Rows'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, typename... CRAs >
struct IsRows< const Rows<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRows type trait for 'volatile Rows'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, typename... CRAs >
struct IsRows< volatile Rows<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRows type trait for 'const volatile Rows'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, typename... CRAs >
struct IsRows< const volatile Rows<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsRows type trait.
// \ingroup math_type_traits
//
// The IsRows_v variable template provides a convenient shortcut to access the nested \a value
// of the IsRows class template. For instance, given the type \a T the following two statements
// are identical:

   \code
   constexpr bool value1 = blaze::IsRows<T>::value;
   constexpr bool value2 = blaze::IsRows_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsRows_v = IsRows<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
