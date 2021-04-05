//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsRow.h
//  \brief Header file for the IsRow type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISROW_H_
#define _BLAZE_MATH_TYPETRAITS_ISROW_H_


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
/*!\brief Compile time check for rows.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a row (i.e. a view on a
// row of a dense or sparse matrix). In case the type is a row, the \a value member constant is
// set to \a true, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class
// derives from \a FalseType.

   \code
   using blaze::aligned;

   using MatrixType1 = blaze::StaticMatrix<int,10UL,16UL>;
   using MatrixType2 = blaze::DynamicMatrix<double>;
   using MatrixType3 = blaze::CompressedMatrix<float>;

   MatrixType1 A;
   MatrixType2 B( 100UL, 200UL );
   MatrixType3 C( 200UL, 250UL );

   using RowType1 = decltype( blaze::row<4UL>( A ) );
   using RowType2 = decltype( blaze::row( B, 16UL ) );
   using RowType3 = decltype( blaze::row( C, 17UL ) );

   blaze::IsRow< RowType1 >::value          // Evaluates to 1
   blaze::IsRow< const RowType2 >::Type     // Results in TrueType
   blaze::IsRow< volatile RowType3 >        // Is derived from TrueType
   blaze::IsRow< MatrixType1 >::value       // Evaluates to 0
   blaze::IsRow< const MatrixType2 >::Type  // Results in FalseType
   blaze::IsRow< volatile MatrixType3 >     // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsRow
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRow type trait for 'Row'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRow< Row<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRow type trait for 'const Row'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRow< const Row<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRow type trait for 'volatile Row'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRow< volatile Row<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsRow type trait for 'const volatile Row'.
// \ingroup math_type_traits
*/
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRow< const volatile Row<MT,SO,DF,SF,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsRow type trait.
// \ingroup math_type_traits
//
// The IsRow_v variable template provides a convenient shortcut to access the nested \a value
// of the IsRow class template. For instance, given the type \a T the following two statements
// are identical:

   \code
   constexpr bool value1 = blaze::IsRow<T>::value;
   constexpr bool value2 = blaze::IsRow_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsRow_v = IsRow<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
