//=================================================================================================
/*!
//  \file blaze/math/traits/RowsTrait.h
//  \brief Header file for the rows trait
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

#ifndef _BLAZE_MATH_TRAITS_ROWSTRAIT_H_
#define _BLAZE_MATH_TRAITS_ROWSTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, size_t > struct RowsTrait;
template< typename, size_t, typename = void > struct RowsTraitEval1;
template< typename, size_t, typename = void > struct RowsTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t M , typename T >
auto evalRowsTrait( const volatile T& ) -> RowsTraitEval1<T,M>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the RowsTrait class.
// \ingroup math_traits
//
// \section rowstrait_general General
//
// The RowsTrait class template offers the possibility to select the resulting data type when
// creating a view on a set of rows of a dense or sparse matrix. In case the given type \a MT
// is a dense or sparse matrix type, RowsTrait defines the nested type \a Type, which represents
// the resulting data type of the rows operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section rowstrait_specializations Creating custom specializations
//
// Per default, RowsTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the RowsTrait template. The
// following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, size_t M >
   struct RowsTrait< DynamicMatrix<T1,SO>, M >
   {
      using Type = DynamicMatrix<T1,false>;
   };
   \endcode

// \n \section rowstrait_examples Examples
//
// The following example demonstrates the use of the RowsTrait template, where depending on
// the given matrix type the resulting rows type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the rows type of a row-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,rowMajor>;
   using ResultType1 = typename blaze::RowsTrait<MatrixType1,0UL>::Type;

   // Definition of the rows type for two rows of a column-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,4UL,3UL,columnMajor>;
   using ResultType2 = typename blaze::RowsTrait<MatrixType2,2UL>::Type;
   \endcode
*/
template< typename MT  // Type of the matrix
        , size_t M >   // Number of compile time indices
struct RowsTrait
   : public decltype( evalRowsTrait<M>( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RowsTrait type trait.
// \ingroup math_traits
//
// The RowsTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the RowsTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::RowsTrait<MT>::Type;
   using Type2 = blaze::RowsTrait_t<MT>;
   \endcode
*/
template< typename MT  // Type of the matrix
        , size_t M >   // Number of compile time indices
using RowsTrait_t = typename RowsTrait<MT,M>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the RowsTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t M     // Number of compile time indices
        , typename >   // Restricting condition
struct RowsTraitEval1
   : public RowsTraitEval2<MT,M>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the RowsTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t M     // Number of compile time indices
        , typename >   // Restricting condition
struct RowsTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
