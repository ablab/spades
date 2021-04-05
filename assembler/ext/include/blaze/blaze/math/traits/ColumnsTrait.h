//=================================================================================================
/*!
//  \file blaze/math/traits/ColumnsTrait.h
//  \brief Header file for the columns trait
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

#ifndef _BLAZE_MATH_TRAITS_COLUMNSTRAIT_H_
#define _BLAZE_MATH_TRAITS_COLUMNSTRAIT_H_


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
template< typename, size_t > struct ColumnsTrait;
template< typename, size_t, typename = void > struct ColumnsTraitEval1;
template< typename, size_t, typename = void > struct ColumnsTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t N, typename T >
auto evalColumnsTrait( const volatile T& ) -> ColumnsTraitEval1<T,N>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the ColumnsTrait class.
// \ingroup math_traits
//
// \section columnstrait_general General
//
// The ColumnsTrait class template offers the possibility to select the resulting data type when
// creating a view on a set of columns of a dense or sparse matrix. In case the given type \a MT
// is a dense or sparse matrix type, ColumnsTrait defines the nested type \a Type, which represents
// the resulting data type of the columns operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section columnstrait_specializations Creating custom specializations
//
// Per default, ColumnsTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the ColumnsTrait template.
// The following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, size_t N >
   struct ColumnsTrait< DynamicMatrix<T1,SO>, N >
   {
      using Type = DynamicMatrix<T1,true>;
   };
   \endcode

// \n \section columnstrait_examples Examples
//
// The following example demonstrates the use of the ColumnsTrait template, where depending on
// the given matrix type the resulting columns type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the columns type of a column-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,columnMajor>;
   using ResultType1 = typename blaze::ColumnsTrait<MatrixType1,0UL>::Type;

   // Definition of the columns type for two columns of a row-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,4UL,3UL,rowMajor>;
   using ResultType2 = typename blaze::ColumnsTrait<MatrixType2,2UL>::Type;
   \endcode
*/
template< typename MT  // Type of the matrix
        , size_t N >   // Number of compile time indices
struct ColumnsTrait
   : public decltype( evalColumnsTrait<N>( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the ColumnsTrait type trait.
// \ingroup math_traits
//
// The ColumnsTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the ColumnsTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::ColumnsTrait<MT>::Type;
   using Type2 = blaze::ColumnsTrait_t<MT>;
   \endcode
*/
template< typename MT  // Type of the matrix
        , size_t N >   // Number of compile time indices
using ColumnsTrait_t = typename ColumnsTrait<MT,N>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the ColumnsTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t N     // Number of compile time indices
        , typename >   // Restricting condition
struct ColumnsTraitEval1
   : public ColumnsTraitEval2<MT,N>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the ColumnsTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t N     // Number of compile time indices
        , typename >   // Restricting condition
struct ColumnsTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
