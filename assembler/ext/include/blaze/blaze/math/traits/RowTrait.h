//=================================================================================================
/*!
//  \file blaze/math/traits/RowTrait.h
//  \brief Header file for the row trait
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

#ifndef _BLAZE_MATH_TRAITS_ROWTRAIT_H_
#define _BLAZE_MATH_TRAITS_ROWTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Infinity.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, size_t... > struct RowTrait;
template< typename, size_t, typename = void > struct RowTraitEval1;
template< typename, size_t, typename = void > struct RowTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t I, typename T >
auto evalRowTrait( const volatile T& ) -> RowTraitEval1<T,I>;

template< typename T >
auto evalRowTrait( const volatile T& ) -> RowTraitEval1<T,inf>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the RowTrait class.
// \ingroup math_traits
//
// \section rowtrait_general General
//
// The RowTrait class template offers the possibility to select the resulting data type when
// creating a view on a specific row of a dense or sparse matrix. In case the given type \a MT
// is a dense or sparse matrix type, RowTrait defines the nested type \a Type, which represents
// the resulting data type of the row operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section rowtrait_specializations Creating custom specializations
//
// Per default, RowTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the RowTrait template. The
// following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, size_t... CRAs >
   struct RowTrait< DynamicMatrix<T1,SO>, CRAs... >
   {
      using Type = DynamicVector<T1,true>;
   };
   \endcode

// \n \section rowtrait_examples Examples
//
// The following example demonstrates the use of the RowTrait template, where depending on
// the given matrix type the resulting row type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the row type of a row-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,rowMajor>;
   using ResultType1 = typename blaze::RowTrait<MatrixType1>::Type;

   // Definition of the row type for the 1st row of a column-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,4UL,3UL,columnMajor>;
   using ResultType2 = typename blaze::RowTrait<MatrixType2,1UL>::Type;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CRAs >  // Compile time row arguments
struct RowTrait
   : public decltype( evalRowTrait<CRAs...>( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RowTrait type trait.
// \ingroup math_traits
//
// The RowTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the RowTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::RowTrait<MT>::Type;
   using Type2 = blaze::RowTrait_t<MT>;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CRAs >  // Compile time row arguments
using RowTrait_t = typename RowTrait<MT,CRAs...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the RowTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t I     // Compile time row index
        , typename >   // Restricting condition
struct RowTraitEval1
   : public RowTraitEval2<MT,I>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the RowTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t I     // Compile time row index
        , typename >   // Restricting condition
struct RowTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
