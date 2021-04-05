//=================================================================================================
/*!
//  \file blaze/math/traits/BandTrait.h
//  \brief Header file for the band trait
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

#ifndef _BLAZE_MATH_TRAITS_BANDTRAIT_H_
#define _BLAZE_MATH_TRAITS_BANDTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Infinity.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, ptrdiff_t... > struct BandTrait;
template< typename, ptrdiff_t, typename = void > struct BandTraitEval1;
template< typename, ptrdiff_t, typename = void > struct BandTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< ptrdiff_t I, typename T >
auto evalBandTrait( const volatile T& ) -> BandTraitEval1<T,I>;

template< typename T >
auto evalBandTrait( const volatile T& ) -> BandTraitEval1<T,inf>;
   /*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the BandTrait class.
// \ingroup math_traits
//
// \section bandtrait_general General
//
// The BandTrait class template offers the possibility to select the resulting data type when
// creating a view on a specific band of a dense or sparse matrix. In case the given type \a MT
// is a dense or sparse matrix type, BandTrait defines the nested type \a Type, which represents
// the resulting data type of the band operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section bandtrait_specializations Creating custom specializations
//
// Per default, BandTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the BandTrait template. The
// following example shows the according specialization for the DynamicMatrix class template:

   \code
   template< typename T1, bool SO, ptrdiff_t... CBAs >
   struct BandTrait< DynamicMatrix<T1,SO>, CBAs... >
   {
      using Type = DynamicVector<T1,true>;
   };
   \endcode

// \n \section bandtrait_examples Examples
//
// The following example demonstrates the use of the BandTrait template, where depending on
// the given matrix type the resulting vector type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the fitting type for any band of a row-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,rowMajor>;
   using ResultType1 = typename blaze::BandTrait<MatrixType1>::Type;

   // Definition of the fitting type of the 3rd band of a column-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,3UL,3UL,columnMajor>;
   using ResultType2 = typename blaze::BandTrait<MatrixType2,3L>::Type;
   \endcode
*/
template< typename MT          // Type of the matrix
        , ptrdiff_t... CBAs >  // Compile time band arguments
struct BandTrait
   : public decltype( evalBandTrait<CBAs...>( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the BandTrait type trait.
// \ingroup math_traits
//
// The BandTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the BandTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::BandTrait<MT>::Type;
   using Type2 = blaze::BandTrait_t<MT>;
   \endcode
*/
template< typename MT          // Type of the matrix
        , ptrdiff_t... CBAs >  // Compile time band arguments
using BandTrait_t = typename BandTrait<MT,CBAs...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the BandTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , ptrdiff_t I  // Compile time band index
        , typename >   // Restricting condition
struct BandTraitEval1
   : public BandTraitEval2<MT,I>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the BandTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , ptrdiff_t I  // Compile time band index
        , typename >   // Restricting condition
struct BandTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
