//=================================================================================================
/*!
//  \file blaze/math/traits/SubmatrixTrait.h
//  \brief Header file for the submatrix trait
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

#ifndef _BLAZE_MATH_TRAITS_SUBMATRIXTRAIT_H_
#define _BLAZE_MATH_TRAITS_SUBMATRIXTRAIT_H_


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
template< typename, size_t... > struct SubmatrixTrait;
template< typename, size_t, size_t, size_t, size_t, typename = void > struct SubmatrixTraitEval1;
template< typename, size_t, size_t, size_t, size_t, typename = void > struct SubmatrixTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t I, size_t J, size_t M, size_t N, typename T >
auto evalSubmatrixTrait( const volatile T& ) -> SubmatrixTraitEval1<T,I,J,M,N>;

template< typename T >
auto evalSubmatrixTrait( const volatile T& ) -> SubmatrixTraitEval1<T,inf,inf,inf,inf>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the SubmatrixTrait class.
// \ingroup math_traits
//
// \section submatrixtrait_general General
//
// The SubmatrixTrait class template offers the possibility to select the resulting data type
// when creating a submatrix of a dense or sparse matrix. In case the given type \a MT is a
// dense or sparse matrix type, SubmatrixTrait defines the nested type \a Type, which represents
// the resulting data type of the submatrix operation. Otherwise therer is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section submatrixtrait_specializations Creating custom specializations
//
// Per default, SubmatrixTrait supports all matrix types of the Blaze library (including views
// and adaptors). For all other data types it is possible to specialize the SubmatrixTrait
// template. The following example shows the according specialization for the DynamicMatrix
// class template:

   \code
   template< typename T1, bool SO >
   struct SubmatrixTrait< DynamicMatrix<T1,SO> >
   {
      using Type = DynamicMatrix<T1,SO>;
   };
   \endcode

// \n \section submatrixtrait_examples Examples
//
// The following example demonstrates the use of the SubmatrixTrait template, where depending
// on the given matrix type the according result type is selected:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the result type of a row-major dynamic matrix
   using MatrixType1 = blaze::DynamicMatrix<int,rowMajor>;
   using ResultType1 = typename blaze::SubmatrixTrait<MatrixType1>::Type;

   // Definition of the result type for the inner four elements of a 4x4 column-major static matrix
   using MatrixType2 = blaze::StaticMatrix<int,4UL,4UL,columnMajor>;
   using ResultType2 = typename blaze::SubmatrixTrait<MatrixType2,1UL,1UL,2UL,2UL>::Type;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CSAs >  // Compile time submatrix arguments
struct SubmatrixTrait
   : public decltype( evalSubmatrixTrait<CSAs...>( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the SubmatrixTrait type trait.
// \ingroup math_traits
//
// The SubmatrixTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the SubmatrixTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::SubmatrixTrait<MT>::Type;
   using Type2 = blaze::SubmatrixTrait_t<MT>;
   \endcode
*/
template< typename MT       // Type of the matrix
        , size_t... CSAs >  // Compile time submatrix arguments
using SubmatrixTrait_t = typename SubmatrixTrait<MT,CSAs...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the SubmatrixTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N     // Number of columns
        , typename >   // Restricting condition
struct SubmatrixTraitEval1
   : public SubmatrixTraitEval2<MT,I,J,M,N>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the SubmatrixTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , size_t I     // Index of the first row
        , size_t J     // Index of the first column
        , size_t M     // Number of rows
        , size_t N     // Number of columns
        , typename >   // Restricting condition
struct SubmatrixTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
