//=================================================================================================
/*!
//  \file blaze/math/traits/ExpandTrait.h
//  \brief Header file for the expand trait
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

#ifndef _BLAZE_MATH_TRAITS_EXPANDTRAIT_H_
#define _BLAZE_MATH_TRAITS_EXPANDTRAIT_H_


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
template< typename, size_t... > struct ExpandTrait;
template< typename, size_t, typename = void > struct ExpandTraitEval1;
template< typename, size_t, typename = void > struct ExpandTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t E, typename T >
auto evalExpandTrait( const volatile T& ) -> ExpandTraitEval1<T,E>;

template< typename T >
auto evalExpandTrait( const volatile T& ) -> ExpandTraitEval1<T,inf>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the ExpandTrait class.
// \ingroup math_traits
//
// \section expandtrait_general General
//
// The ExpandTrait class template offers the possibility to select the resulting data type when
// expanding a dense or sparse vector or matrix. ExpandTrait defines the nested type \a Type,
// which represents the resulting data type of the expand operation. In case the given data type
// is not a dense or sparse vector or matrix type, there is no nested type \a Type. Note that
// \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section expandtrait_specializations Creating custom specializations
//
// Per default, ExpandTrait supports all vector and matrix types of the Blaze library (including
// views and adaptors). For all other data types it is possible to specialize the ExpandTrait
// template. The following example shows the according specialization for the DynamicVector class
// template:

   \code
   template< typename Type, bool TF, size_t... CEAs >
   struct ExpandTrait< DynamicVector<Type,TF>, CEAs... >
   {
      using Type = DynamicMatrix<Type,( TF == columnVector ? columnMajor : rowMajor )>;
   };
   \endcode

// \n \section expandtrait_examples Examples
//
// The following example demonstrates the use of the ExpandTrait template, where depending on
// the given vector or matrix type the resulting type is selected:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Definition of the resulting type of a dynamic column vector
   using VectorType1 = blaze::DynamicVector<int,columnVector>;
   using ResultType1 = typename blaze::ExpandTrait<VectorType1>::Type;

   // Definition of the resulting type of a static row vector
   using VectorType2 = blaze::StaticVector<int,5UL,rowVector>;
   using ResultType2 = typename blaze::ExpandTrait<VectorType2>::Type;
   \endcode
*/
template< typename T        // Type of the operand
        , size_t... CEAs >  // Compile time expansion arguments
struct ExpandTrait
   : public decltype( evalExpandTrait<CEAs...>( std::declval<T&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the ExpandTrait type trait.
// \ingroup math_traits
//
// The ExpandTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the ExpandTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::ExpandTrait<MT>::Type;
   using Type2 = blaze::ExpandTrait_t<MT>;
   \endcode
*/
template< typename T        // Type of the operand
        , size_t... CEAs >  // Compile time expansion arguments
using ExpandTrait_t = typename ExpandTrait<T,CEAs...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the ExpandTrait type trait.
// \ingroup math_traits
*/
template< typename T  // Type of the operand
        , size_t E    // Compile time expansion
        , typename >  // Restricting condition
struct ExpandTraitEval1
   : public ExpandTraitEval2<T,E>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the ExpandTrait type trait.
// \ingroup math_traits
*/
template< typename T  // Type of the operand
        , size_t E    // Compile time expansion
        , typename >  // Restricting condition
struct ExpandTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
