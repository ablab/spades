//=================================================================================================
/*!
//  \file blaze/math/traits/RepeatTrait.h
//  \brief Header file for the repeat trait
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

#ifndef _BLAZE_MATH_TRAITS_REPEATTRAIT_H_
#define _BLAZE_MATH_TRAITS_REPEATTRAIT_H_


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
template< typename, size_t... > struct RepeatTrait;
template< typename, size_t, size_t, size_t, typename = void > struct RepeatTraitEval1;
template< typename, size_t, size_t, size_t, typename = void > struct RepeatTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< size_t R0, size_t R1, size_t R2, typename T >
auto evalRepeatTrait( const volatile T& ) -> RepeatTraitEval1<T,R0,R1,R2>;

template< size_t R0, size_t R1, typename T >
auto evalRepeatTrait( const volatile T& ) -> RepeatTraitEval1<T,R0,R1,inf>;

template< size_t R0, typename T >
auto evalRepeatTrait( const volatile T& ) -> RepeatTraitEval1<T,R0,inf,inf>;

template< typename T >
auto evalRepeatTrait( const volatile T& ) -> RepeatTraitEval1<T,inf,inf,inf>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the RepeatTrait class.
// \ingroup math_traits
//
// \section repeattrait_general General
//
// The RepeatTrait class template offers the possibility to select the resulting data type when
// repeating a dense or sparse vector or matrix. RepeatTrait defines the nested type \a Type,
// which represents the resulting data type of the repeat operation. In case the given data type
// is not a dense or sparse vector or matrix type, there is no nested type \a Type. Note that
// \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section repeattrait_specializations Creating custom specializations
//
// Per default, RepeatTrait supports all vector and matrix types of the Blaze library (including
// views and adaptors). For all other data types it is possible to specialize the RepeatTrait
// template. The following example shows the according specialization for the DynamicVector class
// template:

   \code
   template< typename Type, bool TF, size_t... CRAs >
   struct RepeatTrait< DynamicVector<Type,TF>, CRAs... >
   {
      using Type = DynamicVector<Type,TF>;
   };
   \endcode

// \n \section repeattrait_examples Examples
//
// The following example demonstrates the use of the RepeatTrait template, where depending on
// the given vector or matrix type the resulting type is selected:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   // Definition of the resulting type of a dynamic column vector
   using VectorType1 = blaze::DynamicVector<int,columnVector>;
   using ResultType1 = typename blaze::RepeatTrait<VectorType1>::Type;

   // Definition of the resulting type of a static row vector
   using VectorType2 = blaze::StaticVector<int,5UL,rowVector>;
   using ResultType2 = typename blaze::RepeatTrait<VectorType2>::Type;
   \endcode
*/
template< typename T        // Type of the operand
        , size_t... CRAs >  // Compile time repetition arguments
struct RepeatTrait
   : public decltype( evalRepeatTrait<CRAs...>( std::declval<T&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RepeatTrait type trait.
// \ingroup math_traits
//
// The RepeatTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the RepeatTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::RepeatTrait<MT>::Type;
   using Type2 = blaze::RepeatTrait_t<MT>;
   \endcode
*/
template< typename T        // Type of the operand
        , size_t... CRAs >  // Compile time repeater arguments
using RepeatTrait_t = typename RepeatTrait<T,CRAs...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the RepeatTrait type trait.
// \ingroup math_traits
*/
template< typename T  // Type of the operand
        , size_t R0   // Number of repetitions in the first dimension
        , size_t R1   // Number of repetitions in the second dimension
        , size_t R2   // Number of repetitions in the third dimension
        , typename >  // Restricting condition
struct RepeatTraitEval1
   : public RepeatTraitEval2<T,R0,R1,R2>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the RepeatTrait type trait.
// \ingroup math_traits
*/
template< typename T  // Type of the operand
        , size_t R0   // Number of repetitions in the first dimension
        , size_t R1   // Number of repetitions in the second dimension
        , size_t R2   // Number of repetitions in the third dimension
        , typename >  // Restricting condition
struct RepeatTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
