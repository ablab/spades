//=================================================================================================
/*!
//  \file blaze/math/traits/ReduceTrait.h
//  \brief Header file for the reduce trait
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

#ifndef _BLAZE_MATH_TRAITS_REDUCETRAIT_H_
#define _BLAZE_MATH_TRAITS_REDUCETRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, ReductionFlag... > struct ReduceTrait;
template< typename, typename, typename = void > struct TotalReduceTraitEval1;
template< typename, typename, typename = void > struct TotalReduceTraitEval2;
template< typename, typename, ReductionFlag, typename = void > struct PartialReduceTraitEval1;
template< typename, typename, ReductionFlag, typename = void > struct PartialReduceTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< ReductionFlag RF, typename T, typename OP >
auto evalReduceTrait( const volatile T&, OP ) -> PartialReduceTraitEval1<T,OP,RF>;

template< typename T, typename OP >
auto evalReduceTrait( const volatile T&, OP ) -> TotalReduceTraitEval1<T,OP>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the ReduceTrait class.
// \ingroup math_traits
//
// \section reducetrait_general General
//
// The ReduceTrait class template offers the possibility to select the resulting data type of
// a generic reduction operation on the given type \a T. ReduceTrait defines the nested type
// \a Type, which represents the resulting data type of the reduction operation. In case no
// result type can be determined for the type \a T, there is no nested type \a Type. Note that
// \c const and \c volatile qualifiers and reference modifiers are generally ignored.
//
//
// \n \section reducetrait_specializations Creating custom specializations
//
// ReduceTrait is guaranteed to work for all vector and matrix types of the Blaze library
// (including views and adaptors) and all data types that work in combination with the provided
// reduction operation \a OP. In order to add support for user-defined data types or in order to
// adapt to special cases it is possible to specialize the ReduceTrait template. The following
// examples shows the according specialization for total and partial reduction operations with
// a dynamic matrix, respectively:

   \code
   template< typename T, bool SO, typename OP >
   struct ReduceTrait< DynamicMatrix<T,SO>, OP >
   {
      using Type = T;
   };
   \endcode

   \code
   template< typename T, bool SO, typename OP >
   struct ReduceTrait< DynamicMatrix<T,SO>, OP, columnwise >
   {
      using Type = DynamicVector<T,rowVector>;
   };

   template< typename T, bool SO, typename OP >
   struct ReduceTrait< DynamicMatrix<T,SO>, OP, rowwise >
   {
      using Type = DynamicVector<T,columnVector>;
   };
   \endcode
*/
template< typename T             // Type of the operand
        , typename OP            // Type of the reduction operation
        , ReductionFlag... RF >  // Reduction flag
struct ReduceTrait
   : public decltype( evalReduceTrait<RF...>( std::declval<T&>(), std::declval<OP>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the ReduceTrait class template.
// \ingroup math_traits
//
// The ReduceTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the ReduceTrait class template. For instance, given the type \a T and the custom
// operation type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename blaze::ReduceTrait<T,OP>::Type;
   using Type2 = blaze::ReduceTrait_t<T,OP>;
   \endcode
*/
template< typename T             // Type of the operand
        , typename OP            // Type of the reduction operation
        , ReductionFlag... RF >  // Reduction flag
using ReduceTrait_t = typename ReduceTrait<T,OP,RF...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the ReduceTrait type trait.
// \ingroup math_traits
*/
template< typename T   // Type of the operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct TotalReduceTraitEval1
   : public TotalReduceTraitEval2<T,OP>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the ReduceTrait type trait.
// \ingroup math_traits
*/
template< typename T   // Type of the operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct TotalReduceTraitEval2
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the TotalReduceTraitEval2 class template for vectors and matrices.
// \ingroup math_traits
*/
template< typename T     // Type of the operand
        , typename OP >  // Type of the custom operation
struct TotalReduceTraitEval2< T, OP, EnableIf_t< IsVector_v<T> || IsMatrix_v<T> > >
{
 public:
   //**********************************************************************************************
   using Type = ElementType_t<T>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the ReduceTrait type trait.
// \ingroup math_traits
*/
template< typename T        // Type of the operand
        , typename OP       // Type of the custom operation
        , ReductionFlag RF  // Reduction flag
        , typename >        // Restricting condition
struct PartialReduceTraitEval1
   : public PartialReduceTraitEval2<T,OP,RF>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the ReduceTrait type trait.
// \ingroup math_traits
*/
template< typename T        // Type of the operand
        , typename OP       // Type of the custom operation
        , ReductionFlag RF  // Reduction flag
        , typename >        // Restricting condition
struct PartialReduceTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
