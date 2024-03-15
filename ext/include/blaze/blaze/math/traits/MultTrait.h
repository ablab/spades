//=================================================================================================
/*!
//  \file blaze/math/traits/MultTrait.h
//  \brief Header file for the multiplication trait
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

#ifndef _BLAZE_MATH_TRAITS_MULTTRAIT_H_
#define _BLAZE_MATH_TRAITS_MULTTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/GroupTag.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/CommonType.h>
#include <blaze/util/typetraits/IsNumeric.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blaze/util/typetraits/Void.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, typename = void > struct MultTrait;
template< typename, typename, typename = void > struct MultTraitEval1;
template< typename, typename, typename = void > struct MultTraitEval2;
template< typename, typename, typename = void > struct MultTraitEval3;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
auto evalMultTrait( const volatile T1&, const volatile T2& ) -> MultTraitEval1<T1,T2>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the MultTrait class.
// \ingroup math_traits
//
// \section multtrait_general General
//
// The MultTrait class template offers the possibility to select the resulting data type of
// a generic multiplication operation between the two given types \a T1 and \a T2. MultTrait
// defines the nested type \a Type, which represents the resulting data type of the multiplication.
// In case the two types \a T1 and \a T2 cannot be multiplied, there is no nested type \a Type.
// Note that \c const and \c volatile qualifiers and reference modifiers are generally ignored.
//
//
// \n \section multtrait_specializations Creating custom specializations
//
// MultTrait is guaranteed to work for all built-in data types, complex numbers, all vector
// and matrix types of the Blaze library (including views and adaptors) and all data types that
// provide a multiplication operator (i.e. \c operator*). In order to add support for user-defined
// data types that either don't provide a multiplication operator or whose multiplication operator
// returns a proxy object instead of a concrete type (as it is common in expression template
// libraries) it is possible to specialize the MultTrait template. The following example shows
// the according specialization for the multiplication between two dynamic column vectors:

   \code
   template< typename T1, typename T2 >
   struct MultTrait< DynamicVector<T1,columnVector>, DynamicVector<T2,columnVector> >
   {
      using Type = DynamicVector< typename MultTrait<T1,T2>::Type, columnVector >;
   };
   \endcode

// \n \section multtrait_examples Examples
//
// The following example demonstrates the use of the MultTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >  // The two generic types
   typename MultTrait<T1,T2>::Type       // The resulting generic return type
   mult( const T1& t1, const T2& t2 )    //
   {                                     // The function 'mult' returns the
      return t1 * t2;                    // product of the two given values
   }                                     //
   \endcode
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct MultTrait
   : public decltype( evalMultTrait( std::declval<T1&>(), std::declval<T2&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTrait class template for two numeric types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct MultTrait< T1, T2, EnableIf_t< IsNumeric_v<T1> && IsNumeric_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = CommonType_t<T1,T2>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the MultTrait class template.
// \ingroup math_traits
//
// The MultTrait_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the MultTrait class template. For instance, given the types \a T1 and \a T2 the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::MultTrait<T1,T2>::Type;
   using Type2 = blaze::MultTrait_t<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
using MultTrait_t = typename MultTrait<T1,T2>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the MultTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct MultTraitEval1
   : public MultTraitEval2<T1,T2>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTraitEval1 class template for two 'GroupTag'.
// \ingroup math_traits
*/
template< size_t ID >
struct MultTraitEval1<GroupTag<ID>,GroupTag<ID>,void>
{
 public:
   //**********************************************************************************************
   using Type = GroupTag<ID>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTraitEval1 class template for a 'GroupTag' and a numeric type.
// \ingroup math_traits
*/
template< size_t ID, typename T2 >
struct MultTraitEval1< GroupTag<ID>, T2 >
{
 public:
   //**********************************************************************************************
   using Type = GroupTag<ID>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTraitEval1 class template for a numeric type and a'GroupTag'.
// \ingroup math_traits
*/
template< typename T1, size_t ID >
struct MultTraitEval1< T1, GroupTag<ID> >
{
 public:
   //**********************************************************************************************
   using Type = GroupTag<ID>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the MultTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct MultTraitEval2
   : public MultTraitEval3<T1,T2>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTraitEval2 class template the inner product operation.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsRowVector_v<T1> &&
                                   IsColumnVector_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = MultTrait_t< typename T1::ElementType, typename T2::ElementType >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Third auxiliary helper struct for the MultTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct MultTraitEval3
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the MultTraitEval3 class template for two types supporting multiplication.
// \ingroup math_traits
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
struct MultTraitEval3< T1, T2
                     , Void_t< decltype( std::declval<T1>() * std::declval<T2>() ) > >
{
 public:
   //**********************************************************************************************
   using Type = RemoveCVRef_t< decltype( std::declval<T1>() * std::declval<T2>() ) >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
