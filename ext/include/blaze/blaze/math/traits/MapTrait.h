//=================================================================================================
/*!
//  \file blaze/math/traits/MapTrait.h
//  \brief Header file for the map trait
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

#ifndef _BLAZE_MATH_TRAITS_MAPTRAIT_H_
#define _BLAZE_MATH_TRAITS_MAPTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/GroupTag.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename... > struct MapTrait;
template< typename, typename, typename = void > struct UnaryMapTraitEval1;
template< typename, typename, typename = void > struct UnaryMapTraitEval2;
template< typename, typename, typename = void > struct UnaryMapTraitEval3;
template< typename, typename, typename, typename = void > struct BinaryMapTraitEval1;
template< typename, typename, typename, typename = void > struct BinaryMapTraitEval2;
template< typename, typename, typename, typename = void > struct BinaryMapTraitEval3;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
auto evalMapTrait( const volatile T&, OP ) -> UnaryMapTraitEval1<T,OP>;

template< typename T1, typename T2, typename OP >
auto evalMapTrait( const volatile T1&, const volatile T2&, OP ) -> BinaryMapTraitEval1<T1,T2,OP>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the MapTrait class.
// \ingroup math_traits
//
// \section maptrait_general General
//
// The MapTrait class template offers the possibility to select the resulting data type of a
// generic unary or binary map operation. MapTrait defines the nested type \a Type, which
// represents the resulting data type of the map operation. In case no result type can be
// determined for the type \a T, there is no nested type \a Type. Note that \c const and
// \c volatile qualifiers and reference modifiers are generally ignored.
//
//
// \n \section maptrait_specializations Creating custom specializations
//
// MapTrait is guaranteed to work for all built-in data types, complex numbers, all vector
// and matrix types of the Blaze library (including views and adaptors) and all data types that
// work in combination with the provided custom operation \a OP. In order to add support for
// user-defined data types or in order to adapt to special cases it is possible to specialize
// the MapTrait template. The following example shows the according specialization for map
// operations with a dynamic column vector:

   \code
   template< typename T, typename OP >
   struct MapTrait< DynamicVector<T,columnVector>, OP >
   {
      using Type = DynamicVector< typename MapTrait<T,OP>::Type, columnVector >;
   };
   \endcode
*/
template< typename... Args >  // Types of the map template paramters
struct MapTrait
   : public decltype( evalMapTrait( std::declval<Args&>()... ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the MapTrait class template.
// \ingroup math_traits
//
// The MapTrait_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the MapTrait class template. For instance, given the type \a T and the custom operation
// type \a OP the following two type definitions are identical:

   \code
   using Type1 = typename blaze::MapTrait<T,OP>::Type;
   using Type2 = blaze::MapTrait_t<T,OP>;
   \endcode
*/
template< typename... Args >  // Types of the map template paramters
using MapTrait_t = typename MapTrait<Args...>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T   // Type of the operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct UnaryMapTraitEval1
   : public UnaryMapTraitEval2<T,OP>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UnaryMapTraitEval1 class template for 'GroupTag'.
// \ingroup math_traits
*/
template< size_t ID, typename OP >
struct UnaryMapTraitEval1<GroupTag<ID>,OP,void>
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
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T   // Type of the operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct UnaryMapTraitEval2
   : public UnaryMapTraitEval3<T,OP>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T   // Type of the operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct UnaryMapTraitEval3
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of UnaryMapTraitEval3 for a type supporting the given operation.
// \ingroup math_traits
*/
template< typename T     // Type of the operand
        , typename OP >  // Type of the custom operation
struct UnaryMapTraitEval3< T, OP
                         , Void_t< decltype( std::declval<OP>()( std::declval<T>() ) ) > >
{
 public:
   //**********************************************************************************************
   using Type = RemoveCVRef_t< decltype( std::declval<OP>()( std::declval<T>() ) ) >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTraitEval1
   : public BinaryMapTraitEval2<T1,T2,OP>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the UnaryMapTraitEval1 class template for two 'GroupTag'.
// \ingroup math_traits
*/
template< size_t ID, typename OP >
struct BinaryMapTraitEval1<GroupTag<ID>,GroupTag<ID>,OP,void>
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
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTraitEval2
   : public BinaryMapTraitEval3<T1,T2,OP>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the MapTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename OP  // Type of the custom operation
        , typename >   // Restricting condition
struct BinaryMapTraitEval3
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of BinaryMapTraitEval3 for a type supporting the given operation.
// \ingroup math_traits
*/
template< typename T1    // Type of the left-hand side operand
        , typename T2    // Type of the right-hand side operand
        , typename OP >  // Type of the custom operation
struct BinaryMapTraitEval3< T1, T2, OP
                          , Void_t< decltype( std::declval<OP>()( std::declval<T1>(), std::declval<T2>() ) ) > >
{
 public:
   //**********************************************************************************************
   using Type = RemoveCVRef_t< decltype( std::declval<OP>()( std::declval<T1>(), std::declval<T2>() ) ) >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
