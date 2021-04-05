//=================================================================================================
/*!
//  \file blaze/math/traits/SchurTrait.h
//  \brief Header file for the Schur product trait
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

#ifndef _BLAZE_MATH_TRAITS_SCHURTRAIT_H_
#define _BLAZE_MATH_TRAITS_SCHURTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, typename = void > struct SchurTrait;
template< typename, typename, typename = void > struct SchurTraitEval1;
template< typename, typename, typename = void > struct SchurTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
auto evalSchurTrait( const volatile T1&, const volatile T2& ) -> SchurTraitEval1<T1,T2>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the SchurTrait class.
// \ingroup math_traits
//
// \section schurtrait_general General
//
// The SchurTrait class template offers the possibility to select the resulting data type of
// a generic Schur product operation between the two given types \a T1 and \a T2. SchurTrait
// defines the nested type \a Type, which represents the resulting data type of the Schur
// product. In case \a T1 and \a T2 cannot be combined in a Schur product, there is no nested
// type \a Type. Note that \a const and \a volatile qualifiers and reference modifiers are
// generally ignored.
//
//
// \n \section schurtrait_specializations Creating custom specializations
//
// Per default, SchurTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the SchurTrait template. The
// following example shows the according specialization for the Schur product between two static
// matrices with equal storage order:

   \code
   template< typename T1, size_t M, size_t N, bool SO, typename T2 >
   struct SchurTrait< StaticMatrix<T1,M,N,SO>, StaticMatrix<T2,M,N,SO> >
   {
      using Type = StaticMatrix< MultTrait_t<T1,T2>, M, N, SO >;
   };
   \endcode

// \n \section schurtrait_examples Examples
//
// The following example demonstrates the use of the SchurTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >  // The two generic types
   typename SchurTrait<T1,T2>::Type      // The resulting generic return type
   schur( T1 t1, T2 t2 )                 //
   {                                     // The function 'schur' returns the Schur
      return t1 % t2;                    // product of the two given values
   }                                     //
   \endcode
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SchurTrait
   : decltype( evalSchurTrait( std::declval<T1&>(), std::declval<T2&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the SchurTrait class template.
// \ingroup math_traits
//
// The SchurTrait_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the SchurTrait class template. For instance, given the types \a T1 and \a T2 the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::SchurTrait<T1,T2>::Type;
   using Type2 = blaze::SchurTrait_t<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
using SchurTrait_t = typename SchurTrait<T1,T2>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the SchurTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SchurTraitEval1
   : public SchurTraitEval2<T1,T2>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the SchurTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SchurTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
