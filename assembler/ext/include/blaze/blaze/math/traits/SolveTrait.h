//=================================================================================================
/*!
//  \file blaze/math/traits/SolveTrait.h
//  \brief Header file for the Solver trait
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

#ifndef _BLAZE_MATH_TRAITS_SOLVETRAIT_H_
#define _BLAZE_MATH_TRAITS_SOLVETRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/InvalidType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename, typename, typename = void > struct SolveTrait;
template< typename, typename, typename = void > struct SolveTraitEval1;
template< typename, typename, typename = void > struct SolveTraitEval2;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
auto evalSolveTrait( const volatile T1&, const volatile T2& ) -> SolveTraitEval1<T1,T2>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the SolveTrait class.
// \ingroup math_traits
//
// \section solvetrait_general General
//
// The SolveTrait class template offers the possibility to select the resulting data type for
// solving a linear system of equations (LSE) with the two given types \a T1 and \a T2 via the
// \a solve() operation. SolveTrait defines the nested type \a Type, which represents the
// resulting data type of the \a solve() operation. In case \a T1 and \a T2 cannot be combined
// in a \a solve() operation, there is no nested type \a Type. Note that \a const and \a volatile
// qualifiers and reference modifiers are generally ignored.
//
//
// \n \section solvetrait_specializations Creating custom specializations
//
// Per default, SolveTrait supports all vector and matrix types of the Blaze library (including
// views and adaptors). For all other data types it is possible to specialize the SolveTrait
// template. The following example shows the according specialization for solving an LSE with a
// static matrix and a static vector:

   \code
   template< typename T, size_t N, bool SO, bool TF >
   struct SolveTrait< StaticMatrix<T,N,N,SO>, StaticVector<T,N,TF> >
   {
      using Type = StaticVector<T,N,TF>;
   };
   \endcode

// \n \section solvetrait_examples Examples
//
// The following example demonstrates the use of the SolveTrait template, where depending on
// the two given data types the resulting data type is selected:

   \code
   template< typename T1, typename T2 >    // The two generic types
   typename SolveTrait<T1,T2>::Type        // The resulting generic return type
   solveLSE( const T1& t1, const T2& t2 )  //
   {                                       // The function 'solveLSE' returns the solution
      return solve( t1, t2 );              // for the LSE represented by the two given values
   }                                       //
   \endcode
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SolveTrait
   : public decltype( evalSolveTrait( std::declval<T1&>(), std::declval<T2&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the SolveTrait class template.
// \ingroup math_traits
//
// The SolveTrait_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the SolveTrait class template. For instance, given the types \a T1 and \a T2 the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::SolveTrait<T1,T2>::Type;
   using Type2 = blaze::SolveTrait_t<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
using SolveTrait_t = typename SolveTrait<T1,T2>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief First auxiliary helper struct for the SolveTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SolveTraitEval1
   : public SolveTraitEval2<T1,T2>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Second auxiliary helper struct for the SolveTrait type trait.
// \ingroup math_traits
*/
template< typename T1  // Type of the left-hand side operand
        , typename T2  // Type of the right-hand side operand
        , typename >   // Restricting condition
struct SolveTraitEval2
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
