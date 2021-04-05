//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSMPAssignable.h
//  \brief Header file for the IsSMPAssignable type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSMPASSIGNABLE_H_
#define _BLAZE_MATH_TYPETRAITS_ISSMPASSIGNABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/Void.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the IsSMPAssignable type trait.
// \ingroup math_type_traits
*/
template< typename T, typename = void >
struct IsSMPAssignableHelper
   : public FalseType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T >
struct IsSMPAssignableHelper< T, Void_t< decltype( T::smpAssignable ) > >
   : public BoolConstant< T::smpAssignable >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for data types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is an SMP-assignable data
// type (i.e. if it is a data type that can possibly and efficiently be assigned by several
// threads). In this context, built-in data types as well as complex numbers are non-SMP-assignable,
// whereas several vector and matrix types (as for instance DynamicVector and DynamicMatrix) can be
// SMP-assignable. If the type is SMP-assignable, the \a value member constant is set to \a true,
// the nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class derives from
// \a FalseType.

   \code
   using blaze::StaticVector;
   using blaze::StaticMatrix;
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;

   using VectorType = DynamicVector<int,columnVector>;

   VectorType a( 100UL );

   using SubvectorType = decltype( blaze::subvector( a, 8UL, 16UL ) );

   blaze::IsSMPAssignable< VectorType >::value            // Evaluates to 1
   blaze::IsSMPAssignable< SubvectorType >::Type          // Results in TrueType
   blaze::IsSMPAssignable< DynamicMatrix<int> >           // Is derived from TrueType
   blaze::IsSMPAssignable< int >::value                   // Evaluates to 0
   blaze::IsSMPAssignable< StaticVector<int,3UL> >::Type  // Results in FalseType
   blaze::IsSMPAssignable< StaticMatrix<int,4UL,5UL> >    // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsSMPAssignable
   : public BoolConstant< IsSMPAssignableHelper<T>::value >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSMPAssignable type trait.
// \ingroup math_type_traits
//
// The IsSMPAssignable_v variable template provides a convenient shortcut to access the nested
// \a value of the IsSMPAssignable class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsSMPAssignable<T>::value;
   constexpr bool value2 = blaze::IsSMPAssignable_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsSMPAssignable_v = IsSMPAssignable<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
