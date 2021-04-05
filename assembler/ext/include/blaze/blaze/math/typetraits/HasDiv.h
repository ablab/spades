//=================================================================================================
/*!
//  \file blaze/math/typetraits/HasDiv.h
//  \brief Header file for the HasDiv type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_HASDIV_H_
#define _BLAZE_MATH_TYPETRAITS_HASDIV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/EnableIf.h>
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
/*!\brief Auxiliary helper struct for the IsDivHelper type trait.
// \ingroup math_type_traits
*/
template< typename T1, typename T2, typename = void >
struct HasDivHelper
   : public FalseType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the HasDivHelper type trait for types providing a division operator.
// \ingroup math_type_traits
*/
template< typename T1, typename T2 >
struct HasDivHelper< T1, T2, Void_t< decltype( std::declval<T1>() / std::declval<T2>() ) > >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Availability of a division operator for the given data types.
// \ingroup math_type_traits
//
// This type trait provides the information whether a division operator (i.e. operator/) exists
// for the two given data types \a T1 and \a T2 (taking the cv-qualifiers into account). In case
// the operator is available, the \a value member constant is set to \a true, the nested type
// definition \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value
// is set to \a false, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   blaze::HasDiv< int, int >::value                         // Evaluates to 1
   blaze::HasDiv< complex<float>, complex<float> >::Type    // Results in TrueType
   blaze::HasDiv< DynamicVector<int>, DynamicVector<int> >  // Is derived from TrueType
   blaze::HasDiv< int, complex<float> >::value              // Evaluates to 0
   blaze::HasDiv< complex<int>, complex<float> >::Type      // Results in FalseType
   blaze::HasDiv< DynamicMatrix<int>, DynamicVector<int> >  // Is derived from FalseType
   \endcode
*/
template< typename T1, typename T2, typename = void >
struct HasDiv
   : public HasDivHelper<T1,T2>
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the HasDiv type trait for vector/vector divisions.
// \ingroup math_type_traits
*/
template< typename T1, typename T2 >
struct HasDiv< T1, T2, EnableIf_t< IsVector_v<T1> && IsDenseVector_v<T2> > >
   : public HasDiv< typename T1::ElementType, typename T2::ElementType >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the HasDiv type trait.
// \ingroup math_type_traits
//
// The HasDiv_v variable template provides a convenient shortcut to access the nested \a value
// of the HasDiv class template. For instance, given the types \a T1 and \a T2 the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::HasDiv<T1,T2>::value;
   constexpr bool value2 = blaze::HasDiv_v<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
constexpr bool HasDiv_v = HasDiv<T1,T2>::value;
//*************************************************************************************************

} // namespace blaze

#endif
