//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsCommutative.h
//  \brief Header file for the IsCommutative type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISCOMMUTATIVE_H_
#define _BLAZE_MATH_TYPETRAITS_ISCOMMUTATIVE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typetraits/IsNumeric.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for the commutativity of data types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the two given data types \a T1 and T2 are commutative
// with respect to mathematical operations. If the types are commutative, the \a value member
// constant is set to \a true, the nested type definition \a Type is \a TrueType, and the class
// derives from \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType,
// and the class derives from \a FalseType.

   \code
   using VT = StaticVector<int,3UL>;
   using MT = StaticMatrix<int,3UL,3UL>;

   blaze::IsCommutative< double, complex<double> >::value                 // Evaluates to 1
   blaze::IsCommutative< DynamicVector<int>, DynamicVector<int> >::Type   // Results in TrueType
   blaze::IsCommutative< DynamicMatrix<VT>, DynamicMatrix<VT> >           // Is derived from TrueType
   blaze::IsCommutative< DynamicMatrix<int>, DynamicVector<int> >::value  // Evaluates to 0
   blaze::IsCommutative< DynamicVector<MT>, DynamicVector<VT> >::Type     // Results in FalseType
   blaze::IsCommutative< DynamicMatrix<VT>, DynamicMatrix<MT> >           // Is derived from FalseType
   \endcode
*/
template< typename T1, typename T2, typename = void >
struct IsCommutative
   : public BoolConstant< IsNumeric_v<T1> || IsNumeric_v<T2> >
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsCommutative type trait for vectors and matrices.
// \ingroup math_type_traits
*/
template< typename T1, typename T2 >
struct IsCommutative< T1, T2
                    , EnableIf_t< ( IsVector_v<T1> && IsVector_v<T2> ) ||
                                  ( IsMatrix_v<T1> && IsMatrix_v<T2> ) > >
   : public BoolConstant< IsCommutative< UnderlyingElement_t<T1>, UnderlyingElement_t<T2> >::value >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsCommutative type trait.
// \ingroup math_type_traits
//
// The IsCommutative_v variable template provides a convenient shortcut to access the nested
// \a value of the IsCommutative class template. For instance, given the type1 \a T1 and T2
// the following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsCommutative<T1,T2>::value;
   constexpr bool value2 = blaze::IsCommutative_v<T1,T2>;
   \endcode
*/
template< typename T1, typename T2 >
constexpr bool IsCommutative_v = IsCommutative<T1,T2>::value;
//*************************************************************************************************

} // namespace blaze

#endif
