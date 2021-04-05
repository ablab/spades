//=================================================================================================
/*!
//  \file blaze/util/mpl/Greater.h
//  \brief Header file for the Greater_t alias template
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

#ifndef _BLAZE_UTIL_MPL_GREATER_H_
#define _BLAZE_UTIL_MPL_GREATER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time type comparison.
// \ingroup mpl
//
// The Greater_t alias template compares the two given types using the greater-than operator
// ('>'). In case \a T1::value is greater than \a T2::value, the nested \a value member is set
// to 1. Otherwise it is set to 0.

   \code
   using namespace blaze;

   Greater_t< Int_t<5> , Int_t<2>  >::value      // Evaluates to true
   Greater_t< Int_t<5> , Long_t<2> >::value      // Evaluates to true
   Greater_t< Long_t<2>, Int_t<2>  >::value      // Evaluates to false
   Greater_t< Int_t<2> , Long_t<5> >::value      // Evaluates to false
   Greater_t< Int_t<5> , Int_t<2>  >::ValueType  // Results in bool
   \endcode
*/
template< typename T1    // The type of the left-hand side operand
        , typename T2 >  // The type of the right-hand side operand
using Greater_t = Bool_t< ( T1::value > T2::value ) >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the Greater_t alias.
// \ingroup mpl
//
// The Greater_v variable template provides a convenient shortcut to access the nested \a value of
// the Greater_t alias. For instance, given the types \a T1 and \a T2 the following two statements
// are identical:

   \code
   constexpr bool value1 = Greater_t<T1,T2>::value;
   constexpr bool value2 = Greater_v<T1,T2>;
   \endcode
*/
template< typename T1    // The type of the left-hand side operand
        , typename T2 >  // The type of the right-hand side operand
constexpr bool Greater_v = Greater_t<T1,T2>::value;
//*************************************************************************************************

} // namespace blaze

#endif
