//=================================================================================================
/*!
//  \file blaze/util/mpl/Nor.h
//  \brief Header file for the Nor_t alias template
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

#ifndef _BLAZE_UTIL_MPL_NOR_H_
#define _BLAZE_UTIL_MPL_NOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/Bools.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time logical 'not or' evaluation.
// \ingroup mpl
//
// The Nor_t alias template performs at compile time a logical 'not or' evaluation of at least
// two compile time conditions:

   \code
   using namespace blaze;

   using Type = int;

   Nor_t< IsFloat<Type>   , IsDouble<Type>        >::value  // Evaluates to 1
   Nor_t< IsIntegral<Type>, IsSigned<Type>        >::value  // Evaluates to 0
   Nor_t< IsIntegral<Type>, IsFloatingPoint<Type> >::value  // Evaluates to 0
   \endcode
*/
template< typename T1       // Type of the first mandatory operand
        , typename T2       // Type of the second mandatory operand
        , typename... Ts >  // Types of the optional operands
using Nor_t =
   Bool_t< IsSame< Bools< false, T1::value, T2::value, (Ts::value)... >
                 , Bools< T1::value, T2::value, (Ts::value)..., false > >::value >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the Nor_t alias.
// \ingroup mpl
//
// The Nor_v variable template provides a convenient shortcut to access the nested \a value of
// the Nor_t alias. For instance, given the types \a T1 and \a T2 the following two statements
// are identical:

   \code
   constexpr bool value1 = Nor_t<T1,T2>::value;
   constexpr bool value2 = Nor_v<T1,T2>;
   \endcode
*/
template< typename T1       // Type of the first mandatory operand
        , typename T2       // Type of the second mandatory operand
        , typename... Ts >  // Types of the optional operands
constexpr bool Nor_v = Nor_t<T1,T2,Ts...>::value;
//*************************************************************************************************

} // namespace blaze

#endif
