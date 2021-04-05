//=================================================================================================
/*!
//  \file blaze/util/mpl/Not.h
//  \brief Header file for the Not_t alias template
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

#ifndef _BLAZE_UTIL_MPL_NOT_H_
#define _BLAZE_UTIL_MPL_NOT_H_


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
/*!\brief Compile time type negation.
// \ingroup mpl
//
// The Not_t alias template negates the given compile time condition. In case the given condition
// would evaluate to \a true, the nested member enumeration is set to \a false and vice versa:

   \code
   using namespace blaze;

   Not_t< IsIntegral<int> >::value    // Evaluates to false
   Not_t< IsDouble<int>   >::value    // Evaluates to true
   Not_t< IsSigned<int> >::ValueType  // Results in bool
   \endcode
*/
template< typename C >  // Condition to be negated
using Not_t = Bool_t< !C::value >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the Not_t alias.
// \ingroup mpl
//
// The Not_v variable template provides a convenient shortcut to access the nested \a value of
// the Not_t alias. For instance, given the type \a C the following two statements are identical:

   \code
   constexpr bool value1 = Not_t<C>::value;
   constexpr bool value2 = Not_v<C>;
   \endcode
*/
template< typename C >  // Condition to be negated
constexpr bool Not_v = Not_t<C>::value;
//*************************************************************************************************

} // namespace blaze

#endif
