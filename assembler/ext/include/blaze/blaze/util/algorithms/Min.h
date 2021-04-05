//=================================================================================================
/*!
//  \file blaze/util/algorithms/Min.h
//  \brief Header file for the generic min algorithm
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

#ifndef _BLAZE_UTIL_ALGORITHMS_MIN_H_
#define _BLAZE_UTIL_ALGORITHMS_MIN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/HasLessThan.h>
#include <blaze/util/typetraits/IsSigned.h>
#include <blaze/util/typetraits/IsUnsigned.h>
#include <blaze/util/typetraits/RemoveConst.h>
#include <blaze/util/typetraits/RemoveCVRef.h>
#include <blaze/util/typetraits/RemoveRValueReference.h>


namespace blaze {

//=================================================================================================
//
//  MIN ALGORITHMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Minimum function for two values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \return The minimum of the two values/objects.
//
// This function determines the minimum of the two given values by means of a less-than comparison.
// The return type of the function is determined by the data types of the given arguments.
*/
template< typename T1, typename T2
        , typename R1 = RemoveCVRef_t<T1>
        , typename R2 = RemoveCVRef_t<T2>
        , EnableIf_t< HasLessThan_v<R2,R1> &&
                      !( IsSigned_v<R1> && IsUnsigned_v<R2> ) &&
                      !( IsUnsigned_v<R1> && IsSigned_v<R2> ) >* = nullptr >
constexpr decltype(auto) min( T1&& a, T2&& b )
{
   using Result = decltype( b < a ? std::forward<T2>( b ) : std::forward<T1>( a ) );
   using Return = RemoveConst_t< RemoveRValueReference_t< Result > >;

   return static_cast<Return>( b < a ? std::forward<T2>( b ) : std::forward<T1>( a ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Minimum function for three values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \param c The third value/object.
// \return The minimum of the given values/objects.
//
// This function returns the minimum of the given values/objects. The return type of the function
// is determined by the data types of the given arguments.
*/
template< typename T1, typename T2, typename T3 >
constexpr decltype(auto) min( T1&& a, T2&& b, T3&& c )
{
   using std::forward;
   using blaze::min;

   return min( min( forward<T1>( a ), forward<T2>( b ) ), forward<T3>( c ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Minimum function for at least four values/objects.
// \ingroup algorithms
//
// \param a The first value/object.
// \param b The second value/object.
// \param c The third value/object.
// \param args The pack of additional values/objects.
// \return The minimum of the given values/objects.
//
// This function returns the minimum of the given values/objects. The return type of the function
// is determined by the data types of the given arguments.
*/
template< typename T1, typename T2, typename T3, typename... Ts >
constexpr decltype(auto) min( T1&& a, T2&& b, T3&& c, Ts&&... args )
{
   using std::forward;
   using blaze::min;

   return min( min( min( forward<T1>( a ), forward<T2>( b ) ), forward<T3>( c ) ), forward<Ts>( args )... );
}
//*************************************************************************************************

} // namespace blaze

#endif
