//=================================================================================================
/*!
//  \file blaze/math/shims/Digits.h
//  \brief Header file for the digits shim
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

#ifndef _BLAZE_MATH_SHIMS_DIGITS_H_
#define _BLAZE_MATH_SHIMS_DIGITS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/constraints/Integral.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  DIGITS SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the number of valid digits of an integral value.
// \ingroup math_shims
//
// \param a The integral value.
// \return The number of valid digits.
//
// This function counts the number of valid digits in the given integral value.

   \code
   digits( 100   );  // Returns 3
   digits( 12345 );  // Returns 5
   digits( 0     );  // Returns 0
   \endcode

// The digits function only works for integral built-in data types. The attempt to use any
// other type will result in a compile time error.
*/
template< typename T >
constexpr size_t digits( T a ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_INTEGRAL_TYPE( T );

   size_t count( 0 );

   while( a != 0 ) {
      a /= 10;
      ++count;
   }

   return count;
}
//*************************************************************************************************

} // namespace blaze

#endif
