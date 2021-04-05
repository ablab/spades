//=================================================================================================
/*!
//  \file blaze/math/shims/Clamp.h
//  \brief Header file for the clamp shim
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

#ifndef _BLAZE_MATH_SHIMS_CLAMP_H_
#define _BLAZE_MATH_SHIMS_CLAMP_H_


namespace blaze {

//=================================================================================================
//
//  CLAMP SHIM
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Restricts the given value to the range \f$[min..max]\f$.
// \ingroup math_shims
//
// \param a The given value/object.
// \param min The lower limit of the range.
// \param max The upper limit of the range.
// \return The clamped value.
//
// The \a clamp shim represents an abstract interface for restricting a value to the specified
// range \f$[min..max]\f$:

   \code
   double d1 =  0.5;
   double d2 =  1.5;
   double d3 = -1.5;

   clamp( d1, -1.0, 1.0 );  // Results in 0.5
   clamp( d2, -1.0, 1.0 );  // Results in 1.0
   clamp( d3, -1.0, 1.0 );  // Results in -1.0
   \endcode
*/
template< typename T >
inline const T& clamp( const T& a, const T& min, const T& max ) noexcept
{
   if( a < min )
      return min;
   else if( max < a )
      return max;
   else return a;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
