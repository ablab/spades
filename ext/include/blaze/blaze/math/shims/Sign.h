//=================================================================================================
/*!
//  \file blaze/math/shims/Sign.h
//  \brief Header file for the sign shim
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

#ifndef _BLAZE_MATH_SHIMS_SIGN_H_
#define _BLAZE_MATH_SHIMS_SIGN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/EnableIf.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsSigned.h>


namespace blaze {

//=================================================================================================
//
//  SIGN SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluating the sign of the given value.
// \ingroup math_shims
//
// \param a The given value.
// \return 1 if the value is greater than zero, 0 if it is zero, and -1 if it is less than zero.
//
// The sign function evaluates the sign of the given value \a a of the built-in data type \a T.
// It returns 1 if \a a is greater than zero, 0 if \a a is zero, and -1 if \a a is less than zero.
*/
template< typename T, typename = EnableIf_t< IsIntegral_v<T> > >
constexpr T sign( T a ) noexcept
{
   return ( IsSigned_v<T> )
          ? ( T(0) < a ) - ( a < T(0) )
          : ( T(0) < a );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluating the sign of the given single precision value.
// \ingroup math_shims
//
// \param a The given single precision value.
// \return 1 if the value is greater than zero, 0 if it is zero, and -1 if it is less than zero.
//
// The sign function evaluates the sign of the given single precision value \a a. It returns 1.0
// if \a a is greater than zero, 0.0 if \a a is zero, -1.0 if \a a is less than zero, and \c NaN
// if a is \c NaN.
*/
constexpr float sign( float a ) noexcept
{
   if     ( 0.0F < a ) return  1.0F;
   else if( a < 0.0F ) return -1.0F;
   else                return     a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluating the sign of the given double precision value.
// \ingroup math_shims
//
// \param a The given double precision value.
// \return 1 if the value is greater than zero, 0 if it is zero, and -1 if it is less than zero.
//
// The sign function evaluates the sign of the given double precision value \a a. It returns 1.0
// if \a a is greater than zero, 0.0 if \a a is zero, -1.0 if \a a is less than zero, and \c NaN
// if a is \c NaN.
*/
constexpr double sign( double a ) noexcept
{
   if     ( 0.0 < a ) return  1.0;
   else if( a < 0.0 ) return -1.0;
   else               return    a;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluating the sign of the given extended precision value.
// \ingroup math_shims
//
// \param a The given extended precision value.
// \return 1 if the value is greater than zero, 0 if it is zero, and -1 if it is less than zero.
//
// The sign function evaluates the sign of the given extended precision value \a a. It returns 1.0
// if \a a is greater than zero, 0.0 if \a a is zero, -1.0 if \a a is less than zero, and \c NaN
// if a is \c NaN.
*/
constexpr long double sign( long double a ) noexcept
{
   if     ( 0.0L < a ) return  1.0L;
   else if( a < 0.0L ) return -1.0L;
   else                return     a;
}
//*************************************************************************************************

} // namespace blaze

#endif
