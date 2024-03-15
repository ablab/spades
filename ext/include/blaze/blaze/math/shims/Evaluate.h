//=================================================================================================
/*!
//  \file blaze/math/shims/Evaluate.h
//  \brief Header file for the evaluate shim
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

#ifndef _BLAZE_MATH_SHIMS_EVALUATE_H_
#define _BLAZE_MATH_SHIMS_EVALUATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/HasCompositeType.h>
#include <blaze/math/typetraits/IsProxy.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  EVALUATE SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Formal evaluation of the given argument.
// \ingroup math_shims
//
// \param a The value/object to be evaluated.
// \return The evaluated value/object.
//
// The \a evaluate shim represents an abstract interface for enforcing the evaluation of a given
// value/object of any data type and deducing the correct result type of an operation. For data
// types that don't require an evaluation, as for instance built-in or complex data types, the
// default behavior is not changed.
*/
template< typename T
        , EnableIf_t< !HasCompositeType_v<T> && !IsProxy_v<T> >*  = nullptr >
constexpr T evaluate( const T& a ) noexcept
{
   return a;
}
//*************************************************************************************************

} // namespace blaze

#endif
