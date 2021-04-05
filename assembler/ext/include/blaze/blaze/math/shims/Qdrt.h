//=================================================================================================
/*!
//  \file blaze/math/shims/Qdrt.h
//  \brief Header file for the qdrt shim
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

#ifndef _BLAZE_MATH_SHIMS_QDRT_H_
#define _BLAZE_MATH_SHIMS_QDRT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>


namespace blaze {

//=================================================================================================
//
//  QDRT SHIM
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Computing the quad root (4th root) of the given value/object.
// \ingroup math_shims
//
// \param a The value/object for the computation.
// \return The quad root (4th root) of the given value/object.
//
// The \a qdrt shim represents an abstract interface for computing the quad root (i.e. 4th root)
// of a value/object of any given data type.
*/
template< typename T >
BLAZE_ALWAYS_INLINE decltype(auto) qdrt( const T& a )
{
   return sqrt( sqrt( a ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
