//=================================================================================================
/*!
//  \file blaze/system/Padding.h
//  \brief Header file for the default padding flag for all vectors and matrices of the Blaze library
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

#ifndef _BLAZE_SYSTEM_PADDING_H_
#define _BLAZE_SYSTEM_PADDING_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/config/Padding.h>
#include <blaze/math/PaddingFlag.h>


namespace blaze {

//=================================================================================================
//
//  PADDING FLAG
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default padding flag for all vectors and matrices of the Blaze library.
// \ingroup system
//
// This value specifies the default padding flag for all vectors and matrices of the Blaze
// library. In case no explicit padding flag is specified with the according vector type, this
// setting is used.

   \code
   // Explicit specification of the padding flag => the vector is padded
   StaticVector<double,3UL,columnVector,aligned,padded> a;

   // No explicit specification of the padding flag => use of the default padding flag
   StaticVector<double,3UL,columnVector,aligned> b;
   \endcode

// The default alignment is defined via the BLAZE_DEFAULT_PADDING_FLAG compilation switch (see
// the \ref padding section). Valid settings for this value are blaze::padded and blaze::unpadded.
*/
constexpr PaddingFlag defaultPaddingFlag = BLAZE_DEFAULT_PADDING_FLAG;
//*************************************************************************************************

} // namespace blaze

#endif
