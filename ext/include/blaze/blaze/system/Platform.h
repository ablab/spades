//=================================================================================================
/*!
//  \file blaze/system/Platform.h
//  \brief Platform-specific system settings
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

#ifndef _BLAZE_SYSTEM_PLATFORM_H_
#define _BLAZE_SYSTEM_PLATFORM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>




//=================================================================================================
//
//  PLATFORM MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(_WIN32)
#  define BLAZE_WIN32_PLATFORM 1
#else
#  define BLAZE_WIN32_PLATFORM 0
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(_WIN64)
#  define BLAZE_WIN64_PLATFORM 1
#else
#  define BLAZE_WIN64_PLATFORM 0
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(__MINGW64_VERSION_MAJOR)
#  define BLAZE_MINGW64_PLATFORM 1
#else
#  define BLAZE_MINGW64_PLATFORM 0
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(__MINGW32__)
#  define BLAZE_MINGW32_PLATFORM 1
#else
#  define BLAZE_MINGW32_PLATFORM 0
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( !BLAZE_WIN64_PLATFORM || ( !BLAZE_MINGW32_PLATFORM && !BLAZE_MINGW64_PLATFORM ) );
BLAZE_STATIC_ASSERT( !BLAZE_MINGW32_PLATFORM || ( !BLAZE_WIN64_PLATFORM && !BLAZE_MINGW64_PLATFORM ) );
BLAZE_STATIC_ASSERT( !BLAZE_MINGW64_PLATFORM || ( !BLAZE_WIN64_PLATFORM && !BLAZE_MINGW32_PLATFORM ) );

}
/*! \endcond */
//*************************************************************************************************

#endif
