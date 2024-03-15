//=================================================================================================
/*!
//  \file blaze/system/Compiler.h
//  \brief Compiler-specific system settings
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

#ifndef _BLAZE_SYSTEM_COMPILER_H_
#define _BLAZE_SYSTEM_COMPILER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>




//=================================================================================================
//
//  INTEL COMPILER MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(__INTEL_COMPILER) || defined(__ICL) || defined(__ICC) || defined(__ECC)
#  define BLAZE_INTEL_COMPILER 1
#else
#  define BLAZE_INTEL_COMPILER 0
#endif
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLANG COMPILER MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(__clang__) && !BLAZE_INTEL_COMPILER
#  define BLAZE_CLANG_COMPILER 1
#else
#  define BLAZE_CLANG_COMPILER 0
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#define BLAZE_CLANG_MAJOR_VERSION __clang_major__
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#define BLAZE_CLANG_MINOR_VERSION __clang_minor__
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#define BLAZE_CLANG_PATCH_VERSION __clang_patchlevel__
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GNU COMPILER MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(__GNUC__) && !BLAZE_CLANG_COMPILER && !BLAZE_INTEL_COMPILER
#  define BLAZE_GNU_COMPILER 1
#else
#  define BLAZE_GNU_COMPILER 0
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#define BLAZE_GNU_MAJOR_VERSION __GNUC__
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#define BLAZE_GNU_MINOR_VERSION __GNUC_MINOR__
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MICROSOFT COMPILER MACRO DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if defined(_MSC_VER) && !BLAZE_CLANG_COMPILER && !BLAZE_INTEL_COMPILER
#  define BLAZE_MSC_COMPILER 1
#else
#  define BLAZE_MSC_COMPILER 0
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

BLAZE_STATIC_ASSERT( !BLAZE_GNU_COMPILER || ( !BLAZE_CLANG_COMPILER && !BLAZE_MSC_COMPILER && !BLAZE_INTEL_COMPILER ) );
BLAZE_STATIC_ASSERT( !BLAZE_CLANG_COMPILER || ( !BLAZE_GNU_COMPILER && !BLAZE_MSC_COMPILER && !BLAZE_INTEL_COMPILER ) );
BLAZE_STATIC_ASSERT( !BLAZE_MSC_COMPILER || ( !BLAZE_GNU_COMPILER && !BLAZE_CLANG_COMPILER && !BLAZE_INTEL_COMPILER ) );
BLAZE_STATIC_ASSERT( !BLAZE_INTEL_COMPILER || ( !BLAZE_GNU_COMPILER && !BLAZE_CLANG_COMPILER && !BLAZE_MSC_COMPILER ) );

}
/*! \endcond */
//*************************************************************************************************

#endif
