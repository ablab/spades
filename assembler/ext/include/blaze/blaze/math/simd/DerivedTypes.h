//=================================================================================================
/*!
//  \file blaze/math/simd/DerivedTypes.h
//  \brief Header file for the derived SIMD types
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

#ifndef _BLAZE_MATH_SIMD_DERIVEDTYPES_H_
#define _BLAZE_MATH_SIMD_DERIVEDTYPES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/system/Vectorization.h>


namespace blaze {

//=================================================================================================
//
//  DERIVED SIMD TYPES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The SIMD data type for 'char'.
// \ingroup simd
*/
using SIMDchar = SIMDTrait<char>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'signed char'.
// \ingroup simd
*/
using SIMDschar = SIMDTrait<signed char>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned char'.
// \ingroup simd
*/
using SIMDuchar = SIMDTrait<unsigned char>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'wchar_t'.
// \ingroup simd
*/
using SIMDwchar = SIMDTrait<wchar_t>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<char>'.
// \ingroup simd
*/
using SIMDcchar = SIMDTrait< complex<char> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<signed char>'.
// \ingroup simd
*/
using SIMDcschar = SIMDTrait< complex<signed char> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned char>'.
// \ingroup simd
*/
using SIMDcuchar = SIMDTrait< complex<unsigned char> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<wchar_t>'.
// \ingroup simd
*/
using SIMDcwchar = SIMDTrait< complex<wchar_t> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'short'.
// \ingroup simd
*/
using SIMDshort = SIMDTrait<short>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned short'.
// \ingroup simd
*/
using SIMDushort = SIMDTrait<unsigned short>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<short>'.
// \ingroup simd
*/
using SIMDcshort = SIMDTrait< complex<short> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned short>'.
// \ingroup simd
*/
using SIMDcushort = SIMDTrait< complex<unsigned short> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'int'.
// \ingroup simd
*/
using SIMDint = SIMDTrait<int>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned int'.
// \ingroup simd
*/
using SIMDuint = SIMDTrait<unsigned int>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<int>'.
// \ingroup simd
*/
using SIMDcint = SIMDTrait< complex<int> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned int>'.
// \ingroup simd
*/
using SIMDcuint = SIMDTrait< complex<unsigned int> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'long int'.
// \ingroup simd
*/
using SIMDlong = SIMDTrait<long>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'unsigned long int'.
// \ingroup simd
*/
using SIMDulong = SIMDTrait<unsigned long>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<long int>'.
// \ingroup simd
*/
using SIMDclong = SIMDTrait< complex<long> >::Type;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The SIMD data type for 'complex<unsigned long int>'.
// \ingroup simd
*/
using SIMDculong = SIMDTrait< complex<unsigned long> >::Type;
//*************************************************************************************************

} // namespace blaze

#endif
