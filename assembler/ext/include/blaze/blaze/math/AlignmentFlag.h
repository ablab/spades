//=================================================================================================
/*!
//  \file blaze/math/AlignmentFlag.h
//  \brief Header file for the alignment flag enumeration
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

#ifndef _BLAZE_MATH_ALIGNMENTFLAG_H_
#define _BLAZE_MATH_ALIGNMENTFLAG_H_


namespace blaze {

//=================================================================================================
//
//  ALIGNMENT FLAG
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Alignment flag for (un-)aligned vectors and matrices.
// \ingroup math
//
// Via these flags it is possible to specify subvectors, submatrices, custom vectors and matrices
// as unaligned or aligned. The following example demonstrates the setup of an unaligned subvector:

   \code
   using blaze::columnVector;
   using blaze::unaligned;

   blaze::DynamicVector<int,columnVector> v( 100UL );
   auto sv = subvector<unaligned>( v, 10UL, 20UL );
   \endcode
*/
enum AlignmentFlag : bool
{
   unaligned = false,  //!< Flag for unaligned vectors and matrices.
   aligned   = true    //!< Flag for aligned vectors and matrices.
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Negating the given alignment flag.
// \ingroup math
//
// \param flag The given alignment flag to be negated.
// \return The negated alignment flag.
//
// This logical NOT operator negates the given alignment flag. In case the given flag represents
// \a unaligned, the function returns \a aligned, in case it represents \a aligned it returns
// \a unaligned.
*/
constexpr AlignmentFlag operator!( AlignmentFlag flag ) noexcept
{
   return static_cast<AlignmentFlag>( !static_cast<bool>( flag ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
