//=================================================================================================
/*!
//  \file blaze/math/PaddingFlag.h
//  \brief Header file for the padding flag enumeration
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

#ifndef _BLAZE_MATH_PADDINGFLAG_H_
#define _BLAZE_MATH_PADDINGFLAG_H_


namespace blaze {

//=================================================================================================
//
//  PADDING FLAG VALUES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Padding flag for (un-)padded vectors and matrices.
// \ingroup math
//
// Via these flags it is possible to specify custom vectors and matrices as unpadded or padded.
// The following examples demonstrate the setup of an unaligned, unpadded and aligned, padded
// custom column vector of size 7, respectively:

   \code
   using blaze::CustomVector;
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::columnVector;

   std::vector<int> vec( 7UL );
   CustomVector<int,unaligned,unpadded,columnVector> a( &vec[0], 7UL );
   \endcode

   \code
   using blaze::CustomVector;
   using blaze::ArrayDelete;
   using blaze::aligned;
   using blaze::padded;
   using blaze::columnVector;

   std::vector<int> vec( 16UL );
   CustomVector<int,aligned,padded,columnVector> a( &vec[0], 7UL, 16UL );
   \endcode
*/
enum PaddingFlag : bool
{
   unpadded = false,  //!< Flag for unpadded vectors and matrices.
   padded   = true    //!< Flag for padded vectors and matrices.
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Negating the given padding flag.
// \ingroup math
//
// \param flag The given padding flag to be negated.
// \return The negated padding flag.
//
// This logical NOT operator negates the given padding flag. In case the given flag represents
// \a unpadded, the function returns \a padded, in case it represents \a padded it returns
// \a unpadded.
*/
constexpr PaddingFlag operator!( PaddingFlag flag ) noexcept
{
   return static_cast<PaddingFlag>( !static_cast<bool>( flag ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
