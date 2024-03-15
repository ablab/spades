//=================================================================================================
/*!
//  \file blaze/math/RelaxationFlag.h
//  \brief Header file for the relaxation flag enumeration
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

#ifndef _BLAZE_MATH_RELAXATIONFLAG_H_
#define _BLAZE_MATH_RELAXATIONFLAG_H_


namespace blaze {

//=================================================================================================
//
//  RELAXATION FLAG
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Relaxation flag for strict or relaxed semantics.
// \ingroup math
//
// Via these flags it is possible to specify that according operations should use strict
// semantics instead of relaxed semantics or vice versa. The following example demonstrates this
// by means of the isDefault() function template:

   \code
   using blaze::strict;
   using blaze::relaxed;

   blaze::StaticVector<double,3UL> v{ 0.0, 1E-9, 0.0 };

   isDefault<strict> ( v );  // Returns false
   isDefault<relaxed>( v );  // Returns true
   \endcode
*/
enum RelaxationFlag : bool
{
   strict  = false,  //!< Flag for strict semantics.
   relaxed = true    //!< Flag for relaxed semantics.
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Negating the given relaxation flag.
// \ingroup math
//
// \param flag The given relaxation flag to be negated.
// \return The negated relaxation flag.
//
// This logical NOT operator negates the given relaxation flag. In case the given flag represents
// \a strict, the function returns \a relaxed, in case it represents \a relaxed it returns
// \a strict.
*/
constexpr RelaxationFlag operator!( RelaxationFlag flag ) noexcept
{
   return static_cast<RelaxationFlag>( !static_cast<bool>( flag ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
