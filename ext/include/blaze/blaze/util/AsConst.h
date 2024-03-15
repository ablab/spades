//=================================================================================================
/*!
//  \file blaze/util/AsConst.h
//  \brief Header file for the as_const function template
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

#ifndef _BLAZE_UTIL_ASCONST_H_
#define _BLAZE_UTIL_ASCONST_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/typetraits/AddConst.h>


namespace blaze {

//=================================================================================================
//
//  ASCONST FUNCTIONALITY
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Adding 'const' to the given lvalue.
// \ingroup util
//
// \param v The given lvalue.
// \return The const-qualified lvalue.
//
// This function adds the 'const' qualifier to the given lvalue.
*/
template< typename T >
constexpr AddConst_t<T>& as_const( T& v ) noexcept
{
   return v;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Overload of the as_const() function for rvalues.
// \ingroup util
//
// This overload of the as_const() function disables its use on rvalues. This prevents potential
// misuse as in for instance the following example:

   \code
   for( const auto&& value : as_const( getTemporary() ) )
   {
      // ...
   }
   \endcode
*/
template< typename T >
void as_const( const T&& ) = delete;
//*************************************************************************************************

} // namespace blaze

#endif
