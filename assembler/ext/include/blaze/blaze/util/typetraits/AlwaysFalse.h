//=================================================================================================
/*!
//  \file blaze/util/typetraits/AlwaysFalse.h
//  \brief Header file for the AlwaysFalse type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_ALWAYSFALSE_H_
#define _BLAZE_UTIL_TYPETRAITS_ALWAYSFALSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type dependent compile time \a false.
// \ingroup type_traits
//
// This type trait represents a type dependent compile time \a false. For all possible types,
// the \a value member constant is set to \a false, the nested type definition \a Type is
// \a FalseType, and the class derives from \a FalseType.
*/
template< typename T >
struct AlwaysFalse
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the AlwaysFalse type trait.
// \ingroup type_traits
//
// The AlwaysFalse_v variable template provides a convenient shortcut to access the nested
// \a value of the AlwaysFalse class template. For instance, given the type \a T the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::AlwaysFalse<T>::value;
   constexpr bool value2 = blaze::AlwaysFalse_v<T>;
   \endcode
*/
template< typename T >
constexpr bool AlwaysFalse_v = false;
//*************************************************************************************************

} // namespace blaze

#endif
