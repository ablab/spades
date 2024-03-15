//=================================================================================================
/*!
//  \file blaze/util/typetraits/RemoveCVRef.h
//  \brief Header file for the RemoveCVRef type trait
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

#ifndef _BLAZE_UTIL_TYPETRAITS_REMOVECVREF_H_
#define _BLAZE_UTIL_TYPETRAITS_REMOVECVREF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <type_traits>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Removal of top level cv-qualifiers.
// \ingroup type_traits
//
// The RemoveCVRef type trait removes any reference modifiers and all top level cv-qualifiers
// from the given type \a T.

   \code
   blaze::RemoveCVRef<short>::Type               // Results in 'short'
   blaze::RemoveCVRef<int&>::Type                // Results in 'int'
   blaze::RemoveCVRef<const double>::Type        // Results in 'double'
   blaze::RemoveCVRef<volatile float>::Type      // Results in 'float'
   blaze::RemoveCVRef<const volatile int>::Type  // Results in 'int'
   blaze::RemoveCVRef<const int&>::Type          // Results in 'int'
   blaze::RemoveCVRef<const int*>::Type          // Results in 'int const*'
   blaze::RemoveCVRef<const int* const>::Type    // Results in 'int const*'
   \endcode
*/
template< typename T >
struct RemoveCVRef
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename std::remove_cv< typename std::remove_reference<T>::type >::type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the RemoveCVRef type trait.
// \ingroup type_traits
//
// The RemoveCVRef_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the RemoveCVRef class template. For instance, given the type \a T the following two type
// definitions are identical:

   \code
   using Type1 = typename blaze::RemoveCVRef<T>::Type;
   using Type2 = blaze::RemoveCVRef_t<T>;
   \endcode
*/
template< typename T >
using RemoveCVRef_t = typename RemoveCVRef<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
