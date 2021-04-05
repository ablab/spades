//=================================================================================================
/*!
//  \file blaze/util/typelist/Length.h
//  \brief Header file for the Length class template
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

#ifndef _BLAZE_UTIL_TYPELIST_LENGTH_H_
#define _BLAZE_UTIL_TYPELIST_LENGTH_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>
#include <blaze/util/typelist/TypeList.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculating the length of a type list.
// \ingroup typelist
//
// The Length class can be used to obtain the length of a type list (i.e. the number
// of contained types). In order to obtain the length of a type list, the Length class
// has to be instantiated for a particular type list. The length of the type list can
// be obtained using the member enumeration \a value. The following example gives an
// impression of the use of the Length class:

   \code
   using Floats = blaze::TypeList< float, double, long double >;  // Defining a new type list
   const int length = blaze::Length< Floats >::value;             // The length of the type list
   \endcode
*/
template< typename TL >  // Type of the type list
struct Length;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Spezialization of the Length class for type lists.
// \ingroup typelist
*/
template< typename... Ts >  // Type list elements
struct Length< TypeList<Ts...> >
   : public Size_t< sizeof...( Ts ) >
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the Length type trait.
// \ingroup type_traits
//
// The Length_v variable template provides a convenient shortcut to access the nested \a value of
// the Length class template. For instance, given the type list \a TL the following two statements
// are identical:

   \code
   constexpr size_t value1 = Length<TL>::value;
   constexpr size_t value2 = Length_v<TL>;
   \endcode
*/
template< typename TL >  // Type of the type list
constexpr size_t Length_v = Length<TL>::value;
//*************************************************************************************************

} // namespace blaze

#endif
