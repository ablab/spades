//=================================================================================================
/*!
//  \file blaze/math/views/Check.h
//  \brief Header file for the blaze::checked and blaze::unchecked instances
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

#ifndef _BLAZE_MATH_VIEWS_CHECK_H_
#define _BLAZE_MATH_VIEWS_CHECK_H_


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
/*!\brief Template for the blaze::checked and blaze::unchecked instances.
// \ingroup views
//
// blaze::Check is the template for the blaze::checked and blaze::unchecked instance, which is
// an optional token for the creation of views. It can be used to enforce or skip all runtime
// checks during the creation of a view (subvectors, submatrices, rows, columns, bands, ...).
*/
template< bool C >
struct Check
   : public BoolConstant<C>
{
   //**Constructor*********************************************************************************
   /*!\name Constructor */
   //@{
   constexpr Check() = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  TYPE ALIASES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of the blaze::checked instance.
// \ingroup views
//
// blaze::Checked is the type of the blaze::checked instance, which is an optional token for the
// creation of views. It can be used to enforce runtime checks during the creation of a view
// (subvectors, submatrices, rows, columns, bands, ...).
*/
using Checked = Check<true>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Type of the blaze::unchecked instance.
// \ingroup views
//
// blaze::Unchecked is the type of the blaze::unchecked instance, which is an optional token for
// the creation of views. It can be used to skip all runtime checks during the creation of a view
// (subvectors, submatrices, rows, columns, bands, ...).
*/
using Unchecked = Check<false>;
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL CHECK INSTANCES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Global Checked instance.
// \ingroup views
//
// The blaze::checked instance is an optional token for the creation of views. It can be used
// used to enforce runtime checks during the creation of a view (subvectors, submatrices, rows,
// columns, bands, ...). The following example demonstrates the setup of a checked subvector:

   \code
   blaze::DynamicVector<int> v( 100UL );
   auto sv = subvector( v, 10UL, 20UL, checked );  // Creating an checked subvector
   \endcode
*/
constexpr Checked checked;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global Unchecked instance.
// \ingroup views
//
// The blaze::unchecked instance is an optional token for the creation of views. It can be
// used to skip all runtime checks during the creation of a view (subvectors, submatrices, rows,
// columns, bands, ...). The following example demonstrates the setup of an unchecked subvector:

   \code
   blaze::DynamicVector<int> v( 100UL );
   auto sv = subvector( v, 10UL, 20UL, unchecked );  // Creating an unchecked subvector
   \endcode
*/
constexpr Unchecked unchecked;
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Termination condition for the variadic \c getCheck() function (zero arguments).
// \ingroup views
//
// \return An instance of type blaze::Checked.
*/
constexpr Checked getCheck() noexcept
{
   return Checked{};
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Termination condition for the variadic \c getCheck() function (blaze::Unchecked found).
// \ingroup views
//
// \param a The instance of type blaze::Unchecked.
// \param args The remaining arguments.
// \return An instance of type blaze::Unchecked.
*/
template< typename... Ts >
constexpr Unchecked getCheck( const Unchecked& a, const Ts&... args ) noexcept
{
   MAYBE_UNUSED( a, args... );
   return Unchecked{};
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extracting blaze::Check arguments from a given list of arguments.
// \ingroup views
//
// \param a The first given argument.
// \param args The remaining given arguments.
// \return blaze::Unchecked if at least one blaze::Unchecked is given, blaze::Checked otherwise.
//
// This function extracts any argument of type blaze::Check from the given list of arguments.
// It returns an instance of type blaze::Unchecked if at least one argument of type
// blaze::Unchecked is given, otherwise an instance of blaze::Checked.
*/
template< typename T, typename... Ts >
constexpr auto getCheck( const T& a, const Ts&... args ) noexcept
{
   MAYBE_UNUSED( a );
   return getCheck( args ... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extracting blaze::Check arguments from a given list of arguments.
// \ingroup views
//
// \param args The given arguments.
// \return \a false if at least one blaze::Unchecked is given, \a true otherwise.
//
// This function extracts any argument of type blaze::Check from the given list of arguments.
// It returns \a false if at least one argument of type blaze::Unchecked is given, otherwise
// it returns \a true.
*/
template< typename... Ts >
constexpr bool isChecked( const Ts&... args )
{
   return getCheck( args... ).value;
}
//*************************************************************************************************

} // namespace blaze

#endif
