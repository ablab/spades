//=================================================================================================
/*!
//  \file blaze/util/IntegerSequence.h
//  \brief Header file for the integer_sequence and index_sequence aliases
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

#ifndef _BLAZE_UTIL_INTEGERSEQUENCE_H_
#define _BLAZE_UTIL_INTEGERSEQUENCE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/util/MaybeUnused.h>


namespace blaze {

//=================================================================================================
//
//  TYPE DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\class blaze::integer_sequence
// \brief Integer sequence type of the Blaze library.
// \ingroup util
*/
using std::integer_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::index_sequence
// \brief Index sequence type of the Blaze library.
// \ingroup util
*/
using std::index_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::make_integer_sequence
// \brief Import of the std::make_integer_sequence alias template into the Blaze namespace.
// \ingroup util
*/
using std::make_integer_sequence;
//*************************************************************************************************


//*************************************************************************************************
/*!\class blaze::make_index_sequence
// \brief Import of the std::make_index_sequence alias template into the Blaze namespace.
// \ingroup util
*/
using std::make_index_sequence;
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Equality operator for the comparison of two index sequences.
// \ingroup util
//
// \param lhs The left-hand side index sequence for the comparison.
// \param rhs The right-hand side index sequence for the comparison.
// \return \a true if the two index sequences are equal, \a false if not.
*/
template< size_t... I1s, size_t... I2s >  //
constexpr bool operator==( index_sequence<I1s...> lhs, index_sequence<I2s...> rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality operator for the comparison of two identical index sequences.
// \ingroup util
//
// \param lhs The left-hand side index sequence for the comparison.
// \param rhs The right-hand side index sequence for the comparison.
// \return \a true.
*/
template< size_t... I1s >  //
constexpr bool operator==( index_sequence<I1s...> lhs, index_sequence<I1s...> rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two index sequences.
// \ingroup util
//
// \param lhs The left-hand side index sequence for the comparison.
// \param rhs The right-hand side index sequence for the comparison.
// \return \a true if the two index sequences are not equal, \a false if they are equal.
*/
template< size_t... I1s, size_t... I2s >
constexpr bool operator!=( index_sequence<I1s...> lhs, index_sequence<I2s...> rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );

   return !( lhs == rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Shifts the given index sequence by a given offset.
// \ingroup util
//
// \param sequence The given index sequence
// \return The shifted index sequence.
*/
template< size_t Offset   // The offset for the shift operation
        , size_t... Is >  // The sequence of indices
constexpr decltype(auto) shift( std::index_sequence<Is...> sequence )
{
   MAYBE_UNUSED( sequence );

   return std::index_sequence< ( Is + Offset )... >();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creates a subsequence from the given index sequence.
// \ingroup util
//
// \param sequence The given index sequence
// \return The resulting subsequence.
*/
template< size_t... Is1    // The indices to be selected
        , size_t... Is2 >  // The sequence of indices
constexpr decltype(auto) subsequence( std::index_sequence<Is2...> sequence )
{
   MAYBE_UNUSED( sequence );

   constexpr size_t indices[] = { Is2... };
   return std::index_sequence< indices[Is1]... >();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DEFINITIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the setup of shifted index sequences.
// \ingroup util
//
// The make_shifted_index_sequence alias template provides a convenient way to create index
// sequences with specific initial index and a specific number of indices. The following code
// example demonstrates the use of make_shifted_index_sequence:

   \code
   // Creating the index sequence <2,3,4,5,6>
   using Type = make_shifted_index_sequence<2UL,5UL>;
   \endcode
*/
template< size_t Offset  // The offset of the index sequence
        , size_t N >     // The total number of indices in the index sequence
using make_shifted_index_sequence = decltype( shift<Offset>( make_index_sequence<N>() ) );
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the setup of shifted index subsequences.
// \ingroup util
//
// The make_shifted_index_subsequence alias template provides a convenient way to create a
// subsequence of an index sequences with specific initial index and a specific number of indices.
// The following code example demonstrates the use of make_shifted_index_subsequence:

   \code
   // Creating the subsequence <3,6,8> from the index sequence <2,3,4,5,6,7,8>
   using Type = make_shifted_index_subsequence<2UL,7UL,1UL,4UL,6UL>;
   \endcode
*/
template< size_t Offset    // The offset of the index sequence
        , size_t N         // The total number of indices in the index sequence
        , size_t ... Is >  // The indices to be selected
using make_shifted_index_subsequence =
   decltype( subsequence<Is...>( shift<Offset>( make_index_sequence<N>() ) ) );
//*************************************************************************************************

} // namespace blaze

#endif
