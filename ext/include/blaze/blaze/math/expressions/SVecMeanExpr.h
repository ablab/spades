//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecMeanExpr.h
//  \brief Header file for the sparse vector mean expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECMEANEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECMEANEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c mean() function for general sparse vectors.
// \ingroup sparse_vector
//
// \param sv The given general sparse vector for the mean computation.
// \return The mean of the given vector.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) mean_backend( const SparseVector<VT,TF>& sv, FalseType )
{
   using BT = UnderlyingBuiltin_t<VT>;

   BLAZE_INTERNAL_ASSERT( size( *sv ) > 0UL, "Invalid vector size detected" );

   return sum( *sv ) * inv( BT( size( *sv ) ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c mean() function for uniform sparse vectors.
// \ingroup sparse_vector
//
// \param sv The given uniform sparse vector for the mean computation.
// \return The mean of the given vector.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) mean_backend( const SparseVector<VT,TF>& sv, TrueType )
{
   MAYBE_UNUSED( sv );

   BLAZE_INTERNAL_ASSERT( size( *sv ) > 0UL, "Invalid vector size detected" );

   return ElementType_t<VT>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the (arithmetic) mean for the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector for the mean computation.
// \return The mean of the given vector.
// \exception std::invalid_argument Invalid input vector.
//
// This function computes the
// <a href="https://en.wikipedia.org/wiki/Arithmetic_mean">(arithmetic) mean</a> for the given
// sparse vector \a sv. Both the non-zero and zero elements of the sparse vector are taken into
// account. Example:

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v{ 1, 0, 4, 0, 3, 0, 6, 0, 7, 0 };

   const double m = mean( v );  // Results in 2.1 (i.e. 21/10 )
   \endcode

// In case the size of the given vector is 0, a \a std::invalid_argument is thrown.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline decltype(auto) mean( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   if( size( *sv ) == 0UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input vector" );
   }

   return mean_backend( *sv, IsZero<VT>() );
}
//*************************************************************************************************

} // namespace blaze

#endif
