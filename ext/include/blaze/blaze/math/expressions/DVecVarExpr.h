//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecVarExpr.h
//  \brief Header file for the dense vector variance expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECVAREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECVAREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/dense/UniformVector.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/functors/Pow2.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/IsUniform.h>
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
/*!\brief Backend implementation of the \c var() function for general dense vectors.
// \ingroup dense_vector
//
// \param dv The given general dense vector for the variance computation.
// \return The variance of the given vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) var_backend( const DenseVector<VT,TF>& dv, FalseType )
{
   using BT = UnderlyingBuiltin_t<VT>;

   const size_t n( size( *dv ) );

   BLAZE_INTERNAL_ASSERT( n > 1UL, "Invalid vector size detected" );

   const auto m( uniform<TF>( n, mean( *dv ) ) );

   return sum( map( (*dv) - m, Pow2() ) ) * inv( BT( n-1UL ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c var() function for uniform dense vectors.
// \ingroup dense_vector
//
// \param dv The given uniform dense vector for the variance computation.
// \return The variance of the given vector.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) var_backend( const DenseVector<VT,TF>& dv, TrueType )
{
   MAYBE_UNUSED( dv );

   BLAZE_INTERNAL_ASSERT( size( *dv ) > 1UL, "Invalid vector size detected" );

   return ElementType_t<VT>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the variance for the given dense vector.
// \ingroup dense_vector
//
// \param dv The given dense vector for the variance computation.
// \return The variance of the given vector.
// \exception std::invalid_argument Invalid input vector.
//
// This function computes the <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for
// the given dense vector \a dv. Example:

   \code
   using blaze::DynamicVector;

   DynamicVector<int> v{ 1, 4, 3, 6, 7 };

   const double m = var( v );  // Results in 5.7
   \endcode

// In case the size of the given vector is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
decltype(auto) var( const DenseVector<VT,TF>& dv )
{
   BLAZE_FUNCTION_TRACE;

   const size_t n( size( *dv ) );

   if( n < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input vector" );
   }

   return var_backend( *dv, IsUniform<VT>() );
}
//*************************************************************************************************

} // namespace blaze

#endif
