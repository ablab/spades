//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatVarExpr.h
//  \brief Header file for the dense matrix variance expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATVAREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATVAREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Exception.h>
#include <blaze/math/dense/UniformMatrix.h>
#include <blaze/math/dense/UniformVector.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/functors/Pow2.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c var() function for general dense matrices.
// \ingroup dense_matrix
//
// \param dm The given general dense matrix for the variance computation.
// \return The variance of the given matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) var_backend( const DenseMatrix<MT,SO>& dm, FalseType )
{
   using BT = UnderlyingBuiltin_t<MT>;

   BLAZE_INTERNAL_ASSERT( size( *dm ) > 1UL, "Invalid matrix size detected" );

   const auto m( uniform<SO>( rows( *dm ), columns( *dm ), mean( *dm ) ) );

   return sum( map( (*dm) - m, Pow2() ) ) * inv( BT( size( *dm )-1UL ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c var() function for uniform dense matrices.
// \ingroup dense_matrix
//
// \param dm The given uniform dense matrix for the variance computation.
// \return The var of the given matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) var_backend( const DenseMatrix<MT,SO>& dm, TrueType )
{
   MAYBE_UNUSED( dm );

   BLAZE_INTERNAL_ASSERT( size( *dm ) > 1UL, "Invalid matrix size detected" );

   return ElementType_t<MT>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the variance for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the variance computation.
// \return The variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for
// the given dense matrix \a dm. Example:

   \code
   using blaze::DynamicMatrix;

   DynamicMatrix<int> A{ { 1, 3, 2 }
                       , { 2, 6, 4 }
                       , { 9, 6, 3 } };

   const double v = var( A );  // Results in 6.5
   \endcode

// In case the size of the given matrix is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) var( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( size( *dm ) < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix" );
   }

   return var_backend( *dm, IsUniform<MT>() );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the row-/column-wise \c var() function for general dense matrices.
// \ingroup dense_matrix
//
// \param dm The given general dense matrix for the variance computation.
// \return The row-/column-wise variance of the given matrix.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
decltype(auto) var_backend( const DenseMatrix<MT,SO>& dm, FalseType )
{
   using BT = UnderlyingBuiltin_t<MT>;

   const size_t n( RF == rowwise ? columns( *dm ) : rows( *dm ) );

   BLAZE_INTERNAL_ASSERT( n > 1UL, "Invalid matrix size detected" );

   const auto m( expand( mean<RF>( *dm ), n ) );

   return sum<RF>( map( (*dm) - m, Pow2() ) ) * inv( BT( n-1UL ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the row-/column-wise \c var() function for uniform dense matrices.
// \ingroup dense_matrix
//
// \param dm The given general dense matrix for the variance computation.
// \return The row-/column-wise variance of the given matrix.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
decltype(auto) var_backend( const DenseMatrix<MT,SO>& dm, TrueType )
{
   const size_t n( RF == rowwise ? rows( *dm ) : columns( *dm ) );

   BLAZE_INTERNAL_ASSERT( n > 0UL, "Invalid matrix size detected" );

   constexpr bool TF( ( RF == rowwise ? columnVector : rowVector ) );

   return uniform<TF>( n, ElementType_t<MT>() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the row-/column-wise variance function for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the variance computation.
// \return The row-/column-wise variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the row-/column-wise
// <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for the given dense matrix
// \a dm. In case \a RF is set to \a rowwise, the function returns a column vector containing
// the variance of each row of \a dm. In case \a RF is set to \a columnwise, the function
// returns a row vector containing the variance of each column of \a dm. Example:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::columnVector;
   using blaze::rowVector;

   DynamicMatrix<int> A{ { 1, 3, 2 }
                       , { 2, 6, 4 }
                       , { 9, 6, 3 } };

   DynamicVector<double,columnVector> rv;
   DynamicVector<double,rowVector> cv;

   rv = var<rowwise>( A );     // Results in ( 1  4  9 )
   cv = var<columnwise>( A );  // Results in ( 19  3  1 )
   \endcode

// In case \a RF is set to \a rowwise and the number of columns of the given matrix is smaller
// than 2 or in case \a RF is set to \a columnwise and the number of rows of the given matrix is
// smaller than 2, a \a std::invalid_argument is thrown.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
inline decltype(auto) var( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( RF < 2UL, "Invalid reduction flag" );

   const size_t n( RF == rowwise ? columns( *dm ) : rows( *dm ) );

   if( n < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix" );
   }

   return var_backend<RF>( *dm, IsUniform<MT>() );
}
//*************************************************************************************************

} // namespace blaze

#endif
