//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatSoftmaxExpr.h
//  \brief Header file for the dense matrix softmax expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSOFTMAXEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSOFTMAXEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Column.h>
#include <blaze/math/views/Row.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Computes the softmax function for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the softmax computation.
// \return The resulting matrix.
//
// This function computes the softmax function (i.e. the normalized exponential function) for
// the given dense matrix \a dm (see also https://en.wikipedia.org/wiki/Softmax_function). The
// resulting dense matrix consists of real values in the range (0..1], which add up to 1.

   \code
   // Creating the matrix
   //    ( 1  2  3 )
   //    ( 4  1  2 )
   //    ( 3  4  1 )
   blaze::StaticMatrix<double,3UL,3UL> A{ { 1.0, 2.0, 3.0 }
                                        , { 4.0, 1.0, 2.0 }
                                        , { 3.0, 4.0, 1.0 } };

   // Computing the total softmax of A (sum(B) == 1)
   //    ( 0.0157764  0.0428847  0.116573  )
   //    ( 0.316878   0.0157764  0.0428847 )
   //    ( 0.116573   0.316878   0.0157764 )
   blaze::StaticMatrix<double,3UL,3UL> B;
   B = softmax( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
auto softmax( const DenseMatrix<MT,SO>& dm )
{
   auto tmp( evaluate( exp( *dm ) ) );
   const auto scalar( sum( tmp ) );
   tmp /= scalar;
   return tmp;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the row-/columnwise softmax function for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the softmax computation.
// \return The resulting matrix.
//
// This function computes the row-/columnwise softmax function (i.e. the normalized exponential
// function) for the given dense matrix \a dm (see also https://en.wikipedia.org/wiki/Softmax_function).
// The resulting dense matrix consists of real values in the range (0..1], which add up to the
// numbers of rows or columns, respectively.

   \code
   // Creating the matrix
   //    ( 1  2  3 )
   //    ( 4  1  2 )
   //    ( 3  4  1 )
   blaze::StaticMatrix<double,3UL,3UL> A{ { 1.0, 2.0, 3.0 }
                                        , { 4.0, 1.0, 2.0 }
                                        , { 3.0, 4.0, 1.0 } };

   // Computing the rowwise softmax of A (sum(B) == 3)
   //    ( 0.0900306  0.244728   0.665241 )
   //    ( 0.843795   0.0420101  0.114195 )
   //    ( 0.259496   0.705385   0.035119 )
   blaze::StaticMatrix<double,3UL,3UL> B;
   B = softmax<rowwise>( A );

   // Computing the columnwise softmax of A (sum(C) == 3)
   //    ( 0.035119  0.114195   0.665241  )
   //    ( 0.705385  0.0420101  0.244728  )
   //    ( 0.259496  0.843795   0.0900306 )
   blaze::StaticMatrix<double,3UL,3UL> C;
   C = softmax<columnwise>( A );
   \endcode
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
auto softmax( const DenseMatrix<MT,SO>& dm )
{
   auto tmp( evaluate( exp( *dm ) ) );

   if( RF == rowwise ) {
      for( size_t i=0UL; i<tmp.rows(); ++i ) {
         auto r = row( tmp, i, unchecked );
         const auto scalar( sum( r ) );
         r /= scalar;
      }
   }
   else {
      for( size_t j=0UL; j<tmp.columns(); ++j ) {
         auto c = column( tmp, j, unchecked );
         const auto scalar( sum( c ) );
         c /= scalar;
      }
   }

   return tmp;
}
//*************************************************************************************************

} // namespace blaze

#endif
