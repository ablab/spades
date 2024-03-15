//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatStdDevExpr.h
//  \brief Header file for the dense matrix standard deviation expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSTDDEVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSTDDEVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DMatMapExpr.h>
#include <blaze/math/expressions/DMatVarExpr.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Computes the standard deviation for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the standard deviation computation.
// \return The standard deviation of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the
// <a href="https://en.wikipedia.org/wiki/Standard_deviation">standard deviation</a> for the given
// dense matrix \a dm. Example:

   \code
   using blaze::DynamicMatrix;

   DynamicMatrix<int> A{ { 1, 3, 2 }
                       , { 2, 6, 4 }
                       , { 9, 6, 3 } };

   const double s = stddev( A );  // Results in sqrt(6.5)
   \endcode

// In case the size of the given matrix is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) stddev( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return sqrt( var( *dm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the row-/columnwise standard deviation function for the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The given dense matrix for the standard deviation computation.
// \return The row-/columnwise standard deviation of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the row-/columnwise
// <a href="https://en.wikipedia.org/wiki/Standard_deviation">standard deviation</a> for the
// given dense matrix \a dm. In case \a RF is set to \a rowwise, the function returns a column
// vector containing the standard deviation of each row of \a dm. In case \a RF is set to
// \a columnwise, the function returns a row vector containing the standard deviation of each
// column of \a dm. Example:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::columnVector;
   using blaze::rowVector;

   DynamicMatrix<int> A{ { 1, 3, 2 }
                       , { 2, 6, 4 }
                       , { 9, 6, 3 } };

   DynamicVector<double,columnVector> rs;
   DynamicVector<double,rowVector> cs;

   rs = stddev<rowwise>( A );     // Results in ( 1  2  3 )
   cs = stddev<columnwise>( A );  // Results in ( sqrt(19)  sqrt(3)  1 )
   \endcode

// In case \a RF is set to \a rowwise and the number of columns of the given matrix is smaller
// than 2 or in case \a RF is set to \a columnwise and the number of rows of the given matrix is
// smaller than 2, a \a std::invalid_argument is thrown.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order
inline decltype(auto) stddev( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return sqrt( var<RF>( *dm ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
