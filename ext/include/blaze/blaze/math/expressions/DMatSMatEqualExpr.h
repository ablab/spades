//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatSMatEqualExpr.h
//  \brief Header file for the dense matrix/sparse matrix equality comparison expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSMATEQUALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSMATEQUALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL BINARY RELATIONAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality operator for the comparison of a dense matrix and a row-major sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side row-major sparse matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of a dense matrix and a sparse matrix. Due to the limited
// machine accuracy, a direct comparison of two floating point numbers should be avoided. This
// function offers the possibility to compare two floating-point matrices with a certain accuracy
// margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , bool SO            // Storage order of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side sparse matrix
inline bool equal( const DenseMatrix<MT1,SO>& lhs, const SparseMatrix<MT2,false>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the dense matrix and sparse matrix operand
   CT1 A( *lhs );
   CT2 B( *rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   size_t j( 0UL );

   for( size_t i=0UL; i<B.rows(); ++i ) {
      j = 0UL;
      for( auto element=B.begin(i); element!=B.end(i); ++element, ++j ) {
         for( ; j<element->index(); ++j ) {
            if( !isDefault<RF>( A(i,j) ) ) return false;
         }
         if( !equal<RF>( element->value(), A(i,j) ) ) return false;
      }
      for( ; j<A.columns(); ++j ) {
         if( !isDefault<RF>( A(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality operator for the comparison of a dense matrix and a column-major sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side column-major sparse matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of a dense matrix and a sparse matrix. Due to the limited
// machine accuracy, a direct comparison of two floating point numbers should be avoided. This
// function offers the possibility to compare two floating-point matrices with a certain accuracy
// margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , bool SO            // Storage order of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side sparse matrix
inline bool equal( const DenseMatrix<MT1,SO>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the dense matrix and sparse matrix operand
   CT1 A( *lhs );
   CT2 B( *rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   size_t i( 0UL );

   for( size_t j=0UL; j<B.columns(); ++j ) {
      i = 0UL;
      for( auto element=B.begin(j); element!=B.end(j); ++element, ++i ) {
         for( ; i<element->index(); ++i ) {
            if( !isDefault<RF>( A(i,j) ) ) return false;
         }
         if( !equal<RF>( element->value(), A(i,j) ) ) return false;
      }
      for( ; i<A.rows(); ++i ) {
         if( !isDefault<RF>( A(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality operator for the comparison of a sparse matrix and a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of a sparse matrix and a dense matrix. Due to the limited
// machine accuracy, a direct comparison of two floating point numbers should be avoided. This
// function offers the possibility to compare two floating-point matrices with a certain accuracy
// margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side sparse matrix
        , bool SO1           // Storage order of the left-hand side sparse matrix
        , typename MT2       // Type of the right-hand side dense matrix
        , bool SO2 >         // Storage order of the right-hand side dense matrix
inline bool equal( const SparseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return equal<RF>( rhs, lhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense matrix and a sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator==( const DenseMatrix<MT1,SO1>& lhs, const SparseMatrix<MT2,SO2>& rhs )
{
   return equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a sparse matrix and a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , bool SO1      // Storage order of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator==( const SparseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return equal<relaxed>( rhs, lhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense matrix and a sparse matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side sparse matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SO2 >    // Storage order of the right-hand side sparse matrix
inline bool operator!=( const DenseMatrix<MT1,SO1>& lhs, const SparseMatrix<MT2,SO2>& rhs )
{
   return !equal( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a sparse matrix and a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side sparse matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , bool SO1      // Storage order of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order right-hand side dense matrix
inline bool operator!=( const SparseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return !equal( rhs, lhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
