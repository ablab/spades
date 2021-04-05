//=================================================================================================
/*!
//  \file blaze/math/lapack/getrs.h
//  \brief Header file for the LAPACK general backward substitution functionality (getrs)
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

#ifndef _BLAZE_MATH_LAPACK_GETRS_H_
#define _BLAZE_MATH_LAPACK_GETRS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/Contiguous.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/lapack/clapack/getrs.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK LU-BASED SUBSTITUTION FUNCTIONS (GETRS)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LU-based substitution functions (getrs) */
//@{
template< typename MT, bool SO, typename VT, bool TF >
void getrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char trans,
            const blas_int_t* ipiv );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
void getrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, char trans,
            const blas_int_t* ipiv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the substitution step of solving a general linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup lapack_substitution
//
// \param A The system matrix.
// \param b The right-hand side vector.
// \param trans \c 'N' for \f$ A*x=b \f$, \c 'T' for \f$ A^T*x=b \f$, and \c C for \f$ A^H*x=b \f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Invalid right-hand side vector provided.
// \exception std::invalid_argument Invalid trans argument provided.
//
// This function uses the LAPACK getrs() functions to perform the substitution step to compute
// the solution to the general system of linear equations:
//
//  - \f$ A  *x=b \f$ if \a A is column-major
//  - \f$ A^T*x=b \f$ if \a A is row-major
//
// In this context the general system matrix \a A is a \a n-by-\a n matrix that has already been
// factorized by the getrf() functions and \a x and \a b are n-dimensional vectors. Note that the
// function only works for general, non-adapted matrices with \c float, \c double, \c complex<float>,
// or \c complex<double> element type. The attempt to call the function with adaptors or matrices
// of any other element type results in a compile time error!
//
// If the function exits successfully, the vector \a b contains the solution of the linear system
// of equations. The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a trans argument is neither \c 'N' nor \c 'T' nor \c 'C'.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::columnMajor;
   using blaze::columnVector;
   using blaze::blas_int_t;

   DynamicMatrix<double,columnMajor>  A( 2UL, 2UL );    // The system matrix A
   DynamicVector<double,columnVector> b( 2UL );         // The right-hand side vector b
   DynamicVector<blas_int_t,columnVector> ipiv( 2UL );  // Pivoting indices
   // ... Initialization

   DynamicMatrix<double,columnMajor>  D( A );  // Temporary matrix to be decomposed
   DynamicVector<double,columnVector> x( b );  // Temporary vector for the solution

   getrf( D, ipiv.data() );
   getrs( D, x, 'N', ipiv.data() );

   assert( A * x == b );
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;
   using blaze::blas_int_t;

   DynamicMatrix<double,rowMajor> A( 2UL, 2UL );        // The system matrix A
   DynamicVector<double,columnVector> b( 2UL );         // The right-hand side vector b
   DynamicVector<blas_int_t,columnVector> ipiv( 2UL );  // Pivoting indices
   // ... Initialization

   DynamicMatrix<double,rowMajor>     D( A );  // Temporary matrix to be decomposed
   DynamicVector<double,columnVector> x( b );  // Temporary vector for the solution

   getrf( D, ipiv.data() );
   getrs( D, x, 'N', ipiv.data() );

   assert( trans( A ) * x == b );
   \endcode

// For more information on the getrs() functions (i.e. sgetrs(), dgetrs(), cgetrs(), and zgetrs()),
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT  // Type of the system matrix
        , bool SO      // Storage order of the system matrix
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline void getrs( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& b, char trans,
                   const blas_int_t* ipiv )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT> );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( (*b).size() != (*A).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side vector provided" );
   }

   if( trans != 'N' && trans != 'T' && trans != 'C' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid trans argument provided" );
   }

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t nrhs( 1 );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldb ( numeric_cast<blas_int_t>( (*b).size() ) );
   blas_int_t info( 0 );

   if( n == 0 ) {
      return;
   }

   getrs( trans, n, nrhs, (*A).data(), lda, ipiv, (*b).data(), ldb, &info );

   BLAZE_INTERNAL_ASSERT( info == 0, "Invalid function argument" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the substitution step of solving a general linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_substitution
//
// \param A The system matrix.
// \param B The matrix of right-hand sides.
// \param trans \c 'N' for \f$ A*X=B \f$, \c 'T' for \f$ A^T*X=B \f$, and \c C for \f$ A^H*X=B \f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Invalid right-hand side matrix provided.
// \exception std::invalid_argument Invalid trans argument provided.
//
// This function uses the LAPACK getrs() functions to perform the substitution step to compute
// the solution to the general system of linear equations:
//
//  - \f$ A  *X  =B   \f$ if both \a A and \a B are column-major
//  - \f$ A^T*X  =B   \f$ if \a A is row-major and \a B is column-major
//  - \f$ A  *X^T=B^T \f$ if \a A is column-major and \a B is row-major
//  - \f$ A^T*X^T=B^T \f$ if both \a A and \a B are row-major
//
// In this context the general system matrix \a A is a \a n-by-\a n matrix that has already been
// factorized by the getrf() functions and \a X and \a B are either row-major \a m-by-\a n
// matrices or column-major \a n-by-\a m matrices. Note that the function only works for general,
// non-adapted matrices with \c float, \c double, \c complex<float>, or \c complex<double>
// element type. The attempt to call the function with adaptors or matrices of any other
// element type results in a compile time error!
//
// If the function exits successfully, the matrix \a B contains the solutions of the linear
// system of equations. The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the given \a trans argument is neither \c 'N' nor \c 'T' nor \c 'C'.
//  - ... the sizes of the two given matrices do not match.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::columnMajor;
   using blaze::columnVector;
   using blaze::blas_int_t;

   DynamicMatrix<double,columnMajor> A( 2UL, 2UL );     // The system matrix A
   DynamicMatrix<double,columnMajor> B( 2UL, 4UL );     // The right-hand side matrix B
   DynamicVector<blas_int_t,columnVector> ipiv( 2UL );  // Pivoting indices
   // ... Initialization

   DynamicMatrix<double,columnMajor> D( A );  // Temporary matrix to be decomposed
   DynamicMatrix<double,columnMajor> X( B );  // Temporary matrix for the solution

   getrf( D, ipiv.data() );
   getrs( D, X, 'N', ipiv.data() );

   assert( A * X == B );
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;
   using blaze::blas_int_t;

   DynamicMatrix<double,rowMajor> A( 2UL, 2UL );        // The system matrix A
   DynamicMatrix<double,rowMajor> B( 2UL, 4UL );        // The right-hand side matrix B
   DynamicVector<blas_int_t,columnVector> ipiv( 2UL );  // Pivoting indices
   // ... Initialization

   DynamicMatrix<double,rowMajor> D( A );  // Temporary matrix to be decomposed
   DynamicMatrix<double,rowMajor> X( B );  // Temporary matrix for the solution

   getrf( D, ipiv.data() );
   getrs( D, X, 'N', ipiv.data() );

   assert( trans( A ) * trans( X ) == trans( B ) );
   \endcode

// For more information on the getrs() functions (i.e. sgetrs(), dgetrs(), cgetrs(), and zgetrs()),
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline void getrs( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                   char trans, const blas_int_t* ipiv )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( trans != 'N' && trans != 'T' && trans != 'C' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid trans argument provided" );
   }

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows()    ) );
   blas_int_t mrhs( numeric_cast<blas_int_t>( SO2 ? (*B).rows() : (*B).columns() ) );
   blas_int_t nrhs( numeric_cast<blas_int_t>( SO2 ? (*B).columns() : (*B).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldb ( numeric_cast<blas_int_t>( (*B).spacing() ) );
   blas_int_t info( 0 );

   if( n != mrhs ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side matrix provided" );
   }

   if( n == 0 ) {
      return;
   }

   getrs( trans, n, nrhs, (*A).data(), lda, ipiv, (*B).data(), ldb, &info );

   BLAZE_INTERNAL_ASSERT( info == 0, "Invalid function argument" );
}
//*************************************************************************************************

} // namespace blaze

#endif
