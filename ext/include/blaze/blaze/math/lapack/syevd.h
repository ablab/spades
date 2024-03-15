//=================================================================================================
/*!
//  \file blaze/math/lapack/syevd.h
//  \brief Header file for the LAPACK symmetric matrix eigenvalue functions (syevd)
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

#ifndef _BLAZE_MATH_LAPACK_SYEVD_H_
#define _BLAZE_MATH_LAPACK_SYEVD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Contiguous.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/lapack/clapack/syevd.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK SYMMETRIC MATRIX EIGENVALUE FUNCTIONS (SYEVD)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK symmetric matrix eigenvalue functions (syevd) */
//@{
template< typename MT, bool SO, typename VT, bool TF >
void syevd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense symmetric matrix.
// \ingroup lapack_eigenvalue
//
// \param A The given symmetric matrix.
// \param jobz \c 'V' to compute the eigenvectors of \a A, \c 'N' to only compute the eigenvalues.
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param w The resulting vector of eigenvalues.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Invalid jobz argument provided.
// \exception std::invalid_argument Invalid uplo argument provided.
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function computes the eigenvalues of a symmetric \a n-by-\a n matrix based on the LAPACK
// syevd() functions, which use a divide-and-conquer strategy. Optionally, it computes the left
// or right eigenvectors.
//
// The real eigenvalues are returned in ascending order in the given vector \a w. \a w is resized
// to the correct size (if possible and necessary). In case \a A is a row-major matrix, the left
// eigenvectors are returned in the rows of \a A, in case \a A is a column-major matrix, the right
// eigenvectors are returned in the columns of \a A.
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given \a jobz argument is neither \c 'V' nor \c 'N';
//  - ... the given \a uplo argument is neither \c 'L' nor \c 'U';
//  - ... the eigenvalue computation fails.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The symmetric matrix A
   // ... Initialization

   DynamicVector<double,columnVector> w( 5UL );  // The vector for the real eigenvalues

   syevd( A, w, 'V', 'L' );
   \endcode

// For more information on the syevd() functions (i.e. ssyevd() and dsyevd()) see the LAPACK
// online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector w
        , bool TF >    // Transpose flag of the vector w
inline void syevd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w, char jobz, char uplo )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ElementType_t<MT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ElementType_t<VT> );

   using ET = ElementType_t<MT>;

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( jobz != 'V' && jobz != 'N' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid jobz argument provided" );
   }

   if( uplo != 'L' && uplo != 'U' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uplo argument provided" );
   }

   resize( *w, (*A).rows(), false );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows()    ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   if( n == 0 ) {
      return;
   }

   blas_int_t lwork( 2*n*n + 6*n + 3 );
   const std::unique_ptr<ET[]> work( new ET[lwork] );

   blas_int_t liwork( 3 + 5*n );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[liwork] );

   if( IsRowMajorMatrix_v<MT> ) {
      ( uplo == 'L' )?( uplo = 'U' ):( uplo = 'L' );
   }

   syevd( jobz, uplo, n, (*A).data(), lda, (*w).data(),
          work.get(), lwork, iwork.get(), liwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
