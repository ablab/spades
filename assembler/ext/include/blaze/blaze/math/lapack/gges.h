//=================================================================================================
/*!
//  \file blaze/math/lapack/gges.h
//  \brief Header file for the LAPACK general matrix eigenvalue functions (gges)
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

#ifndef _BLAZE_MATH_LAPACK_GGES_H_
#define _BLAZE_MATH_LAPACK_GGES_H_


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
#include <blaze/math/lapack/clapack/gges.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/constraints/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK GENERALIZED MATRIX EIGENVALUE FUNCTIONS (GGES)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK generalized Schur factorization functions (gges) */
//@{
template< typename MT1, bool SO1, typename MT2, bool SO2
        , typename VT1, bool TF1, typename VT2, bool TF2 >
void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
           DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta );

template< typename MT1, bool SO1, typename MT2, bool SO2
        , typename VT1, bool TF1, typename VT2, bool TF2
        , typename Select >
void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
           DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta,
           Select select );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3
        , typename VT1, bool TF1, typename VT2, bool TF2, typename MT4, bool SO4 >
void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, DenseMatrix<MT3,SO3>& VSL,
           DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3
        , typename VT1, bool TF1, typename VT2, bool TF2, typename MT4, bool SO4
        , typename Select >
void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B, DenseMatrix<MT3,SO3>& VSL,
           DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR,
           Select select );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gges kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of real general matrices.
// \param B The second matrix of the given pair of real general matrices.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting real vector of eigenvalue denominators.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::runtime_error Schur factorization computation failed.
//
// This function is the backend implementation for computing the generalized Schur factorization
// of the given pair of real dense general matrices.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gges() function.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename Select >  // Type of the eigenvalue selector
inline auto gges_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                          DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta,
                          Select select )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *B ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*alpha).size() == (*A).rows(), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*beta).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT1>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).rows() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb( numeric_cast<blas_int_t>( (*B).spacing() ) );
   blas_int_t info( 0 );
   blas_int_t sdim( 0 );

   blas_int_t lwork = max( 8*n, 6*n + 16 );
   const std::unique_ptr<BT[]> alphar( new BT[n] );
   const std::unique_ptr<BT[]> alphai( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[max(1,lwork)] );
   const std::unique_ptr<blas_int_t[]> bwork( select ? new blas_int_t[n] : nullptr );

   gges( 'N', 'N', ( select ? 'S' : 'N' ), select, n, (*A).data(), lda, (*B).data(), ldb, &sdim,
         alphar.get(), alphai.get(), (*beta).data(), nullptr, 1, nullptr, 1,
         work.get(), lwork, bwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for generalized eigenvalue decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Generalized eigenvalue decomposition failed" );
   }

   for( size_t i=0UL; i<(*A).rows(); ++i ) {
      (*alpha)[i] = CT( alphar[i], alphai[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gges kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::runtime_error Schur factorization computation failed.
//
// This function is the backend implementation for computing the generalized Schur factorization
// of the given pair of complex dense general matrices.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gges() function.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename Select >  // Type of the eigenvalue selector
inline auto gges_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                          DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta,
                          Select select )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *B ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*alpha).size() == (*A).rows(), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*beta).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT1>;
   using BT = typename CT::value_type;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).rows() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb( numeric_cast<blas_int_t>( (*B).spacing() ) );
   blas_int_t info( 0 );
   blas_int_t sdim( 0 );

   blas_int_t lwork = max( 1, 2*n );
   const std::unique_ptr<CT[]> work( new CT[max(1,lwork)] );
   const std::unique_ptr<BT[]> rwork( new BT[8*n] );
   const std::unique_ptr<blas_int_t[]> bwork( select ? new blas_int_t[n] : nullptr );

   gges( 'N', 'N', ( select ? 'S' : 'N' ), select, n, (*A).data(), lda, (*B).data(), ldb, &sdim,
         (*alpha).data(), (*beta).data(), nullptr, 1, nullptr, 1,
         work.get(), lwork, rwork.get(), bwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for generalized eigenvalue decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Generalized eigenvalue decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization of the given pair of
//        dense general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of \a n-by-\a n real non-symmetric matrices (\a A,\a B), the
// generalized eigenvalues and the generalized real Schur form (\a S,\a T), which on exit replaces
// (\a A, \a B).
//
// A generalized eigenvalue for a pair of matrices (\a A,\a B) is a scalar \a w or a ratio
// \f$ alpha/beta = w \f$, such that \f$ A - w*B \f$ is singular. It is usually represented as
// the pair (\a alpha,\a beta), as there is a reasonable interpretation for \a beta=0 or both
// being zero. The complex eigenvalues are returned as numerators and denominators in the given
// vectors \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the size of the given matrices \a A and \a B don't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The first general matrix
   DynamicMatrix<double,rowMajor> B( 5UL, 5UL );  // The second general matrix
   // ... Initialization

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The vector for the eigenvalue numerators
   DynamicVector<double,columnVector> beta( 5UL );            // The vector for the eigenvalue denominators

   gges( A, B, alpha, beta );
   \endcode

// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix B
        , bool SO2      // Storage order of the matrix B
        , typename VT1  // Type of the vector alpha
        , bool TF1      // Transpose flag of the vector alpha
        , typename VT2  // Type of the vector beta
        , bool TF2 >    // Transpose flag of the vector beta
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                  DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *B, N, N, false );
   resize( *alpha, N, false );
   resize( *beta, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( *A, *B, *alpha, *beta, nullptr );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization of the given pair of
//        dense general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of \a n-by-\a n real non-symmetric matrices (\a A,\a B), the
// generalized eigenvalues and the generalized real Schur form (\a S,\a T), which on exit replaces
// (\a A, \a B).
//
// A generalized eigenvalue for a pair of matrices (\a A,\a B) is a scalar \a w or a ratio
// \f$ alpha/beta = w \f$, such that \f$ A - w*B \f$ is singular. It is usually represented as
// the pair (\a alpha,\a beta), as there is a reasonable interpretation for \a beta=0 or both
// being zero. The complex eigenvalues are returned as numerators and denominators in the given
// vectors \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// It also orders the eigenvalues such that a selected cluster of eigenvalues appears in the
// leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The according selection of eigenvalues is performed via the given \a select
// function. In case the given pair of matrices is a pair of real matrices the function is
// required to accept three pointers to real arguments, in case the given pair of matrices
// is a pair of complex matrices the function is required to accept two pointers to complex
// arguments.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the size of the given matrices \a A and \a B don't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The first general matrix
   DynamicMatrix<double,rowMajor> B( 5UL, 5UL );  // The second general matrix
   // ... Initialization

   const auto select = []( T* alphar, T* alphai, T* beta ) -> blas_int_t {
      return *alphar > 0.0;
   };

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The vector for the eigenvalue numerators
   DynamicVector<double,columnVector> beta( 5UL );            // The vector for the eigenvalue denominators

   gges( A, B, alpha, beta, select );
   \endcode

// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename Select >  // Type of the eigenvalue selector
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                  DenseVector<VT1,TF1>& alpha, DenseVector<VT2,TF2>& beta,
                  Select select )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *B, N, N, false );
   resize( *alpha, N, false );
   resize( *beta, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( *A, *B, *alpha, *beta, select );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gges kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of real general matrices.
// \param B The second matrix of the given pair of real general matrices.
// \param VSL The resulting real matrix of left Schur vectors.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting real vector of eigenvalue denominators.
// \param VSR The resulting real matrix of right Schur vectors.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::runtime_error Schur factorization computation failed.
//
// This function is the backend implementation for computing the generalized Schur factorization
// of the given pair of real dense general matrices.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gges() function.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename MT3       // Type of the matrix VSL
        , bool SO3           // Storage order of the matrix VSL
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename MT4       // Type of the matrix VSR
        , bool SO4           // Storage order of the matrix VSR
        , typename Select >  // Type of the eigenvalue selector
inline auto gges_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                          DenseMatrix<MT3,SO3>& VSL, DenseVector<VT1,TF1>& alpha,
                          DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR,
                          Select select )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *B ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VSL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VSR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VSL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VSR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*alpha).size() == (*A).rows(), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*beta).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT1>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   const blas_int_t n    ( numeric_cast<blas_int_t>( (*A).rows() ) );
   const blas_int_t lda  ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb  ( numeric_cast<blas_int_t>( (*B).spacing() ) );
   const blas_int_t ldvsl( numeric_cast<blas_int_t>( (*VSL).spacing() ) );
   const blas_int_t ldvsr( numeric_cast<blas_int_t>( (*VSR).spacing() ) );
   blas_int_t info( 0 );
   blas_int_t sdim( 0 );

   blas_int_t lwork = max( 8*n, 6*n + 16 );
   const std::unique_ptr<BT[]> alphar( new BT[n] );
   const std::unique_ptr<BT[]> alphai( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[max(1, lwork)] );
   const std::unique_ptr<blas_int_t[]> bwork( select ? new blas_int_t[n] : nullptr );

   gges( 'V', 'V', ( select ? 'S' : 'N' ), select, n, (*A).data(), lda, (*B).data(), ldb, &sdim,
         alphar.get(), alphai.get(), (*beta).data(), (*VSL).data(), ldvsl, (*VSR).data(), ldvsr,
         work.get(), lwork, bwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for generalized eigenvalue decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Generalized eigenvalue decomposition failed" );
   }

   for( size_t i=0UL; i<(*A).rows(); ++i ) {
      (*alpha)[i] = CT( alphar[i], alphai[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gges kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param VSL The resulting complex matrix of left Schur vectors.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \param VSR The resulting complex matrix of right Schur vectors.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::runtime_error Schur factorization computation failed.
//
// This function is the backend implementation for computing the generalized Schur factorization
// of the given pair of complex dense general matrices.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gges() function.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename MT3       // Type of the matrix VSL
        , bool SO3           // Storage order of the matrix VSL
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename MT4       // Type of the matrix VSR
        , bool SO4           // Storage order of the matrix VSR
        , typename Select >  // Type of the eigenvalue selector
inline auto gges_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                          DenseMatrix<MT3,SO3>& VSL, DenseVector<VT1,TF1>& alpha,
                          DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR,
                          Select select )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *B ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VSL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VSR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VSL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VSR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*alpha).size() == (*A).rows(), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*beta).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT1>;
   using BT = typename CT::value_type;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   const blas_int_t n    ( numeric_cast<blas_int_t>( (*A).rows() ) );
   const blas_int_t lda  ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb  ( numeric_cast<blas_int_t>( (*B).spacing() ) );
   const blas_int_t ldvsl( numeric_cast<blas_int_t>( (*VSL).spacing() ) );
   const blas_int_t ldvsr( numeric_cast<blas_int_t>( (*VSR).spacing() ) );
   blas_int_t info( 0 );
   blas_int_t sdim( 0 );

   blas_int_t lwork = max( 1, 2*n );
   const std::unique_ptr<CT[]> work( new CT[max(1,lwork)] );
   const std::unique_ptr<BT[]> rwork( new BT[8*n] );
   const std::unique_ptr<blas_int_t[]> bwork( select ? new blas_int_t[n] : nullptr );

   gges( 'V', 'V', ( select ? 'S' : 'N' ), select, n, (*A).data(), lda, (*B).data(), ldb, &sdim,
         (*alpha).data(), (*beta).data(), (*VSL).data(), ldvsl, (*VSR).data(), ldvsr,
         work.get(), lwork, rwork.get(), bwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for generalized eigenvalue decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Generalized eigenvalue decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization of the given pair of
//        dense general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param VSL The resulting complex matrix of left Schur vectors.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \param VSR The resulting complex matrix of right Schur vectors.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of \a n-by-\a n real non-symmetric matrices (\a A,\a B), the
// generalized eigenvalues, the generalized real Schur form (\a S,\a T), which on exit replaces
// (\a A, \a B), and the left and right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

         \f[ (A^FA,B^FB) = ( (VSL^FL)*(S^FA)*(VSR^FR)^T, (VSL^FL)*T*(VSR^FR)^T ) \f]

// where \a FA, \a FB, \a FL, \a FR are transposition flags:
//  - \a FA == 1 if \a A is column-major and \a FA == T (transpose) if \a A is row-major;
//  - \a FB == 1 if \a B is column-major and \a FB == T (transpose) if \a B is row-major;
//  - \a FL == 1 if \a VSL is column-major and \a FL == T (transpose) if \a VSL is row-major;
//  - \a FR == 1 if \a VSR is column-major and \a FR == T (transpose) if \a VSR is row-major.
//
// A generalized eigenvalue for a pair of matrices (\a A,\a B) is a scalar \a w or a ratio
// \f$ alpha/beta = w \f$, such that \f$ A - w*B \f$ is singular. It is usually represented as
// the pair (\a alpha,\a beta), as there is a reasonable interpretation for \a beta=0 or both
// being zero. The complex eigenvalues are returned as numerators and denominators in the given
// vectors \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the size of the given matrices \a A and \a B don't match;
//  - ... the given matrix \a VSL is a fixed size matrix and the size doesn't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VSR is a fixed size matrix and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The first general matrix
   DynamicMatrix<double,rowMajor> B( 5UL, 5UL );  // The second general matrix
   // ... Initialization

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The vector for the eigenvalue numerators
   DynamicVector<double,columnVector> beta( 5UL );            // The vector for the eigenvalue denominators
   DynamicMatrix<double,columnVector> VSL( 5UL, 5UL );        // The matrix for the left Schur vectors
   DynamicMatrix<double,columnVector> VSR( 5UL, 5UL );        // The matrix for the right Schur vectors

   gges( A, B, VSL, alpha, beta, VSR );
   \endcode

// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix B
        , bool SO2      // Storage order of the matrix B
        , typename MT3  // Type of the matrix VSL
        , bool SO3      // Storage order of the matrix VSL
        , typename VT1  // Type of the vector alpha
        , bool TF1      // Transpose flag of the vector alpha
        , typename VT2  // Type of the vector beta
        , bool TF2      // Transpose flag of the vector beta
        , typename MT4  // Type of the matrix VSR
        , bool SO4 >    // Storage order of the matrix VSR
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                  DenseMatrix<MT3,SO3>& VSL, DenseVector<VT1,TF1>& alpha,
                  DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT4> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   if( (*A).rows() != (*B).rows() || (*A).columns() != (*B).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   resize( *VSL, N, N, false );
   resize( *alpha, N, false );
   resize( *beta, N, false );
   resize( *VSR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( *A, *B, *VSL, *alpha, *beta, *VSR, nullptr );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur factorization of the given pair of
//        dense general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The first matrix of the given pair of complex general matrices.
// \param B The second matrix of the given pair of complex general matrices.
// \param VSL The resulting complex matrix of left Schur vectors.
// \param alpha The resulting complex vector of eigenvalue numerators.
// \param beta The resulting complex vector of eigenvalue denominators.
// \param VSR The resulting complex matrix of right Schur vectors.
// \param select Logical function for the eigenvalue selection.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Vector or matrix cannot be resized.
// \exception std::runtime_error Schur factorization computation failed.
//
// This function computes for a pair of \a n-by-\a n real non-symmetric matrices (\a A,\a B), the
// generalized eigenvalues, the generalized real Schur form (\a S,\a T), which on exit replaces
// (\a A, \a B), and the left and right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

         \f[ (A^FA,B^FB) = ( (VSL^FL)*(S^FA)*(VSR^FR)^T, (VSL^FL)*T*(VSR^FR)^T ) \f]

// where \a FA, \a FB, \a FL, \a FR are transposition flags:
//  - \a FA == 1 if \a A is column-major and \a FA == T (transpose) if \a A is row-major;
//  - \a FB == 1 if \a B is column-major and \a FB == T (transpose) if \a B is row-major;
//  - \a FL == 1 if \a VSL is column-major and \a FL == T (transpose) if \a VSL is row-major;
//  - \a FR == 1 if \a VSR is column-major and \a FR == T (transpose) if \a VSR is row-major.
//
// A generalized eigenvalue for a pair of matrices (\a A,\a B) is a scalar \a w or a ratio
// \f$ alpha/beta = w \f$, such that \f$ A - w*B \f$ is singular. It is usually represented as
// the pair (\a alpha,\a beta), as there is a reasonable interpretation for \a beta=0 or both
// being zero. The complex eigenvalues are returned as numerators and denominators in the given
// vectors \a alpha, \a beta, which are resized to the correct size (if possible and necessary).
//
// It also orders the eigenvalues such that a selected cluster of eigenvalues appears in the
// leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The leading columns of \a VSL and \a VSR then form an orthonormal basis for the
// corresponding left and right eigenspaces (deflating subspaces). The according selection of
// eigenvalues is performed via the given \a select function. In case the given pair of matrices
// is a pair of real matrices the function is required to accept three pointers to real arguments,
// in case the given pair of matrices is a pair of complex matrices the function is required to
// accept two pointers to complex arguments.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the size of the given matrices \a A and \a B don't match;
//  - ... the given matrix \a VSL is a fixed size matrix and the size doesn't match;
//  - ... the given vector \a alpha is a fixed size vector and the size doesn't match;
//  - ... the given vector \a beta is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VSR is a fixed size matrix and the size doesn't match;
//  - ... the Schur factorization computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The first general matrix
   DynamicMatrix<double,rowMajor> B( 5UL, 5UL );  // The second general matrix
   // ... Initialization

   const auto select = []( T* alphar, T* alphai, T* beta ) -> blas_int_t {
      return *alphar > 0.0;
   };

   DynamicVector<complex<double>,columnVector> alpha( 5UL );  // The vector for the eigenvalue numerators
   DynamicVector<double,columnVector> beta( 5UL );            // The vector for the eigenvalue denominators
   DynamicMatrix<double,columnVector> VSL( 5UL, 5UL );        // The matrix for the left Schur vectors
   DynamicMatrix<double,columnVector> VSR( 5UL, 5UL );        // The matrix for the right Schur vectors

   gges( A, B, VSL, alpha, beta, VSR, select );
   \endcode

// For more information on the gges() functions (i.e. sgges(), dgges(), cgges(), and zgges())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1       // Type of the matrix A
        , bool SO1           // Storage order of the matrix A
        , typename MT2       // Type of the matrix B
        , bool SO2           // Storage order of the matrix B
        , typename MT3       // Type of the matrix VSL
        , bool SO3           // Storage order of the matrix VSL
        , typename VT1       // Type of the vector alpha
        , bool TF1           // Transpose flag of the vector alpha
        , typename VT2       // Type of the vector beta
        , bool TF2           // Transpose flag of the vector beta
        , typename MT4       // Type of the matrix VSR
        , bool SO4           // Storage order of the matrix VSR
        , typename Select >  // Type of the eigenvalue selector
inline void gges( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& B,
                  DenseMatrix<MT3,SO3>& VSL, DenseVector<VT1,TF1>& alpha,
                  DenseVector<VT2,TF2>& beta, DenseMatrix<MT4,SO4>& VSR,
                  Select select )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT4> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *B, N, N, false );
   resize( *VSL, N, N, false );
   resize( *alpha, N, false );
   resize( *beta, N, false );
   resize( *VSR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   gges_backend( *A, *B, *VSL, *alpha, *beta, *VSR, select );
}
//*************************************************************************************************

} // namespace blaze

#endif
