//=================================================================================================
/*!
//  \file blaze/math/lapack/geev.h
//  \brief Header file for the LAPACK general matrix eigenvalue functions (geev)
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

#ifndef _BLAZE_MATH_LAPACK_GEEV_H_
#define _BLAZE_MATH_LAPACK_GEEV_H_


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
#include <blaze/math/lapack/clapack/geev.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
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
//  LAPACK GENERAL MATRIX EIGENVALUE FUNCTIONS (GEEV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK general matrix eigenvalue functions (geev) */
//@{
template< typename MT, bool SO, typename VT, bool TF >
void geev( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename VT, bool TF >
void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w );

template< typename MT1, bool SO1, typename VT, bool TF, typename MT2, bool SO2 >
void geev( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& VR );

template< typename MT1, bool SO1, typename MT2, bool SO2, typename VT, bool TF, typename MT3, bool SO3 >
void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL,
           DenseVector<VT,TF>& w, DenseMatrix<MT3,SO3>& VR );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given real general matrix.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues of the given real
// dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector w
        , bool TF >    // Transpose flag of the vector w
inline auto geev_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT>;
   using BT = ElementType_t<MT>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 3*n + 2 );
   const std::unique_ptr<BT[]> wr  ( new BT[n] );
   const std::unique_ptr<BT[]> wi  ( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[lwork] );

   geev( 'N', 'N', n, (*A).data(), lda, wr.get(), wi.get(),
         nullptr, 1, nullptr, 1, work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }

   for( size_t i=0UL; i<(*A).rows(); ++i ) {
      (*w)[i] = CT( wr[i], wi[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given complex general matrix.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues of the given
// complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector w
        , bool TF >    // Transpose flag of the vector w
inline auto geev_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size() == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT>;
   using BT = UnderlyingElement_t<CT>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 2*n + 2 );
   const std::unique_ptr<CT[]> work ( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[2*n] );

   geev( 'N', 'N', n, (*A).data(), lda, (*w).data(),
         nullptr, 1, nullptr, 1, work.get(), lwork, rwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general matrix.
// \ingroup lapack_eigenvalue
//
// \param A The given general matrix.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function computes the eigenvalues of a non-symmetric, non-Hermitian \a n-by-\a n matrix
// based on the LAPACK geev() functions. The complex eigenvalues are returned in the given vector
// \a w, which is resized to the correct size (if possible and necessary). Please note that no
// order of eigenvalues can be assumed, except that complex conjugate pairs of eigenvalues appear
// consecutively with the eigenvalue having the positive imaginary part first.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the eigenvalue computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The general matrix A
   // ... Initialization

   DynamicVector<complex<double>,columnVector> w( 5UL );  // The vector for the complex eigenvalues

   geev( A, w );
   \endcode

// For more information on the geev() functions (i.e. sgeev(), dgeev(), cgeev(), and zgeev())
// see the LAPACK online documentation browser:
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
inline void geev( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& w )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *w, N, false );

   if( N == 0UL ) {
      return;
   }

   geev_backend( *A, *w );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given real general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and left
// eigenvectors of the given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF >     // Transpose flag of the vector w
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 4*n + 2 );
   const std::unique_ptr<BT[]> vl  ( new BT[n*n] );
   const std::unique_ptr<BT[]> wr  ( new BT[n] );
   const std::unique_ptr<BT[]> wi  ( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[lwork] );

   geev( ( SO1 ? 'V' : 'N' ), ( SO1 ? 'N' : 'V' ), n, (*A).data(), lda, wr.get(), wi.get(),
         ( SO1 ? vl.get() : nullptr ), ( SO1 ? n : 1 ),
         ( SO1 ? nullptr : vl.get() ), ( SO1 ? 1 : n ),
         work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }

   const size_t N( (*A).rows() );

   for( size_t j=0UL; j<N; ++j ) {
      (*w)[j] = CT( wr[j], wi[j] );
   }

   for( size_t j=0UL; j<N; ++j ) {
      if( j+1UL < N && equal( (*w)[j], conj( (*w)[j+1UL] ) ) ) {
         for( size_t i=0UL; i<N; ++i )
         {
            const size_t j1( SO1 ? j : j+1UL );
            const size_t j2( SO1 ? j+1UL : j );

            const BT vl1( vl[i+j*N] );
            const BT vl2( vl[i+(j+1UL)*N] );

            ( SO2 ? (*VL)(i,j1) : (*VL)(j1,i) ) = CT( vl1, ( SO2 ?  vl2 : -vl2 ) );
            ( SO2 ? (*VL)(i,j2) : (*VL)(j2,i) ) = CT( vl1, ( SO2 ? -vl2 :  vl2 ) );
         }

         ++j;
      }
      else {
         for( size_t i=0UL; i<N; ++i ) {
            ( SO2 ? (*VL)(i,j) : (*VL)(j,i) ) = CT( vl[i+j*N] );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given complex general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and left
// of the given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF >     // Transpose flag of the vector w
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldvl( numeric_cast<blas_int_t>( (*VL).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 2*n + 2 );
   const std::unique_ptr<CT[]> work ( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[2*n] );

   geev( ( SO1 ? 'V' : 'N' ), ( SO1 ? 'N' : 'V' ), n, (*A).data(), lda, (*w).data(),
         ( SO1 ? (*VL).data() : nullptr ), ( SO1 ? ldvl : 1 ),
         ( SO1 ? nullptr : (*VL).data() ), ( SO1 ? 1 : ldvl ),
         work.get(), lwork, rwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general matrix.
// \ingroup lapack_eigenvalue
//
// \param A The given general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function computes the eigenvalues of a non-symmetric, non-Hermitian \a n-by-\a n matrix
// based on the LAPACK geev() functions. Additionally, it computes the left eigenvectors. The
// left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// The complex eigenvalues are returned in the given vector \a w. The left eigenvectors are
// returned in the rows of \a VL in case \a VL is a row-major matrix and in the columns of
// \a VL in case \a VL is a column-major matrix. Both \a w and \a VL are resized to the correct
// dimensions (if possible and necessary). Please note that no order of eigenvalues can be
// assumed, except that complex conjugate pairs of eigenvalues appear consecutively with the
// eigenvalue having the positive imaginary part first. The computed eigenvectors are normalized
// to have Euclidean norm equal to 1 and largest component real.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given matrix \a VL is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the eigenvalue computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<complex<double>,rowMajor> VL( 5UL, 5UL );  // The matrix for the left eigenvectors
   DynamicVector<complex<double>,columnVector> w( 5UL );    // The vector for the complex eigenvalues

   geev( A, VL, w );
   \endcode

// For more information on the geev() functions (i.e. sgeev(), dgeev(), cgeev(), and zgeev())
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
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF >     // Transpose flag of the vector w
inline void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL, DenseVector<VT,TF>& w )
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
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *w, N, false );
   resize( *VL, N, N, false );

   if( N == 0UL ) {
      return;
   }

   geev_backend( *A, *VL, *w );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given real general matrix.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and right
// eigenvectors of the given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT2  // Type of the matrix VR
        , bool SO2 >    // Storage order of the matrix VR
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& VR )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 4*n + 2 );
   const std::unique_ptr<BT[]> vr  ( new BT[n*n] );
   const std::unique_ptr<BT[]> wr  ( new BT[n] );
   const std::unique_ptr<BT[]> wi  ( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[lwork] );

   geev( ( SO1 ? 'N' : 'V' ), ( SO1 ? 'V' : 'N' ), n, (*A).data(), lda, wr.get(), wi.get(),
         ( SO1 ? nullptr : vr.get() ), ( SO1 ? 1 : n ),
         ( SO1 ? vr.get() : nullptr ), ( SO1 ? n : 1 ),
         work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }

   const size_t N( (*A).rows() );

   for( size_t j=0UL; j<N; ++j ) {
      (*w)[j] = CT( wr[j], wi[j] );
   }

   for( size_t j=0UL; j<N; ++j ) {
      if( j+1UL < N && equal( (*w)[j], conj( (*w)[j+1UL] ) ) ) {
         for( size_t i=0UL; i<N; ++i )
         {
            const size_t j1( SO1 ? j : j+1UL );
            const size_t j2( SO1 ? j+1UL : j );

            const BT vr1( vr[i+j*N] );
            const BT vr2( vr[i+(j+1UL)*N] );

            ( SO2 ? (*VR)(i,j1) : (*VR)(j1,i) ) = CT( vr1, ( SO2 ?  vr2 : -vr2 ) );
            ( SO2 ? (*VR)(i,j2) : (*VR)(j2,i) ) = CT( vr1, ( SO2 ? -vr2 :  vr2 ) );
         }

         ++j;
      }
      else {
         for( size_t i=0UL; i<N; ++i ) {
            ( SO2 ? (*VR)(i,j) : (*VR)(j,i) ) = CT( vr[i+j*N] );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given complex general matrix.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and right
// eigenvectors of the given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT2  // Type of the matrix VR
        , bool SO2 >    // Storage order of the matrix VR
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& VR )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldvr( numeric_cast<blas_int_t>( (*VR).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 2*n + 2 );
   const std::unique_ptr<CT[]> work ( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[2*n] );

   geev( ( SO1 ? 'N' : 'V' ), ( SO1 ? 'V' : 'N' ), n, (*A).data(), lda, (*w).data(),
         ( SO1 ? nullptr : (*VR).data() ), ( SO1 ? 1 : ldvr ),
         ( SO1 ? (*VR).data() : nullptr ), ( SO1 ? ldvr : 1 ),
         work.get(), lwork, rwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general matrix.
// \ingroup lapack_eigenvalue
//
// \param A The given general matrix.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function computes the eigenvalues of a non-symmetric, non-Hermitian \a n-by-\a n matrix
// based on the LAPACK geev() functions. Additionally, it computes the right eigenvectors. The
// right eigenvector \f$v[j]\f$ of \a A satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue.
//
// The complex eigenvalues are returned in the given vector \a w. The right eigenvectors are
// returned in the rows of \a VR in case \a VR is a row-major matrix and in the columns of
// \a VR in case \a VR is a column-major matrix. Both \a w and \a VL are resized to the correct
// dimensions (if possible and necessary). Please note that no order of eigenvalues can be
// assumed, except that complex conjugate pairs of eigenvalues appear consecutively with the
// eigenvalue having the positive imaginary part first. The computed eigenvectors are normalized
// to have Euclidean norm equal to 1 and largest component real.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VR is a fixed size matrix and the dimensions don't match;
//  - ... the eigenvalue computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The general matrix A
   // ... Initialization

   DynamicVector<complex<double>,columnVector> w( 5UL );    // The vector for the complex eigenvalues
   DynamicMatrix<complex<double>,rowMajor> VR( 5UL, 5UL );  // The matrix for the right eigenvectors

   geev( A, w, VR );
   \endcode

// For more information on the geev() functions (i.e. sgeev(), dgeev(), cgeev(), and zgeev())
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
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT2  // Type of the matrix VR
        , bool SO2 >    // Storage order of the matrix VR
inline void geev( DenseMatrix<MT1,SO1>& A, DenseVector<VT,TF>& w, DenseMatrix<MT2,SO2>& VR )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<MT2> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *w, N, false );
   resize( *VR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   geev_backend( *A, *w, *VR );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for real general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given real general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and left and right
// eigenvectors of the given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT3  // Type of the matrix VR
        , bool SO3 >    // Storage order of the matrix VR
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL,
                          DenseVector<VT,TF>& w, DenseMatrix<MT3,SO3>& VR )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<VT>;
   using BT = ElementType_t<MT1>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 4*n + 2 );
   const std::unique_ptr<BT[]> vl  ( new BT[n*n] );
   const std::unique_ptr<BT[]> vr  ( new BT[n*n] );
   const std::unique_ptr<BT[]> wr  ( new BT[n] );
   const std::unique_ptr<BT[]> wi  ( new BT[n] );
   const std::unique_ptr<BT[]> work( new BT[lwork] );

   geev( 'V', 'V', n, (*A).data(), lda, wr.get(), wi.get(),
         ( SO1 ? vl.get() : vr.get() ), n, ( SO1 ? vr.get() : vl.get() ), n,
         work.get(), lwork, &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }

   const size_t N( (*A).rows() );

   for( size_t j=0UL; j<N; ++j ) {
      (*w)[j] = CT( wr[j], wi[j] );
   }

   for( size_t j=0UL; j<N; ++j ) {
      if( j+1UL < N && equal( (*w)[j], conj( (*w)[j+1UL] ) ) ) {
         for( size_t i=0UL; i<N; ++i )
         {
            const size_t j1( SO1 ? j : j+1UL );
            const size_t j2( SO1 ? j+1UL : j );

            const BT vl1( vl[i+j*N] );
            const BT vl2( vl[i+(j+1UL)*N] );
            const BT vr1( vr[i+j*N] );
            const BT vr2( vr[i+(j+1UL)*N] );

            ( SO2 ? (*VL)(i,j1) : (*VL)(j1,i) ) = CT( vl1, ( SO2 ?  vl2 : -vl2 ) );
            ( SO2 ? (*VL)(i,j2) : (*VL)(j2,i) ) = CT( vl1, ( SO2 ? -vl2 :  vl2 ) );
            ( SO3 ? (*VR)(i,j1) : (*VR)(j1,i) ) = CT( vr1, ( SO3 ?  vr2 : -vr2 ) );
            ( SO3 ? (*VR)(i,j2) : (*VR)(j2,i) ) = CT( vr1, ( SO3 ? -vr2 :  vr2 ) );
         }

         ++j;
      }
      else {
         for( size_t i=0UL; i<N; ++i ) {
            ( SO2 ? (*VL)(i,j) : (*VL)(j,i) ) = CT( vl[i+j*N] );
            ( SO3 ? (*VR)(i,j) : (*VR)(j,i) ) = CT( vr[i+j*N] );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK geev kernel for complex general matrices.
// \ingroup lapack_eigenvalue
//
// \param A The given complex general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function is the backend implementation for computing the eigenvalues and left and right
// of the given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// geev() function.
*/
template< typename MT1  // Type of the matrix A
        , bool SO1      // Storage order of the matrix A
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT3  // Type of the matrix VR
        , bool SO3 >    // Storage order of the matrix VR
inline auto geev_backend( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL,
                          DenseVector<VT,TF>& w, DenseMatrix<MT3,SO3>& VR )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( isSquare( *A ) , "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VL ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( isSquare( *VR ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*VL).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*VR).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*w).size()  == (*A).rows(), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( CT );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( BT );

   blas_int_t n   ( numeric_cast<blas_int_t>( (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldvl( numeric_cast<blas_int_t>( (*VL).spacing() ) );
   blas_int_t ldvr( numeric_cast<blas_int_t>( (*VR).spacing() ) );
   blas_int_t info( 0 );

   blas_int_t lwork( 2*n + 2 );
   const std::unique_ptr<CT[]> work ( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[2*n] );

   geev( 'V', 'V', n, (*A).data(), lda, (*w).data(),
         ( SO1 ? (*VL).data() : (*VR).data() ), ( SO1 ? ldvl : ldvr ),
         ( SO1 ? (*VR).data() : (*VL).data() ), ( SO1 ? ldvr : ldvl ),
         work.get(), lwork, rwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for eigenvalue computation" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Eigenvalue computation failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general matrix.
// \ingroup lapack_eigenvalue
//
// \param A The given general matrix.
// \param VL The resulting matrix of left eigenvectors.
// \param w The resulting vector of complex eigenvalues.
// \param VR The resulting matrix of right eigenvectors.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Eigenvalue computation failed.
//
// This function computes the eigenvalues of a non-symmetric, non-Hermitian \a n-by-\a n matrix
// based on the LAPACK geev() functions. Additionally, it computes the left and right eigenvectors.
// The right eigenvector \f$v[j]\f$ of \a A satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// The complex eigenvalues are returned in the given vector \a w. The left eigenvectors are
// returned in the rows of \a VL in case \a VL is a row-major matrix and in the columns of \a VL
// in case \a VL is a column-major matrix. The right eigenvectors are returned in the rows of
// \a VR in case \a VR is a row-major matrix and in the columns of \a VR in case \a VR is a
// column-major matrix. \a w, \a VL, and \a VR are resized to the correct dimensions (if
// possible and necessary). Please note that no order of eigenvalues can be assumed, except
// that complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having
// the positive imaginary part first. The computed eigenvectors are normalized to have Euclidean
// norm equal to 1 and largest component real.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a A is not a square matrix;
//  - ... the given matrix \a VL is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a w is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a VR is a fixed size matrix and the dimensions don't match;
//  - ... the eigenvalue computation fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 5UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<complex<double>,rowMajor> VL( 5UL, 5UL );  // The matrix for the left eigenvectors
   DynamicVector<complex<double>,columnVector> w( 5UL );    // The vector for the complex eigenvalues
   DynamicMatrix<complex<double>,rowMajor> VR( 5UL, 5UL );  // The matrix for the right eigenvectors

   geev( A, VL, w, VR );
   \endcode

// For more information on the geev() functions (i.e. sgeev(), dgeev(), cgeev(), and zgeev())
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
        , typename MT2  // Type of the matrix VL
        , bool SO2      // Storage order of the matrix VL
        , typename VT   // Type of the vector w
        , bool TF       // Transpose flag of the vector w
        , typename MT3  // Type of the matrix VR
        , bool SO3 >    // Storage order of the matrix VR
inline void geev( DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& VL,
                  DenseVector<VT,TF>& w, DenseMatrix<MT3,SO3>& VR )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );
   BLAZE_CONSTRAINT_MUST_BE_COMPLEX_TYPE( ElementType_t<MT3> );

   const size_t N( (*A).rows() );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   resize( *w, N, false );
   resize( *VL, N, N, false );
   resize( *VR, N, N, false );

   if( N == 0UL ) {
      return;
   }

   geev_backend( *A, *VL, *w, *VR );
}
//*************************************************************************************************

} // namespace blaze

#endif
