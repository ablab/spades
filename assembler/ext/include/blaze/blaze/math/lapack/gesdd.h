//=================================================================================================
/*!
//  \file blaze/math/lapack/gesdd.h
//  \brief Header file for the LAPACK singular value decomposition functions (gesdd)
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

#ifndef _BLAZE_MATH_LAPACK_GESDD_H_
#define _BLAZE_MATH_LAPACK_GESDD_H_


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
#include <blaze/math/lapack/clapack/gesdd.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK SVD FUNCTIONS (GESDD)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK SVD functions (gesdd) */
//@{
template< typename MT, bool SO, typename VT, bool TF >
void gesdd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s, char jobz );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
void gesdd( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V, char jobz );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3 >
void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
            DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobz );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param s The resulting vector of singular values.
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector s
        , bool TF >    // Transpose flag of the vector s
inline auto gesdd_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT> > >
{
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using ET = ElementType_t<MT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( 3*minimum + max( maximum, 7*minimum ) + 2 );
   const std::unique_ptr<ET[]> work( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( 'N', m, n, (*A).data(), lda, (*s).data(),
          nullptr, 1, nullptr, 1, work.get(), lwork, iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param s The resulting vector of singular values.
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector s
        , bool TF >    // Transpose flag of the vector s
inline auto gesdd_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT> > >
{
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT>;
   using BT = UnderlyingElement_t<CT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( 2*minimum + maximum + 2 );
   const std::unique_ptr<CT[]> work( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[7*minimum] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( 'N', m, n, (*A).data(), lda, (*s).data(),
          nullptr, 1, nullptr, 1, work.get(), lwork, rwork.get(), iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \return void
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix
// based on the LAPACK gesdd() functions, which use a divide-and-conquer strategy.
//
// The complete decomposition has the form

                           \f[ A = U \cdot S \cdot V, \f]

// where \c S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \c S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A.
//
// The resulting min(\a m,\a n) real and non-negative singular values are returned in descending
// order in the vector \a s, which is resized to the correct size (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the singular value decomposition fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s( 5UL );  // The vector for the singular values

   gesdd( A, s );
   \endcode

// For more information on the gesdd() functions (i.e. sgesdd(), dgesdd(), cgesdd(), and zgesdd())
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
        , typename VT  // Type of the vector s
        , bool TF >    // Transpose flag of the vector s
inline void gesdd( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s )
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

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   resize( *s, mindim, false );

   if( M == 0UL || N == 0 ) {
      return;
   }

   gesdd_backend( A, s );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF >       // Transpose flag of the vector s
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                           DenseVector<VT,TF>& s, char jobz )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'O' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || isSquare( *U ), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*U).rows() == (*A).rows(), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using ET = ElementType_t<MT1>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( ( jobz == 'O' )
                     ?( 3*minimum + max( maximum, 5*minimum*minimum + 4*maximum ) + 2 )
                     :( 3*minimum + max( maximum, 7*minimum ) + 2 ) );
   const std::unique_ptr<ET[]> work( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? (*U).data() : nullptr ), ( SO ? ldu : 1 ),
          ( SO ? nullptr : (*U).data() ), ( SO ? 1 : ldu ),
          work.get(), lwork, iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF >       // Transpose flag of the vector s
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                           DenseVector<VT,TF>& s, char jobz )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'O' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || isSquare( *U ), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*U).rows() == (*A).rows(), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( ( jobz == 'O' )
                     ?( 2*minimum*minimum + 2*minimum + maximum + 2 )
                     :( 2*minimum + maximum + 2 ) );
   const blas_int_t lrwork( ( jobz == 'O' )
                            ?( max( 5*minimum*minimum + 5*minimum,
                                    2*maximum*minimum + 2*minimum*minimum + minimum ) )
                            :( 7*minimum ) );
   const std::unique_ptr<CT[]> work( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[lrwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? (*U).data() : nullptr ), ( SO ? ldu : 1 ),
          ( SO ? nullptr : (*U).data() ), ( SO ? 1 : ldu ),
          work.get(), lwork, rwork.get(), iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::invalid_argument Invalid input matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid jobz argument provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix
// (\a m < \a n) based on the LAPACK gesdd() functions, which uses a divide-and-conquer strategy.
// Optionally, it computes the left or right singular vectors.
//
// The complete decomposition has the form

                           \f[ A = U \cdot S \cdot V, \f]

// where \c S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \c S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A.
//
// The complete decomposition has the form

                           \f[ A = U \cdot S \cdot V, \f]

// where \c S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \c S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A.
//
// The resulting min(\a m,\a n) real and non-negative singular values are returned in descending
// order in the vector \a s, which is resized to the correct size (if possible and necessary).
//
// The parameter \a jobz specifies the computation of the left and right singular vectors:
//
//   - \c 'O': If \a m < \a n, the \a m left singular vectors are returned in the columns of \a U
//             and the \a m right singular vectors are returned in the rows of \a A; \a U is either
//             resized (if possible) or expected to be a \a m-by-\a m matrix, the dimensions of
//             \a A are not adapted. Otherwise a \a std::invalid_argument exception is thrown.
//   - \c 'N': No left and right signular vectors are computed; \a U is not referenced.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... \a jobz == \c 'O' and \a m >= \a n;
//  - ... the given matrix \a U is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given \a jobz argument is neither \c 'O' nor \c 'N';
//  - ... the singular value decomposition fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U( 5UL, 5UL );  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s( 5UL );       // The vector for the singular values

   gesdd( A, U, s, 'O' );
   \endcode

// For more information on the gesdd() functions (i.e. sgesdd(), dgesdd(), cgesdd(), and zgesdd())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF >       // Transpose flag of the vector s
inline void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                   DenseVector<VT,TF>& s, char jobz )
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

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   if( jobz != 'O' && jobz != 'N' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid jobz argument provided" );
   }

   if( jobz == 'O' && M >= N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix provided" );
   }

   resize( *s, mindim, false );

   if( jobz == 'O' ) {
      resize( *U, M, M, false );
   }

   if( (*A).rows() == 0UL || (*A).columns() == 0UL ) {
      return;
   }

   gesdd_backend( *A, *U, *s, jobz );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT2 >  // Type of the matrix V
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s,
                           DenseMatrix<MT2,SO>& V, char jobz )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'O' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || isSquare( *V ), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*V).rows() == (*A).columns(), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using ET = ElementType_t<MT1>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( ( jobz == 'O' )
                     ?( 3*minimum + max( maximum, 5*minimum*minimum + 4*maximum + 2 ) )
                     :( 3*minimum + max( maximum, 7*minimum ) + 2 ) );
   const std::unique_ptr<ET[]> work( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? nullptr : (*V).data() ), ( SO ? 1 : ldv ),
          ( SO ? (*V).data() : nullptr ), ( SO ? ldv : 1 ),
          work.get(), lwork, iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT2 >  // Type of the matrix V
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s,
                           DenseMatrix<MT2,SO>& V, char jobz )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'O' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || isSquare( *V ), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*V).rows() == (*A).columns(), "Invalid matrix dimensions detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( ( jobz == 'O' )
                     ?( 2*minimum*minimum + 2*minimum + maximum + 2 )
                     :( 2*minimum + maximum + 2 ) );
   const blas_int_t lrwork( ( jobz == 'O' )
                            ?( max( 5*minimum*minimum + 5*minimum,
                                    2*maximum*minimum + 2*minimum*minimum + minimum ) )
                            :( 7*minimum ) );
   const std::unique_ptr<CT[]> work( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[lrwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? nullptr : (*V).data() ), ( SO ? 1 : ldv ),
          ( SO ? (*V).data() : nullptr ), ( SO ? ldv : 1 ),
          work.get(), lwork, rwork.get(), iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'O' or \c 'N').
// \return void
// \exception std::invalid_argument Invalid input matrix provided.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid jobz argument provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix
// (\a m >= \a n) based on the LAPACK gesdd() functions, which uses a divide-and-conquer strategy.
// Optionally, it computes the left or right singular vectors.
//
// The complete decomposition has the form

                           \f[ A = U \cdot S \cdot V, \f]

// where \c S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \c S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A.
//
// The resulting min(\a m,\a n) real and non-negative singular values are returned in descending
// order in the vector \a s, which is resized to the correct size (if possible and necessary).
//
// The parameter \a jobz specifies the computation of the left and right singular vectors:
//
//   - \c 'O': If \a m >= \a n, the \a n left singular vectors are returned in the columns of \a A
//             and the \a n right singular vectors are returned in the rows of \a V; \a V is either
//             resized (if possible) or expected to be a \a n-by-\a n matrix, the dimensions of
//             \a A are not adapted. Otherwise a \a std::invalid_argument exception is thrown.
//   - \c 'N': No left and right signular vectors are computed; \a V is not referenced.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... \a jobz == \c 'O' and \a m < \a n;
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the given \a jobz argument is neither \c 'O' nor \c 'N';
//  - ... the singular value decomposition fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s( 5UL );       // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V( 5UL, 8UL );  // The matrix for the right singular vectors

   gesdd( A, s, V, 'O' );
   \endcode

// For more information on the gesdd() functions (i.e. sgesdd(), dgesdd(), cgesdd(), and zgesdd())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF >       // Transpose flag of the vector s
inline void gesdd( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s,
                   DenseMatrix<MT2,SO>& V, char jobz )
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

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   if( jobz != 'O' && jobz != 'N' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid jobz argument provided" );
   }

   if( jobz == 'O' && M < N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix provided" );
   }

   resize( *s, mindim, false );

   if( jobz == 'O' ) {
      resize( *V, N, N, false );
   }

   if( (*A).rows() == 0UL || (*A).columns() == 0UL ) {
      return;
   }

   gesdd_backend( *A, *s, *V, jobz );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'A', \c 'S', or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT3 >  // Type of the matrix V
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                           DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobz )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'A' || jobz == 'S' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*U).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'A' || isSquare( *U ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'S' || (*U).columns() == min( (*A).rows(), (*A).columns() ), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'A' || isSquare( *V ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'S' || (*V).rows() == min( (*A).rows(), (*A).columns() ), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using ET = ElementType_t<MT1>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( 4*minimum*minimum + 6*minimum + maximum + 2 );
   const std::unique_ptr<ET[]> work( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? (*U).data() : (*V).data() ), ( SO ? ldu : ldv ),
          ( SO ? (*V).data() : (*U).data() ), ( SO ? ldv : ldu ),
          work.get(), lwork, iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesdd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'A', \c 'S', or \c 'N').
// \return void
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesdd() function.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT3 >  // Type of the matrix V
inline auto gesdd_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                           DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobz )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> > >
{
   BLAZE_INTERNAL_ASSERT( jobz == 'A' || jobz == 'S' || jobz == 'N', "Invalid jobz flag detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*U).rows() == (*A).rows(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'A' || isSquare( *U ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'S' || (*U).columns() == min( (*A).rows(), (*A).columns() ), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz == 'N' || (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'A' || isSquare( *V ), "Invalid non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( jobz != 'S' || (*V).rows() == min( (*A).rows(), (*A).columns() ), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t info( 0 );

   const blas_int_t minimum( min( m, n ) );
   const blas_int_t maximum( max( m, n ) );

   blas_int_t lwork( 4*minimum*minimum + 6*minimum + maximum + 2 );
   const blas_int_t lrwork( max( 5*minimum*minimum + 5*minimum,
                                 2*maximum*minimum + 2*minimum*minimum + minimum ) );
   const std::unique_ptr<CT[]> work( new CT[lwork] );
   const std::unique_ptr<BT[]> rwork( new BT[lrwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[8*minimum] );

   gesdd( jobz, m, n, (*A).data(), lda, (*s).data(),
          ( SO ? (*U).data() : (*V).data() ), ( SO ? ldu : ldv ),
          ( SO ? (*V).data() : (*U).data() ), ( SO ? ldv : ldu ),
          work.get(), lwork, rwork.get(), iwork.get(), &info );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param jobz Specifies the computation of the singular vectors (\c 'A', \c 'S', or \c 'N').
// \return void
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid jobz argument provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix
// based on the LAPACK gesdd() functions, which uses a divide-and-conquer strategy. Optionally,
// it computes the left or right singular vectors.
//
// The complete decomposition has the form

                           \f[ A = U \cdot S \cdot V, \f]

// where \c S is a \a m-by-\a n matrix, which is zero except for its min(\a m,\a n) diagonal
// elements, \a U is an \a m-by-\a m orthogonal matrix, and \a V is a \a n-by-\a n orthogonal
// matrix. The diagonal elements of \c S are the singular values of \a A, the first min(\a m,\a n)
// columns of \a U and rows of \a V are the left and right singular vectors of \a A.
//
// The resulting min(\a m,\a n) real and non-negative singular values are returned in descending
// order in the vector \a s, which is resized to the correct size (if possible and necessary).
//
// The parameter \a jobz specifies the computation of the left and right singular vectors:
//
//   - \c 'A': All \a m columns of \a U and all \a n rows of \a V are returned in \a U and
//             \a V; \a U is either resized (if possible) or expected to be a \a m-by-\a m matrix
//             and \a V is either resized (if possible) or expected to be a \a n-by-\a n matrix.
//   - \c 'S': The min(\a m,\a n) left singular vectors are returned in the columns of \a U;
//             \a U is either resized (if possible) or expected to be a \a m-by-min(\a m,\a n)
//             matrix. The min(\a m,\a n) right singular vectors are returned in the rows of \a V;
//             \a V is either resized (if possible) or expected to be a min(\a m,\a n)-by-\a n
//             matrix.
//   - \c 'N': No columns of \a U or rows of \a V are computed; both \a U and \a V are not
//             referenced.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given matrix \a U is a fixed size matrix and the dimensions don't match;
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the given \a jobz argument is neither \c 'A', \c 'S', nor \c 'N';
//  - ... the singular value decomposition fails.
//
// In all failure cases an exception is thrown.
//
// Examples:

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U( 5UL, 5UL );  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s( 5UL );       // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V( 8UL, 8UL );  // The matrix for the right singular vectors

   gesdd( A, U, s, V, 'A' );
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U( 5UL, 5UL );  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s( 5UL );       // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V( 5UL, 8UL );  // The matrix for the right singular vectors

   gesdd( A, U, s, V, 'S' );
   \endcode

// For more information on the gesdd() functions (i.e. sgesdd(), dgesdd(), cgesdd(), and zgesdd())
// see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename MT2    // Type of the matrix U
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT3 >  // Type of the matrix V
inline void gesdd( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                   DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, char jobz )
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

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   if( jobz != 'A' && jobz != 'S' && jobz != 'N' ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid jobz argument provided" );
   }

   resize( *s, mindim, false );

   if( jobz != 'N' ) {
      resize( *U, M, ( jobz == 'A' ? M : mindim ), false );
      resize( *V, ( jobz == 'A' ? N : mindim ), N, false );
   }

   if( (*A).rows() == 0UL || (*A).columns() == 0UL ) {
      return;
   }

   gesdd_backend( *A, *U, *s, *V, jobz );
}
//*************************************************************************************************

} // namespace blaze

#endif
