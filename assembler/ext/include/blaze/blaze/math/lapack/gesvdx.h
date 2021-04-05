//=================================================================================================
/*!
//  \file blaze/math/lapack/gesvdx.h
//  \brief Header file for the LAPACK singular value decomposition functions (gesvdx)
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

#ifndef _BLAZE_MATH_LAPACK_GESVDX_H_
#define _BLAZE_MATH_LAPACK_GESVDX_H_


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
#include <blaze/math/lapack/clapack/gesvdx.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/math/views/Subvector.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK SVD FUNCTIONS (GESVDX)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK SVD functions (gesvdx) */
//@{
template< typename MT, bool SO, typename VT, bool TF >
size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s );

template< typename MT, bool SO, typename VT, bool TF, typename ST >
size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, ST low, ST upp );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename ST >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
               DenseVector<VT,TF>& s, ST low, ST upp );

template< typename MT1, bool SO, typename VT, bool TF, typename MT2 >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V );

template< typename MT1, bool SO, typename VT, bool TF, typename MT2, typename ST >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s,
               DenseMatrix<MT2,SO>& V, ST low, ST upp );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3 >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
               DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V );

template< typename MT1, bool SO, typename MT2, typename VT, bool TF, typename MT3, typename ST >
size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
               DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, ST low, ST upp );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param s The resulting vector of singular values.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using ET = ElementType_t<MT>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   ET* sptr( (*s).data() );
   std::unique_ptr<ET[]> stmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new ET[2UL*mindim] );
      sptr = stmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<ET[]>  work ( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( 'N', 'N', range, m, n, (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           nullptr, 1, nullptr, 1, work.get(), lwork, iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired ) {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }
   }

   return num;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param s The resulting vector of singular values.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using CT = ElementType_t<MT>;
   using BT = UnderlyingElement_t<CT>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   BT* sptr( (*s).data() );
   std::unique_ptr<BT[]> stmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new BT[2UL*mindim] );
      sptr = stmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<CT[]>  work ( new CT[lwork] );
   const std::unique_ptr<BT[]>  rwork( new BT[17*minimum*minimum] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( 'N', 'N', range, m, n, (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           nullptr, 1, nullptr, 1, work.get(), lwork, rwork.get(), iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired ) {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }
   }

   return num;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes all min(\a m,\a n) singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The resulting real and non-negative
// singular values are returned in descending order in the vector \a s, which is resized to the
// correct size (if possible and necessary).
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

   gesvdx( A, s );
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
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
inline size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s )
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

   using ET = ElementType_t<MT>;
   using UT = UnderlyingElement_t<ET>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   resize( *s, mindim, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   return gesvdx_backend( *A, *s, 'A', UT(), UT(), 0, 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The number of singular values to be
// computed is specified by the lower bound \a low and the upper bound \a upp, which either form
// an integral or a floating point range.
//
// In case \a low and \a upp are of integral type, the function computes all singular values in
// the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given scalar values don't form a proper range;
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

   gesvdx( A, s, 1.0, 2.0 );  // Computes all eigenvalues in the interval (1..2]
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s( 3UL );  // The vector for the singular values

   gesvdx( A, s, 0, 2 );  // Computes the first three eigenvalues
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline size_t gesvdx( DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, ST low, ST upp )
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

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   if( IsFloatingPoint_v<ST> && low >= upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid value range provided" );
   }

   if( !IsFloatingPoint_v<ST> && low > upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );
   const size_t expected( IsFloatingPoint_v<ST> ? mindim : size_t( upp - low ) + 1UL );

   if( !IsFloatingPoint_v<ST> && expected > mindim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   resize( *s, expected, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   const char       range( IsFloatingPoint_v<ST> ? 'V' : 'I' );
   const ST         vl   ( IsFloatingPoint_v<ST> ? low : ST() );
   const ST         vu   ( IsFloatingPoint_v<ST> ? upp : ST() );
   const blas_int_t il   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( low ) );
   const blas_int_t iu   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( upp ) );

   const size_t actual( gesvdx_backend( *A, *s, range, vl, vu, il, iu ) );

   if( IsResizable_v<VT> ) {
      resize( *s, actual, true );
   }

   return actual;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).rows() == (*A).rows()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).columns() == (*s).size(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using ET = ElementType_t<MT1>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   ET* sptr( (*s).data() );
   ET* uptr( (*U).data() );
   std::unique_ptr<ET[]> stmp;
   std::unique_ptr<ET[]> utmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new ET[2UL*mindim] );
      utmp.reset( new ET[M*mindim] );
      sptr = stmp.get();
      uptr = utmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<ET[]>  work ( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( ( SO ? 'V' : 'N' ), ( SO ? 'N' : 'V' ), range, m, n,
           (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? uptr : nullptr ), ( tmpRequired ? m : ( SO ? ldu : 1 ) ),
           ( SO ? nullptr : uptr ), ( tmpRequired ? mindim : ( SO ? 1 : ldu ) ),
           work.get(), lwork, iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*U).columns(); ++j ) {
            for( size_t i=0UL; i<(*U).rows(); ++i ) {
               (*U)(i,j) = utmp[i+j*M];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*U).rows(); ++i ) {
            for( size_t j=0UL; j<(*U).columns(); ++j ) {
               (*U)(i,j) = utmp[i*mindim+j];
            }
         }
      }
   }

   return num;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).rows() == (*A).rows()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).columns() == (*s).size(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   BT* sptr( (*s).data() );
   CT* uptr( (*U).data() );
   std::unique_ptr<BT[]> stmp;
   std::unique_ptr<CT[]> utmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new BT[2UL*mindim] );
      utmp.reset( new CT[M*mindim] );
      sptr = stmp.get();
      uptr = utmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<CT[]>  work ( new CT[lwork] );
   const std::unique_ptr<BT[]>  rwork( new BT[17*minimum*minimum] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( ( SO ? 'V' : 'N' ), ( SO ? 'N' : 'V' ), range, m, n,
           (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? uptr : nullptr ), ( tmpRequired ? m : ( SO ? ldu : 1 ) ),
           ( SO ? nullptr : uptr ), ( tmpRequired ? mindim : ( SO ? 1 : ldu ) ),
           work.get(), lwork, rwork.get(), iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*U).columns(); ++j ) {
            for( size_t i=0UL; i<(*U).rows(); ++i ) {
               (*U)(i,j) = utmp[i+j*M];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*U).rows(); ++i ) {
            for( size_t j=0UL; j<(*U).columns(); ++j ) {
               (*U)(i,j) = utmp[i*mindim+j];
            }
         }
      }
   }

   return num;
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
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes all min(\a m,\a n) singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The resulting real and non-negative
// singular values are stored in descending order in the given vector \a s, the resulting left
// singular vectors are stored in the given matrix \a U. Both \a s and \a U are resized to the
// correct dimensions (if possible and necessary).
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

   gesvdx( A, U, s );
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
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
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s )
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

   using ET = ElementType_t<MT1>;
   using UT = UnderlyingElement_t<ET>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   resize( *s, mindim, false );
   resize( *U, M, mindim, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   return gesvdx_backend( *A, *U, *s, 'A', UT(), UT(), 0, 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The number of singular values to be
// computed is specified by the lower bound \a low and the upper bound \a upp, which either form
// an integral or a floating point range.
//
// In case \a low and \a upp are of integral type, the function computes all singular values in
// the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector. The resulting left singular vectors are stored
// in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-\a num matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector. The resulting left singular vectors are
// stored in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-min(\a m,\a n) matrix.
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
//  - ... the given scalar values don't form a proper range;
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

   gesvdx( A, U, s, 1.0, 2.0 );  // Computes all eigenvalues in the interval (1..2]
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U( 5UL, 3UL );  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s( 3UL );       // The vector for the singular values

   gesvdx( A, U, s, 0, 2 );  // Computes the first three eigenvalues
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                      DenseVector<VT,TF>& s, ST low, ST upp )
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

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   if( IsFloatingPoint_v<ST> && low >= upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid value range provided" );
   }

   if( !IsFloatingPoint_v<ST> && low > upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );
   const size_t expected( IsFloatingPoint_v<ST> ? mindim : size_t( upp - low ) + 1UL );

   if( !IsFloatingPoint_v<ST> && expected > mindim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   resize( *s, expected, false );
   resize( *U, M, expected, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   const char       range( IsFloatingPoint_v<ST> ? 'V' : 'I' );
   const ST         vl   ( IsFloatingPoint_v<ST> ? low : ST() );
   const ST         vu   ( IsFloatingPoint_v<ST> ? upp : ST() );
   const blas_int_t il   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( low ) );
   const blas_int_t iu   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( upp ) );

   const size_t actual( gesvdx_backend( *A, *U, *s, range, vl, vu, il, iu ) );

   if( IsResizable_v<VT> ) {
      resize( *s, actual, true );
   }

   if( IsResizable_v<MT2> ) {
      resize( *U, M, actual, true );
   }

   return actual;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT2   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).rows()    == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using ET = ElementType_t<MT1>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   ET* sptr( (*s).data() );
   ET* vptr( (*V).data() );
   std::unique_ptr<ET[]> stmp;
   std::unique_ptr<ET[]> vtmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new ET[2UL*mindim] );
      vtmp.reset( new ET[mindim*N] );
      sptr = stmp.get();
      vptr = vtmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<ET[]>  work ( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( ( SO ? 'N' : 'V' ), ( SO ? 'V' : 'N' ), range, m, n,
           (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? nullptr : vptr ), ( tmpRequired ? m : ( SO ? 1 : ldv ) ),
           ( SO ? vptr : nullptr ), ( tmpRequired ? mindim : ( SO ? ldv : 1 ) ),
           work.get(), lwork, iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*V).columns(); ++j ) {
            for( size_t i=0UL; i<(*V).rows(); ++i ) {
               (*V)(i,j) = vtmp[i+j*mindim];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*V).rows(); ++i ) {
            for( size_t j=0UL; j<(*V).columns(); ++j ) {
               (*V)(i,j) = vtmp[i*N+j];
            }
         }
      }
   }

   return num;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT2   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V,
                            char range, ST vl, ST vu, blas_int_t il, blas_int_t iu )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).rows()    == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   BT* sptr( (*s).data() );
   CT* vptr( (*V).data() );
   std::unique_ptr<BT[]> stmp;
   std::unique_ptr<CT[]> vtmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new BT[2UL*mindim] );
      vtmp.reset( new CT[mindim*N] );
      sptr = stmp.get();
      vptr = vtmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<CT[]>  work ( new CT[lwork] );
   const std::unique_ptr<BT[]>  rwork( new BT[17*minimum*minimum] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( ( SO ? 'N' : 'V' ), ( SO ? 'V' : 'N' ), range, m, n,
           (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? nullptr : vptr ), ( tmpRequired ? m : ( SO ? 1 : ldv ) ),
           ( SO ? vptr : nullptr ), ( tmpRequired ? mindim : ( SO ? ldv : 1 ) ),
           work.get(), lwork, rwork.get(), iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*V).columns(); ++j ) {
            for( size_t i=0UL; i<(*V).rows(); ++i ) {
               (*V)(i,j) = vtmp[i+j*mindim];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*V).rows(); ++i ) {
            for( size_t j=0UL; j<(*V).columns(); ++j ) {
               (*V)(i,j) = vtmp[i*N+j];
            }
         }
      }
   }

   return num;
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
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes all min(\a m,\a n) singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The resulting real and non-negative
// singular values are stored in descending order in the given vector \a s. The resulting right
// singular vectors are stored in the given min(\a m,\a n)-by-\a n matrix \a V. Both \a s and
// \a V are resized to the correct dimensions (if possible and necessary).
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
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

   gesvdx( A, s, V );
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT2 >  // Type of the matrix V
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s, DenseMatrix<MT2,SO>& V )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   using ET = ElementType_t<MT1>;
   using UT = UnderlyingElement_t<ET>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   resize( *s, mindim, false );
   resize( *V, mindim, N, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   return gesvdx_backend( *A, *s, *V, 'A', UT(), UT(), 0, 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The number of singular values to be
// computed is specified by the lower bound \a low and the upper bound \a upp, which either form
// an integral or a floating point range.
//
// In case \a low and \a upp are of integral type, the function computes all singular values in
// the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector. The resulting right singular vectors are stored
// in the given matrix \a V, which is either resized (if possible) or expected to be a
// \a num-by-\a n matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector. The resulting right singular vectors are
// stored in the given matrix \a V, which is either resized (if possible) or expected to be a
// min(\a m,\a n)-by-\a n matrix.
//
// Note that this function can only be used for general, non-adapted matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with any adapted matrix or matrices of any other element type results in a compile
// time error!
//
// The function fails if ...
//
//  - ... the given vector \a s is a fixed size vector and the size doesn't match;
//  - ... the given matrix \a V is a fixed size matrix and the dimensions don't match;
//  - ... the given scalar values don't form a proper range;
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

   gesvdx( A, s, V, 1.0, 2.0 );  // Computes all eigenvalues in the interval (1..2]
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s( 3UL );       // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V( 3UL, 8UL );  // The matrix for the right singular vectors

   gesvdx( A, s, V, 0, 2 );  // Computes the first three eigenvalues
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT2   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseVector<VT,TF>& s,
                      DenseMatrix<MT2,SO>& V, ST low, ST upp )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   if( IsFloatingPoint_v<ST> && low >= upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid value range provided" );
   }

   if( !IsFloatingPoint_v<ST> && low > upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );
   const size_t expected( IsFloatingPoint_v<ST> ? mindim : size_t( upp - low ) + 1UL );

   if( !IsFloatingPoint_v<ST> && expected > mindim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   resize( *s, expected, false );
   resize( *V, expected, N, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   const char       range( IsFloatingPoint_v<ST> ? 'V' : 'I' );
   const ST         vl   ( IsFloatingPoint_v<ST> ? low : ST() );
   const ST         vu   ( IsFloatingPoint_v<ST> ? upp : ST() );
   const blas_int_t il   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( low ) );
   const blas_int_t iu   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( upp ) );

   const size_t actual( gesvdx_backend( *A, *s, *V, range, vl, vu, il, iu ) );

   if( IsResizable_v<VT> ) {
      resize( *s, actual, true );
   }

   if( IsResizable_v<MT2> ) {
      resize( *V, actual, N, true );
   }

   return actual;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for real general matrices.
// \ingroup lapack_singular_value
//
// \param A The given real general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given real dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT3   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s,
                            DenseMatrix<MT3,SO>& V, char range, ST vl, ST vu,
                            blas_int_t il, blas_int_t iu )
   -> DisableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).rows()    == (*A).rows()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).columns() == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).rows()    == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using ET = ElementType_t<MT1>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   ET* sptr( (*s).data() );
   ET* uptr( (*U).data() );
   ET* vptr( (*V).data() );
   std::unique_ptr<ET[]> stmp;
   std::unique_ptr<ET[]> utmp;
   std::unique_ptr<ET[]> vtmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new ET[2UL*mindim] );
      utmp.reset( new ET[M*mindim] );
      vtmp.reset( new ET[mindim*N] );
      sptr = stmp.get();
      uptr = utmp.get();
      vptr = vtmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<ET[]>  work ( new ET[lwork] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( 'V', 'V', range, m, n, (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? uptr : vptr ), ( tmpRequired ? m : ( SO ? ldu : ldv ) ),
           ( SO ? vptr : uptr ), ( tmpRequired ? mindim : ( SO ? ldv : ldu ) ),
           work.get(), lwork, iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*U).columns(); ++j ) {
            for( size_t i=0UL; i<(*U).rows(); ++i ) {
               (*U)(i,j) = utmp[i+j*M];
            }
         }

         for( size_t j=0UL; j<(*V).columns(); ++j ) {
            for( size_t i=0UL; i<(*V).rows(); ++i ) {
               (*V)(i,j) = vtmp[i+j*mindim];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*U).rows(); ++i ) {
            for( size_t j=0UL; j<(*U).columns(); ++j ) {
               (*U)(i,j) = utmp[i*mindim+j];
            }
         }

         for( size_t i=0UL; i<(*V).rows(); ++i ) {
            for( size_t j=0UL; j<(*V).columns(); ++j ) {
               (*V)(i,j) = vtmp[i*N+j];
            }
         }
      }
   }

   return num;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the LAPACK gesvd kernel for complex general matrices.
// \ingroup lapack_singular_value
//
// \param A The given complex general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param range Specifies the range of singular values to find (\c 'A', \c 'V', or \c 'I').
// \param vl The lower bound of the interval to be searched for singular values (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for singular values (\a vl < \a vu).
// \param il The index of the smallest singular value to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest singular value to be returned (0 <= \a il <= \a iu).
// \return The total number of singular values found.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function is the backend implementation for the singular value decomposition of the
// given complex dense general matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the dispatch to
// the correct LAPACK function. Calling this function explicitly might result in erroneous
// results and/or in compilation errors. Instead of using this function use the according
// gesvdx() function.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT3   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline auto gesvdx_backend( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U, DenseVector<VT,TF>& s,
                            DenseMatrix<MT3,SO>& V, char range, ST vl, ST vu,
                            blas_int_t il, blas_int_t iu )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT1> >, size_t >
{
   BLAZE_INTERNAL_ASSERT( range == 'A' || range == 'V' || range == 'I', "Invalid range flag detected" );
   BLAZE_INTERNAL_ASSERT( range != 'A' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'V' || (*s).size() == min( (*A).rows(), (*A).columns() ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( range != 'I' || (*s).size() == size_t( iu-il+1 ), "Invalid vector dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).rows()    == (*A).rows()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*U).columns() == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).rows()    == (*s).size()   , "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( (*V).columns() == (*A).columns(), "Invalid matrix dimension detected" );
   BLAZE_INTERNAL_ASSERT( vl <= vu, "Invalid floating point range detected" );
   BLAZE_INTERNAL_ASSERT( il <= iu, "Invalid integral range detected" );

   using CT = ElementType_t<MT1>;
   using BT = UnderlyingElement_t<CT>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? M : N ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? N : M ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t ldu ( numeric_cast<blas_int_t>( (*U).spacing() ) );
   blas_int_t ldv ( numeric_cast<blas_int_t>( (*V).spacing() ) );
   blas_int_t ns  ( 0 );
   blas_int_t info( 0 );

   BT* sptr( (*s).data() );
   CT* uptr( (*U).data() );
   CT* vptr( (*V).data() );
   std::unique_ptr<BT[]> stmp;
   std::unique_ptr<CT[]> utmp;
   std::unique_ptr<CT[]> vtmp;

   const bool tmpRequired( (*s).size() < mindim );

   if( tmpRequired ) {
      stmp.reset( new BT[2UL*mindim] );
      utmp.reset( new CT[M*mindim] );
      vtmp.reset( new CT[mindim*N] );
      sptr = stmp.get();
      uptr = utmp.get();
      vptr = vtmp.get();
   }

   const blas_int_t minimum( min( m, n ) );

   blas_int_t lwork( minimum*( minimum*3 + 20 ) + 2 );
   const std::unique_ptr<CT[]>  work ( new CT[lwork] );
   const std::unique_ptr<BT[]>  rwork( new BT[17*minimum*minimum] );
   const std::unique_ptr<blas_int_t[]> iwork( new blas_int_t[12*minimum] );

   gesvdx( 'V', 'V', range, m, n, (*A).data(), lda, vl, vu, il, iu, &ns, sptr,
           ( SO ? uptr : vptr ), ( tmpRequired ? m : ( SO ? ldu : ldv ) ),
           ( SO ? vptr : uptr ), ( tmpRequired ? mindim : ( SO ? ldv : ldu ) ),
           work.get(), lwork, rwork.get(), iwork.get(), &info );

   const size_t num( numeric_cast<size_t>( ns ) );

   BLAZE_INTERNAL_ASSERT( info >= 0, "Invalid argument for singular value decomposition" );
   BLAZE_INTERNAL_ASSERT( num <= (*s).size(), "Invalid number of singular values detected" );

   if( info > 0 ) {
      BLAZE_THROW_LAPACK_ERROR( "Singular value decomposition failed" );
   }

   if( tmpRequired )
   {
      for( size_t i=0UL; i<(*s).size(); ++i ) {
         (*s)[i] = stmp[i];
      }

      if( SO )
      {
         for( size_t j=0UL; j<(*U).columns(); ++j ) {
            for( size_t i=0UL; i<(*U).rows(); ++i ) {
               (*U)(i,j) = utmp[i+j*M];
            }
         }

         for( size_t j=0UL; j<(*V).columns(); ++j ) {
            for( size_t i=0UL; i<(*V).rows(); ++i ) {
               (*V)(i,j) = vtmp[i+j*mindim];
            }
         }
      }
      else
      {
         for( size_t i=0UL; i<(*U).rows(); ++i ) {
            for( size_t j=0UL; j<(*U).columns(); ++j ) {
               (*U)(i,j) = utmp[i*mindim+j];
            }
         }

         for( size_t i=0UL; i<(*V).rows(); ++i ) {
            for( size_t j=0UL; j<(*V).columns(); ++j ) {
               (*V)(i,j) = vtmp[i*N+j];
            }
         }
      }
   }

   return num;
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
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes all min(\a m,\a n) singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The resulting real and non-negative
// singular values are stored in descending order in the given vector \a s, the resulting left
// singular vectors are stored in the given matrix \a U, and the resulting right singular vectors
// are stred in the given matrix \a V. \a s, \a U, and \a V are resized to the correct dimensions
// (if possible and necessary).
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
   DynamicMatrix<double,rowMajor>     V( 5UL, 8UL );  // The matrix for the right singular vectors

   gesvdx( A, U, s, V );
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
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
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                      DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V )
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

   using ET = ElementType_t<MT1>;
   using UT = UnderlyingElement_t<ET>;

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );

   resize( *s, mindim, false );
   resize( *U, M, mindim, false );
   resize( *V, mindim, N, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   return gesvdx_backend( *A, *U, *s, *V, 'A', UT(), UT(), 0, 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the singular value decomposition (SVD) of the given dense general matrix.
// \ingroup lapack_singular_value
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Vector cannot be resized.
// \exception std::invalid_argument Matrix cannot be resized.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values of the given general \a m-by-\a n
// matrix \a A by means of the LAPACK gesvdx() functions. The number of singular values to be
// computed is specified by the lower bound \a low and the upper bound \a upp, which either form
// an integral or a floating point range.
//
// In case \a low and \a upp form are of integral type, the function computes all singular values
// in the index range \f$[low..upp]\f$. The \a num resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a \a num-dimensional vector. The resulting left singular vectors are stored
// in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-\a num matrix. The resulting right singular vectors are stored in the given matrix \a V,
// which is either resized (if possible) or expected to be a \a num-by-\a n matrix.
//
// In case \a low and \a upp are of floating point type, the function computes all singular values
// in the half-open interval \f$(low..upp]\f$. The resulting real and non-negative singular values
// are stored in descending order in the given vector \a s, which is either resized (if possible)
// or expected to be a min(\a m,\a n)-dimensional vector. The resulting left singular vectors are
// stored in the given matrix \a U, which is either resized (if possible) or expected to be a
// \a m-by-min(\a m,\a n) matrix. The resulting right singular vectors are stored in the given
// matrix \a V, which is either resized (if possible) or expected to be a min(\a m,\a n)-by-\a n
// matrix.
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
//  - ... the given scalar values don't form a proper range;
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
   DynamicMatrix<double,rowMajor>     V( 5UL, 8UL );  // The matrix for the right singular vectors

   gesvdx( A, U, s, V, 1.0, 2.0 );  // Computes all eigenvalues in the interval (1..2]
   \endcode

   \code
   using blaze::DynamicMatrix;
   using blaze::DynamicVector;
   using blaze::rowMajor;
   using blaze::columnVector;

   DynamicMatrix<double,rowMajor> A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U( 5UL, 3UL );  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s( 3UL );       // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V( 3UL, 8UL );  // The matrix for the right singular vectors

   gesvdx( A, U, s, V, 0, 2 );  // Computes the first three eigenvalues
   \endcode

// For more information on the gesvdx() functions (i.e. sgesvdx(), dgesvdx(), cgesvdx(), and
// zgesvdx()) see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename MT2   // Type of the matrix U
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT3   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline size_t gesvdx( DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                      DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, ST low, ST upp )
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

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   if( IsFloatingPoint_v<ST> && low >= upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid value range provided" );
   }

   if( !IsFloatingPoint_v<ST> && low > upp ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   const size_t M( (*A).rows() );
   const size_t N( (*A).columns() );
   const size_t mindim( min( M, N ) );
   const size_t expected( IsFloatingPoint_v<ST> ? mindim : size_t( upp - low ) + 1UL );

   if( !IsFloatingPoint_v<ST> && expected > mindim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range provided" );
   }

   resize( *s, expected, false );
   resize( *U, M, expected, false );
   resize( *V, expected, N, false );

   if( M == 0UL || N == 0UL ) {
      return 0;
   }

   const char       range( IsFloatingPoint_v<ST> ? 'V' : 'I' );
   const ST         vl   ( IsFloatingPoint_v<ST> ? low : ST() );
   const ST         vu   ( IsFloatingPoint_v<ST> ? upp : ST() );
   const blas_int_t il   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( low ) );
   const blas_int_t iu   ( IsFloatingPoint_v<ST> ? 0 : numeric_cast<blas_int_t>( upp ) );

   const size_t actual( gesvdx_backend( *A, *U, *s, *V, range, vl, vu, il, iu ) );

   if( IsResizable_v<VT> ) {
      resize( *s, actual, true );
   }

   if( IsResizable_v<MT2> ) {
      resize( *U, M, actual, true );
   }

   if( IsResizable_v<MT3> ) {
      resize( *V, actual, N, true );
   }

   return actual;
}
//*************************************************************************************************

} // namespace blaze

#endif
