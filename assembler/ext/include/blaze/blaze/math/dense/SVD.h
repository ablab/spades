//=================================================================================================
/*!
//  \file blaze/math/dense/SVD.h
//  \brief Header file for the dense matrix singular value decomposition (SVD)
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

#ifndef _BLAZE_MATH_DENSE_SVD_H_
#define _BLAZE_MATH_DENSE_SVD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/lapack/gesdd.h>
#include <blaze/math/lapack/gesvdx.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/util/mpl/If.h>


namespace blaze {

//=================================================================================================
//
//  SINGULAR VALUE DECOMPOSITION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Singular value decomposition functions */
//@{
template< typename MT, bool SO, typename VT, bool TF >
inline void svd( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s );

template< typename MT1, bool SO, typename VT, bool TF, typename MT2, typename MT3 >
inline void svd( const DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                 DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V );

template< typename MT, bool SO, typename VT, bool TF, typename ST >
inline size_t svd( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, ST low, ST upp );

template< typename MT1, bool SO, typename VT, bool TF, typename MT2, typename MT3, typename ST >
inline size_t svd( const DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                   DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, ST low, ST upp );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Singular value decomposition (SVD) of the given dense general matrix.
// \ingroup dense_matrix
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \return void
// \exception std::invalid_argument Size of fixed size vector does not match.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix.
// The resulting min(\a m,\a n) singular values are stored in the given vector \a s, which is
// resized to the correct size (if possible and necessary).
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

   DynamicMatrix<double,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s;  // The vector for the singular values

   svd( A, s );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
//
// \note Further options for computing singular values and singular vectors are available via the
// gesvd(), gesdd(), and gesvdx() functions.
*/
template< typename MT  // Type of the matrix A
        , bool SO      // Storage order of the matrix A
        , typename VT  // Type of the vector s
        , bool TF >    // Transpose flag of the vector s
inline void svd( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   using ATmp = ResultType_t< RemoveAdaptor_t<MT> >;
   using STmp = If_t< IsContiguous_v<VT>, VT&, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( ATmp );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<ATmp> );

   ATmp Atmp( *A );
   STmp stmp( *s );

   gesdd( Atmp, stmp );

   if( !IsContiguous_v<VT> ) {
      (*s) = stmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Singular value decomposition (SVD) of the given dense general matrix.
// \ingroup dense_matrix
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \return void
// \exception std::invalid_argument Dimensions of fixed size matrix U do not match.
// \exception std::invalid_argument Size of fixed size vector does not match.
// \exception std::invalid_argument Dimensions of fixed size matrix V do not match.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function performs the singular value decomposition of a general \a m-by-\a n matrix. The
// resulting min(\a m,\a n) singular values are stored in the given vector \a s, the resulting
// left singular vectors are stored in the given matrix \a U, and the resulting right singular
// vectors are stored in the given matrix \a V. \a s, \a U and \a V are resized to the correct
// dimensions (if possible and necessary).
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

   DynamicMatrix<double,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U;  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s;  // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V;  // The matrix for the right singular vectors

   svd( A, U, s, V );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
//
// \note Further options for computing singular values and singular vectors are available via the
// gesvd(), gesdd(), and gesvdx() functions.
*/
template< typename MT1    // Type of the matrix A
        , bool SO         // Storage order of all matrices
        , typename VT     // Type of the vector s
        , bool TF         // Transpose flag of the vector s
        , typename MT2    // Type of the matrix U
        , typename MT3 >  // Type of the matrix V
inline void svd( const DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                 DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   using ATmp = ResultType_t< RemoveAdaptor_t<MT1> >;
   using UTmp = If_t< IsContiguous_v<MT2>, MT2&, ResultType_t<MT2> >;
   using STmp = If_t< IsContiguous_v<VT>, VT&, ResultType_t<VT> >;
   using VTmp = If_t< IsContiguous_v<MT3>, MT3&, ResultType_t<MT3> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( ATmp );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<ATmp> );

   ATmp Atmp( *A );
   UTmp Utmp( *U );
   STmp stmp( *s );
   VTmp Vtmp( *V );

   gesdd( Atmp, Utmp, stmp, Vtmp, 'S' );

   if( !IsContiguous_v<MT2> ) {
      (*U) = Utmp;
   }

   if( !IsContiguous_v<VT> ) {
      (*s) = stmp;
   }

   if( !IsContiguous_v<MT3> ) {
      (*V) = Vtmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Singular value decomposition (SVD) of the given dense general matrix.
// \ingroup dense_matrix
//
// \param A The given general matrix.
// \param s The resulting vector of singular values.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Size of fixed size vector does not match.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values of the given general \a m-by-\a n
// matrix. The number of singular values to be computed is specified by the lower bound \a low and
// the upper bound \a upp, which either form an integral or a floating point range.
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

   DynamicMatrix<double,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicVector<double,columnVector> s;  // The vector for the singular values

   svd( A, s, 0, 2 );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
//
// \note Further options for computing singular values and singular vectors are available via the
// gesvd(), gesdd(), and gesvdx() functions.
*/
template< typename MT    // Type of the matrix A
        , bool SO        // Storage order of the matrix A
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename ST >  // Type of the scalar boundary values
inline size_t svd( const DenseMatrix<MT,SO>& A, DenseVector<VT,TF>& s, ST low, ST upp )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   using ATmp = ResultType_t< RemoveAdaptor_t<MT> >;
   using STmp = If_t< IsContiguous_v<VT>, VT&, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( ATmp );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<ATmp> );

   ATmp Atmp( *A );
   STmp stmp( *s );

   const auto num = gesvdx( Atmp, stmp, low, upp );

   if( !IsContiguous_v<VT> ) {
      (*s) = stmp;
   }

   return num;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Singular value decomposition (SVD) of the given dense general matrix.
// \ingroup dense_matrix
//
// \param A The given general matrix.
// \param U The resulting matrix of left singular vectors.
// \param s The resulting vector of singular values.
// \param V The resulting matrix of right singular vectors.
// \param low The lower bound of the interval to be searched for singular values.
// \param upp The upper bound of the interval to be searched for singular values.
// \return The total number of singular values found.
// \exception std::invalid_argument Dimensions of fixed size matrix U do not match.
// \exception std::invalid_argument Size of fixed size vector does not match.
// \exception std::invalid_argument Dimensions of fixed size matrix V do not match.
// \exception std::invalid_argument Invalid value range provided.
// \exception std::invalid_argument Invalid index range provided.
// \exception std::runtime_error Singular value decomposition failed.
//
// This function computes a specified number of singular values and singular vectors of the given
// general \a m-by-\a n matrix. The number of singular values and vectors to be computed is
// specified by the lower bound \a low and the upper bound \a upp, which either form an integral
// or a floating point range.
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
// matrix \a V, which is either resized (if possible) or expected to be a \a min(\a m,\a n)-by-\a n
// matrix.
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

   DynamicMatrix<double,rowMajor>  A( 5UL, 8UL );  // The general matrix A
   // ... Initialization

   DynamicMatrix<double,rowMajor>     U;  // The matrix for the left singular vectors
   DynamicVector<double,columnVector> s;  // The vector for the singular values
   DynamicMatrix<double,rowMajor>     V;  // The matrix for the right singular vectors

   svd( A, U, s, V, 0, 2 );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
//
// \note Further options for computing singular values and singular vectors are available via the
// gesvd(), gesdd(), and gesvdx() functions.
*/
template< typename MT1   // Type of the matrix A
        , bool SO        // Storage order of all matrices
        , typename VT    // Type of the vector s
        , bool TF        // Transpose flag of the vector s
        , typename MT2   // Type of the matrix U
        , typename MT3   // Type of the matrix V
        , typename ST >  // Type of the scalar boundary values
inline size_t svd( const DenseMatrix<MT1,SO>& A, DenseMatrix<MT2,SO>& U,
                   DenseVector<VT,TF>& s, DenseMatrix<MT3,SO>& V, ST low, ST upp )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT3 );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ST );

   using ATmp = ResultType_t< RemoveAdaptor_t<MT1> >;
   using UTmp = If_t< IsContiguous_v<MT2>, MT2&, ResultType_t<MT2> >;
   using STmp = If_t< IsContiguous_v<VT>, VT&, ResultType_t<VT> >;
   using VTmp = If_t< IsContiguous_v<MT3>, MT3&, ResultType_t<MT3> >;

   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( ATmp );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( ATmp );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<ATmp> );

   ATmp Atmp( *A );
   UTmp Utmp( *U );
   STmp stmp( *s );
   VTmp Vtmp( *V );

   const auto num = gesvdx( Atmp, Utmp, stmp, Vtmp, low, upp );

   if( !IsContiguous_v<MT2> ) {
      (*U) = Utmp;
   }

   if( !IsContiguous_v<VT> ) {
      (*s) = stmp;
   }

   if( !IsContiguous_v<MT3> ) {
      (*V) = Vtmp;
   }

   return num;
}
//*************************************************************************************************

} // namespace blaze

#endif
