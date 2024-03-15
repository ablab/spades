//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/getrf.h
//  \brief Header file for the CLAPACK getrf wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GETRF_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GETRF_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION) && !defined(BLAS_H)
extern "C" {

void sgetrf_( blaze::blas_int_t* m, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, blaze::blas_int_t* info );
void dgetrf_( blaze::blas_int_t* m, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, blaze::blas_int_t* info );
void cgetrf_( blaze::blas_int_t* m, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, blaze::blas_int_t* info );
void zgetrf_( blaze::blas_int_t* m, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, blaze::blas_int_t* info );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK LU DECOMPOSITION FUNCTIONS (GETRF)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LU decomposition functions (getrf) */
//@{
void getrf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda,
            blas_int_t* ipiv, blas_int_t* info );

void getrf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda,
            blas_int_t* ipiv, blas_int_t* info );

void getrf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda,
            blas_int_t* ipiv, blas_int_t* info );

void getrf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda,
            blas_int_t* ipiv, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general single precision
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \a m-by-\a n single
// precision column-major matrix based on the LAPACK sgetrf() function, which uses partial
// pivoting with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the sgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void getrf( blas_int_t m, blas_int_t n, float* A, blas_int_t lda,
                   blas_int_t* ipiv, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgetrf_( &m, &n, A, &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general double precision
//        column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \a m-by-\a n double
// precision column-major matrix based on the LAPACK dgetrf() function, which uses partial
// pivoting with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the dgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void getrf( blas_int_t m, blas_int_t n, double* A, blas_int_t lda,
                   blas_int_t* ipiv, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgetrf_( &m, &n, A, &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general single precision
//        complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \a m-by-\a n single
// precision complex column-major matrix based on the LAPACK cgetrf() function, which uses
// partial pivoting with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the cgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void getrf( blas_int_t m, blas_int_t n, complex<float>* A, blas_int_t lda,
                   blas_int_t* ipiv, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cgetrf_( &m, &n, reinterpret_cast<ET*>( A ), &lda, ipiv, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the LU decomposition of the given dense general double precision
//        complex column-major matrix.
// \ingroup lapack_decomposition
//
// \param m The number of rows of the given matrix \f$[0..\infty)\f$.
// \param n The number of columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array for the pivot indices; size >= min( \a m, \a n ).
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix LU decomposition of a general \a m-by-\a n double
// precision complex column-major matrix based on the LAPACK zgetrf() function, which uses
// partial pivoting with row interchanges. The resulting decomposition has the form

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix (lower trapezoidal if
// \a m > \a n), and \c U is an upper triangular matrix (upper trapezoidal if \a m < \a n). The
// resulting decomposition is stored within the matrix \a A: \c L is stored in the lower part of
// \a A and \c U is stored in the upper part. The unit diagonal elements of \c L are not stored.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The decomposition finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but the factor U(i,i) is singular.
//
// For more information on the zgetrf() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void getrf( blas_int_t m, blas_int_t n, complex<double>* A, blas_int_t lda,
                   blas_int_t* ipiv, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zgetrf_( &m, &n, reinterpret_cast<ET*>( A ), &lda, ipiv, info );
}
//*************************************************************************************************

} // namespace blaze

#endif
