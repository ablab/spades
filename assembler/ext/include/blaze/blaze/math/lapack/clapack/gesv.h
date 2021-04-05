//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/gesv.h
//  \brief Header file for the CLAPACK gesv wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GESV_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GESV_H_


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

void sgesv_( blaze::blas_int_t* n, blaze::blas_int_t* nrhs, float* A, blaze::blas_int_t* lda,
             blaze::blas_int_t* ipiv, float* b, blaze::blas_int_t* ldb, blaze::blas_int_t* info );
void dgesv_( blaze::blas_int_t* n, blaze::blas_int_t* nrhs, double* A, blaze::blas_int_t* lda,
             blaze::blas_int_t* ipiv, double* b, blaze::blas_int_t* ldb, blaze::blas_int_t* info );
void cgesv_( blaze::blas_int_t* n, blaze::blas_int_t* nrhs, float* A, blaze::blas_int_t* lda,
             blaze::blas_int_t* ipiv, float* b, blaze::blas_int_t* ldb, blaze::blas_int_t* info );
void zgesv_( blaze::blas_int_t* n, blaze::blas_int_t* nrhs, double* A, blaze::blas_int_t* lda,
             blaze::blas_int_t* ipiv, double* b, blaze::blas_int_t* ldb, blaze::blas_int_t* info );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK GENERAL LINEAR SYSTEM FUNCTIONS (GESV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK general linear system functions (gesv) */
//@{
void gesv( blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda,
           blas_int_t* ipiv, float* B, blas_int_t ldb, blas_int_t* info );

void gesv( blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda,
           blas_int_t* ipiv, double* B, blas_int_t ldb, blas_int_t* info );

void gesv( blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda,
           blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info );

void gesv( blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda,
           blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a general single precision linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK sgesv() function to compute the solution to the general system
// of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n matrix and \a X and \a B are
// \a n-by-\a nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular the solution could not be computed.
//
// For more information on the sgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gesv( blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda,
                  blas_int_t* ipiv, float* B, blas_int_t ldb, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgesv_( &n, &nrhs, A, &lda, ipiv, B, &ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a general double precision linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK dgesv() function to compute the solution to the general system
// of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n matrix and \a X and \a B are
// \a n-by-\a nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular the solution could not be computed.
//
// For more information on the dgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gesv( blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda,
                  blas_int_t* ipiv, double* B, blas_int_t ldb, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgesv_( &n, &nrhs, A, &lda, ipiv, B, &ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a general single precision complex linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK cgesv() function to compute the solution to the general system
// of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n matrix and \a X and \a B are
// \a n-by-\a nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular the solution could not be computed.
//
// For more information on the cgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gesv( blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda,
                  blas_int_t* ipiv, complex<float>* B, blas_int_t ldb, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cgesv_( &n, &nrhs, reinterpret_cast<ET*>( A ), &lda, ipiv,
           reinterpret_cast<ET*>( B ), &ldb, info );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a general double precision complex linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK zgesv() function to compute the solution to the general system
// of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n matrix and \a X and \a B are
// \a n-by-\a nrhs matrices.
//
// The LU decomposition with partial pivoting and row interchanges is used to factor \a A as

                          \f[ A = P \cdot L \cdot U, \f]

// where \c P is a permutation matrix, \c L is a lower unitriangular matrix, and \c U is an upper
// triangular matrix. The resulting decomposition is stored within \a A: \c L is stored in the
// lower part of \a A and \c U is stored in the upper part. The unit diagonal elements of \c L
// are not stored. The factored form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the decomposition has been completed, but since factor U(i,i) is exactly
//          singular the solution could not be computed.
//
// For more information on the zgesv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gesv( blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda,
                  blas_int_t* ipiv, complex<double>* B, blas_int_t ldb, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zgesv_( &n, &nrhs, reinterpret_cast<ET*>( A ), &lda, ipiv,
           reinterpret_cast<ET*>( B ), &ldb, info );
}
//*************************************************************************************************

} // namespace blaze

#endif
