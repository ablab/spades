//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/posv.h
//  \brief Header file for the CLAPACK posv wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_POSV_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_POSV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


//=================================================================================================
//
//  LAPACK FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if !defined(INTEL_MKL_VERSION)
extern "C" {

void sposv_( char* uplo, blaze::blas_int_t* n, blaze::blas_int_t* nrhs, float* A,
             blaze::blas_int_t* lda, float* b, blaze::blas_int_t* ldb,
             blaze::blas_int_t* info, blaze::fortran_charlen_t nuplo );
void dposv_( char* uplo, blaze::blas_int_t* n, blaze::blas_int_t* nrhs, double* A,
             blaze::blas_int_t* lda, double* b, blaze::blas_int_t* ldb,
             blaze::blas_int_t* info, blaze::fortran_charlen_t nuplo );
void cposv_( char* uplo, blaze::blas_int_t* n, blaze::blas_int_t* nrhs, float* A,
             blaze::blas_int_t* lda, float* b, blaze::blas_int_t* ldb,
             blaze::blas_int_t* info, blaze::fortran_charlen_t nuplo );
void zposv_( char* uplo, blaze::blas_int_t* n, blaze::blas_int_t* nrhs, double* A,
             blaze::blas_int_t* lda, double* b, blaze::blas_int_t* ldb,
             blaze::blas_int_t* info, blaze::fortran_charlen_t nuplo );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK POSITIVE DEFINITE LINEAR SYSTEM FUNCTIONS (POSV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK positive definite linear system functions (posv) */
//@{
void posv( char uplo, blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda,
           float* B, blas_int_t ldb, blas_int_t* info );

void posv( char uplo, blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda,
           double* B, blas_int_t ldb, blas_int_t* info );

void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda,
           complex<float>* B, blas_int_t ldb, blas_int_t* info );

void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda,
           complex<double>* B, blas_int_t ldb, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a positive definite single precision linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK sposv() function to compute the solution to the positive definite
// system of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n positive definite matrix
// and \a X and \a B are \a n-by-\a nrhs matrices.
//
// The Cholesky decomposition is used to factor \a A as

                      \f[ A = U^{T} U \texttt{ (if uplo = 'U'), or }
                          A = L L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U is an upper triangular matrix and \c L is a lower triangular matrix. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched. The factored
// form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the leading minor of order i is not positive definite, so the decomposition
//          could not be completed and the solution has not been computed.
//
// For more information on the sposv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void posv( char uplo, blas_int_t n, blas_int_t nrhs, float* A, blas_int_t lda,
                  float* B, blas_int_t ldb, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sposv_( &uplo, &n, &nrhs, A, &lda, B, &ldb, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a positive definite double precision linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK dposv() function to compute the solution to the positive definite
// system of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n positive definite matrix
// and \a X and \a B are \a n-by-\a nrhs matrices.
//
// The Cholesky decomposition is used to factor \a A as

                      \f[ A = U^{T} U \texttt{ (if uplo = 'U'), or }
                          A = L L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U is an upper triangular matrix and \c L is a lower triangular matrix. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched. The factored
// form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the leading minor of order i is not positive definite, so the decomposition
//          could not be completed and the solution has not been computed.
//
// For more information on the dposv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void posv( char uplo, blas_int_t n, blas_int_t nrhs, double* A, blas_int_t lda,
                  double* B, blas_int_t ldb, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dposv_( &uplo, &n, &nrhs, A, &lda, B, &ldb, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a positive definite single precision complex linear system of
//        equations (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK cposv() function to compute the solution to the positive definite
// system of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n positive definite matrix
// and \a X and \a B are \a n-by-\a nrhs matrices.
//
// The Cholesky decomposition is used to factor \a A as

                      \f[ A = U^{T} U \texttt{ (if uplo = 'U'), or }
                          A = L L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U is an upper triangular matrix and \c L is a lower triangular matrix. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched. The factored
// form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the leading minor of order i is not positive definite, so the decomposition
//          could not be completed and the solution has not been computed.
//
// For more information on the cposv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<float>* A, blas_int_t lda,
                  complex<float>* B, blas_int_t ldb, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cposv_( &uplo, &n, &nrhs, reinterpret_cast<ET*>( A ), &lda,
           reinterpret_cast<ET*>( B ), &ldb, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for solving a positive definite double precision complex linear system of
//        equations (\f$ A*X=B \f$).
// \ingroup lapack_solver
//
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param nrhs The number of right-hand side vectors \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major square matrix.
// \param lda The total number of elements between two columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of the column-major matrix.
// \param ldb The total number of elements between two columns of matrix \a B \f$[0..\infty)\f$.
// \param info Return code of the function call.
// \return void
//
// This function uses the LAPACK zposv() function to compute the solution to the positive definite
// system of linear equations \f$ A*X=B \f$, where \a A is a \a n-by-\a n positive definite matrix
// and \a X and \a B are \a n-by-\a nrhs matrices.
//
// The Cholesky decomposition is used to factor \a A as

                      \f[ A = U^{T} U \texttt{ (if uplo = 'U'), or }
                          A = L L^{T} \texttt{ (if uplo = 'L'), } \f]

// where \c U is an upper triangular matrix and \c L is a lower triangular matrix. The resulting
// decomposition is stored within \a A: In case \a uplo is set to \c 'L' the result is stored in
// the lower part of the matrix and the upper part remains untouched, in case \a uplo is set to
// \c 'U' the result is stored in the upper part and the lower part remains untouched. The factored
// form of \a A is then used to solve the system of equations.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The function finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the leading minor of order i is not positive definite, so the decomposition
//          could not be completed and the solution has not been computed.
//
// For more information on the zposv() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void posv( char uplo, blas_int_t n, blas_int_t nrhs, complex<double>* A, blas_int_t lda,
                  complex<double>* B, blas_int_t ldb, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zposv_( &uplo, &n, &nrhs, reinterpret_cast<ET*>( A ), &lda,
           reinterpret_cast<ET*>( B ), &ldb, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************

} // namespace blaze

#endif
