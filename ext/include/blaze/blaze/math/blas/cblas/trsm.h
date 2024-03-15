//=================================================================================================
/*!
//  \file blaze/math/blas/cblas/trsm.h
//  \brief Header file for the CBLAS trsm wrapper functions
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

#ifndef _BLAZE_MATH_BLAS_CBLAS_TRSM_H_
#define _BLAZE_MATH_BLAS_CBLAS_TRSM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/Types.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/Complex.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  BLAS TRIANGULAR LSE SOLVER FUNCTIONS (TRSM)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS triangular LSE solver functions (trsm) */
//@{
#if BLAZE_BLAS_MODE

void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
           CBLAS_DIAG diag, blas_int_t m, blas_int_t n, float alpha, const float* A,
           blas_int_t lda, float* B, blas_int_t ldb );

void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
           CBLAS_DIAG diag, blas_int_t m, blas_int_t n, double alpha, const double* A,
           blas_int_t lda, double* B, blas_int_t ldb );

void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
           CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<float> alpha,
           const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb );

void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
           CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<double> alpha,
           const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for solving a triangular system of equations with single precision matrices
//        (\f$ A*X=\alpha*B \f$ or \f$ X*A=\alpha*B \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a B (\a CblasRowMajor or \a CblasColMajor).
// \param side \a CblasLeft to compute \f$ A*X=\alpha*B \f$, \a CblasRight to compute \f$ X*A=\alpha*B \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param m The number of rows of matrix \a B \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ B \f$.
// \param A Pointer to the first element of the triangular matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \return void
//
// This function solves a triangular system of equations with multiple values for the right side
// based on the cblas_strsm() function. During the solution process, matrix \a B is overwritten
// with the resulting matrix \a X. Note that matrix \a A is expected to be a square matrix.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                  CBLAS_DIAG diag, blas_int_t m, blas_int_t n, float alpha, const float* A,
                  blas_int_t lda, float* B, blas_int_t ldb )
{
   cblas_strsm( order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for solving a triangular system of equations with double precision matrices
//        (\f$ A*X=\alpha*B \f$ or \f$ X*A=\alpha*B \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a B (\a CblasRowMajor or \a CblasColMajor).
// \param side \a CblasLeft to compute \f$ A*X=\alpha*B \f$, \a CblasRight to compute \f$ X*A=\alpha*B \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param m The number of rows of matrix \a B \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ B \f$.
// \param A Pointer to the first element of the triangular matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \return void
//
// This function solves a triangular system of equations with multiple values for the right side
// based on the cblas_dtrsm() function. During the solution process, matrix \a B is overwritten
// with the resulting matrix \a X. Note that matrix \a A is expected to be a square matrix.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                  CBLAS_DIAG diag, blas_int_t m, blas_int_t n, double alpha, const double* A,
                  blas_int_t lda, double* B, blas_int_t ldb )
{
   cblas_dtrsm( order, side, uplo, transA, diag, m, n, alpha, A, lda, B, ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for solving a triangular system of equations with single precision complex
//        matrices (\f$ A*X=\alpha*B \f$ or \f$ X*A=\alpha*B \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a B (\a CblasRowMajor or \a CblasColMajor).
// \param side \a CblasLeft to compute \f$ A*X=\alpha*B \f$, \a CblasRight to compute \f$ X*A=\alpha*B \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param m The number of rows of matrix \a B \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ B \f$.
// \param A Pointer to the first element of the triangular matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \return void
//
// This function solves a triangular system of equations with multiple values for the right side
// based on the cblas_ctrsm() function. During the solution process, matrix \a B is overwritten
// with the resulting matrix \a X. Note that matrix \a A is expected to be a square matrix.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                  CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<float> alpha,
                  const complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cblas_ctrsm( order, side, uplo, transA, diag, m, n, reinterpret_cast<const float*>( &alpha ),
                reinterpret_cast<const float*>( A ), lda, reinterpret_cast<float*>( B ), ldb );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for solving a triangular system of equations with double precision complex
//        matrices (\f$ A*X=\alpha*B \f$ or \f$ X*A=\alpha*B \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a B (\a CblasRowMajor or \a CblasColMajor).
// \param side \a CblasLeft to compute \f$ A*X=\alpha*B \f$, \a CblasRight to compute \f$ X*A=\alpha*B \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param diag Specifies whether \a A is unitriangular (\a CblasNonUnit or \a CblasUnit).
// \param m The number of rows of matrix \a B \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ B \f$.
// \param A Pointer to the first element of the triangular matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \return void
//
// This function solves a triangular system of equations with multiple values for the right side
// based on the cblas_ztrsm() function. During the solution process, matrix \a B is overwritten
// with the resulting matrix \a X. Note that matrix \a A is expected to be a square matrix.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void trsm( CBLAS_ORDER order, CBLAS_SIDE side, CBLAS_UPLO uplo, CBLAS_TRANSPOSE transA,
                  CBLAS_DIAG diag, blas_int_t m, blas_int_t n, complex<double> alpha,
                  const complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   cblas_ztrsm( order, side, uplo, transA, diag, m, n, reinterpret_cast<const double*>( &alpha ),
                reinterpret_cast<const double*>( A ), lda, reinterpret_cast<double*>( B ), ldb );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
