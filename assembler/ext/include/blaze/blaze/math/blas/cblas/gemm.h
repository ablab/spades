//=================================================================================================
/*!
//  \file blaze/math/blas/cblas/gemm.h
//  \brief Header file for the CBLAS gemm wrapper functions
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

#ifndef _BLAZE_MATH_BLAS_CBLAS_GEMM_H_
#define _BLAZE_MATH_BLAS_CBLAS_GEMM_H_


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
//  BLAS GENERAL MATRIX MULTIPLICATION FUNCTIONS (GEMM)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS general matrix multiplication functions (gemm) */
//@{
#if BLAZE_BLAS_MODE

void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
           blas_int_t m, blas_int_t n, blas_int_t k, float alpha, const float* A,
           blas_int_t lda, const float* B, blas_int_t ldb, float beta, float* C,
           blas_int_t ldc );

void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
           blas_int_t m, blas_int_t n, blas_int_t k, double alpha, const double* A,
           blas_int_t lda, const double* B, blas_int_t ldb, double beta, float* C,
           blas_int_t ldc );

void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
           blas_int_t m, blas_int_t n, blas_int_t k, complex<float> alpha,
           const complex<float>* A, blas_int_t lda, const complex<float>* B,
           blas_int_t ldb, complex<float> beta, complex<float>* C, blas_int_t ldc );

void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
           blas_int_t m, blas_int_t n, blas_int_t k, complex<double> alpha,
           const complex<double>* A, blas_int_t lda, const complex<double>* B,
           blas_int_t ldb, complex<double> beta, complex<double>* C, blas_int_t ldc );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with single precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param transB Specifies whether to transpose matrix \a B (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A and \a C \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B and \a C \f$[0..\infty)\f$.
// \param k The number of columns of matrix \a A and rows in matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param beta The scaling factor for \f$ C \f$.
// \param C Pointer to the first element of matrix \a C.
// \param ldc The total number of elements between two rows/columns of matrix \a C \f$[0..\infty)\f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for general single
// precision matrices based on the cblas_sgemm() function (\f$ C=\alpha*A*B+\beta*C \f$).
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                  blas_int_t m, blas_int_t n, blas_int_t k, float alpha, const float* A,
                  blas_int_t lda, const float* B, blas_int_t ldb, float beta, float* C,
                  blas_int_t ldc )
{
   cblas_sgemm( order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with double precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param transB Specifies whether to transpose matrix \a B (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A and \a C \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B and \a C \f$[0..\infty)\f$.
// \param k The number of columns of matrix \a A and rows in matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param beta The scaling factor for \f$ C \f$.
// \param C Pointer to the first element of matrix \a C.
// \param ldc The total number of elements between two rows/columns of matrix \a C \f$[0..\infty)\f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for double precision
// matrices based on the cblas_dgemm() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                  blas_int_t m, blas_int_t n, blas_int_t k, double alpha, const double* A,
                  blas_int_t lda, const double* B, blas_int_t ldb, double beta, double* C,
                  blas_int_t ldc )
{
   cblas_dgemm( order, transA, transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with single precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param transB Specifies whether to transpose matrix \a B (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A and \a C \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B and \a C \f$[0..\infty)\f$.
// \param k The number of columns of matrix \a A and rows in matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param beta The scaling factor for \f$ C \f$.
// \param C Pointer to the first element of matrix \a C.
// \param ldc The total number of elements between two rows/columns of matrix \a C \f$[0..\infty)\f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for single precision
// complex matrices based on the cblas_cgemm() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                  blas_int_t m, blas_int_t n, blas_int_t k, complex<float> alpha,
                  const complex<float>* A, blas_int_t lda, const complex<float>* B,
                  blas_int_t ldb, complex<float> beta, complex<float>* C, blas_int_t ldc )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

   cblas_cgemm( order, transA, transB, m, n, k, reinterpret_cast<const float*>( &alpha ),
                reinterpret_cast<const float*>( A ), lda, reinterpret_cast<const float*>( B ),
                ldb, reinterpret_cast<const float*>( &beta ), reinterpret_cast<float*>( C ), ldc );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication with double precision
//        matrices (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup blas
//
// \param order Specifies the storage order of matrix \a A (\a CblasColMajor or \a CblasColMajor).
// \param transA Specifies whether to transpose matrix \a A (\a CblasNoTrans or \a CblasTrans).
// \param transB Specifies whether to transpose matrix \a B (\a CblasNoTrans or \a CblasTrans).
// \param m The number of rows of matrix \a A and \a C \f$[0..\infty)\f$.
// \param n The number of columns of matrix \a B and \a C \f$[0..\infty)\f$.
// \param k The number of columns of matrix \a A and rows in matrix \a B \f$[0..\infty)\f$.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param A Pointer to the first element of matrix \a A.
// \param lda The total number of elements between two rows/columns of matrix \a A \f$[0..\infty)\f$.
// \param B Pointer to the first element of matrix \a B.
// \param ldb The total number of elements between two rows/columns of matrix \a B \f$[0..\infty)\f$.
// \param beta The scaling factor for \f$ C \f$.
// \param C Pointer to the first element of matrix \a C.
// \param ldc The total number of elements between two rows/columns of matrix \a C \f$[0..\infty)\f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication for double precision
// complex matrices based on the cblas_zgemm() function.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gemm( CBLAS_ORDER order, CBLAS_TRANSPOSE transA, CBLAS_TRANSPOSE transB,
                  blas_int_t m, blas_int_t n, blas_int_t k, complex<double> alpha,
                  const complex<double>* A, blas_int_t lda, const complex<double>* B,
                  blas_int_t ldb, complex<double> beta, complex<double>* C, blas_int_t ldc )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

   cblas_zgemm( order, transA, transB, m, n, k, reinterpret_cast<const double*>( &alpha ),
                reinterpret_cast<const double*>( A ), lda, reinterpret_cast<const double*>( B ),
                ldb, reinterpret_cast<const double*>( &beta ), reinterpret_cast<double*>( C ), ldc );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
