//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/gges.h
//  \brief Header file for the CLAPACK gges wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GGES_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GGES_H_


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

void sgges_( char* jobvsl, char* jobvsr, char* sort,
             blaze::blas_int_t (*selctg)( const float*, const float*, const float* ),
             blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda, float* B, blaze::blas_int_t* ldb,
             blaze::blas_int_t* sdim, float* alphar, float* alphai, float* beta, float* VSL,
             blaze::blas_int_t* ldvsl, float* VSR, blaze::blas_int_t* ldvsr, float* work,
             blaze::blas_int_t* lwork, blaze::blas_int_t* bwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr,
             blaze::fortran_charlen_t nsort );
void dgges_( char* jobvsl, char* jobvsr, char* sort,
             blaze::blas_int_t (*selctg)( const double*, const double*, const double* ),
             blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda, double* B, blaze::blas_int_t* ldb,
             blaze::blas_int_t* sdim, double* alphar, double* alphai, double* beta, double* VSL,
             blaze::blas_int_t* ldvsl, double* VSR, blaze::blas_int_t* ldvsr, double* work,
             blaze::blas_int_t* lwork, blaze::blas_int_t* bwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr,
             blaze::fortran_charlen_t nsort );
void cgges_( char* jobvsl, char* jobvsr, char* sort,
             blaze::blas_int_t (*selctg)( const float*, const float* ),
             blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda, float* B, blaze::blas_int_t* ldb,
             blaze::blas_int_t* sdim, float* alpha, float* beta, float* VSL, blaze::blas_int_t* ldvsl,
             float* VSR, blaze::blas_int_t* ldvsr, float* work, blaze::blas_int_t* lwork,
             float* rwork, blaze::blas_int_t* bwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr,
             blaze::fortran_charlen_t nsort );
void zgges_( char* jobvsl, char* jobvsr, char* sort,
             blaze::blas_int_t (*selctg)( const double*, const double* ),
             blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda, double* B, blaze::blas_int_t* ldb,
             blaze::blas_int_t* sdim, double* alpha, double* beta, double* VSL, blaze::blas_int_t* ldvsl,
             double* VSR, blaze::blas_int_t* ldvsr, double* work, blaze::blas_int_t* lwork,
             double* rwork, blaze::blas_int_t* bwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvsl, blaze::fortran_charlen_t njobvsr,
             blaze::fortran_charlen_t nsort );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK GENERALIZED SCHUR DECOMPOSITION FUNCTIONS (GGES)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK generalized Schur decomposition functions (gges) */
//@{
void gges( char jobvsl, char jobvsr, char sort,
           blas_int_t (*selctg)( const float*, const float*, const float* ), blas_int_t n,
           float* A, blas_int_t lda, float* B, blas_int_t ldb,
           blas_int_t* sdim, float* alphar, float* alphai, float* beta,
           float* VSL, blas_int_t ldvsl, float* VSR,
           blas_int_t ldvsr, float* work, blas_int_t lwork,
           blas_int_t* bwork, blas_int_t* info );

void gges( char jobvsl, char jobvsr, char sort,
           blas_int_t (*selctg)( const double*, const double*, const double* ), blas_int_t n,
           double* A, blas_int_t lda, double* B, blas_int_t ldb,
           blas_int_t* sdim, double* alphar, double* alphai, double* beta,
           double* VSL, blas_int_t ldvsl, double* VSR,
           blas_int_t ldvsr, double* work, blas_int_t lwork,
           blas_int_t* bwork, blas_int_t* info );

void gges( char jobvsl, char jobvsr, char sort,
           blas_int_t (*selctg)( const complex<float>*, const complex<float>* ), blas_int_t n,
           complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb,
           blas_int_t* sdim, complex<float>* alpha, complex<float>* beta,
           complex<float>* VSL, blas_int_t ldvsl, complex<float>* VSR,
           blas_int_t ldvsr, complex<float>* work, blas_int_t lwork,
           float* rwork, blas_int_t* bwork, blas_int_t* info );

void gges( char jobvsl, char jobvsr, char sort,
           blas_int_t (*selctg)( const complex<double>*, const complex<double>* ), blas_int_t n,
           complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb,
           blas_int_t* sdim, complex<double>* alpha, complex<double>* beta,
           complex<double>* VSL, blas_int_t ldvsl, complex<double>* VSR,
           blas_int_t ldvsr, complex<double>* work, blas_int_t lwork,
           double* rwork, blas_int_t* bwork, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition for a pair of
//        non-symmetric single precision column-major matrices.
// \ingroup lapack_eigenvalue
//
// \param jobvsl \c 'V' to compute the left Schur vectors, \c 'N' to not compute them.
// \param jobvsr \c 'V' to compute the right Schur vectors, \c 'N' to not compute them.
// \param sort \c 'S' to order the eigenvalues on the diagonal, \c 'N' to not order them (see \a selctg).
// \param selctg Logical function of three single precision arguments.
// \param n The order of the matrices \a A, \a B, \a VSL, and \a VSR \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix \a A.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param B Pointer to the first element of the single precision column-major matrix \a B.
// \param ldb The total number of elements between two columns of the matrix B; size >= max( 1, \a n ).
// \param sdim Returns the number of eigenvalues (after sorting) for which \a selctg is \a true.
// \param alphar Pointer to the first element of the real part of the eigenvalue numerator; size == n.
// \param alphai Pointer to the first element of the imaginary part of the eigenvalue numerator; size == n.
// \param beta Pointer to the first element of the eigenvalue denominator; size == n.
// \param VSL Pointer to the first element of the column-major matrix for the left Schur vectors.
// \param ldvsl The total number of elements between two columns of the matrix VSL \f$[0..\infty)\f$.
// \param VSR Pointer to the first element of the column-major matrix for the right Schur vectors.
// \param ldvsr The total number of elements between two columns of the matrix VSR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param bwork Auxiliary array; size == n.
// \param info Return code of the function call.
//
// This function computes for a pair of \a n-by-\a n single precision non-symmetric matrices
// (\a A,\a B) the generalized eigenvalues, the generalized real Schur form (\a S, \a T), and
// optionally the left and/or right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

                      \f[ (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T ). \f]

// Optionally, it also orders the eigenvalues such that a selected cluster of eigenvalues appears
// in the leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The leading columns of \a VSL and \a VSR then form an orthonormal basis for the
// corresponding left and right eigenspaces (deflating subspaces).
//
// The parameter \a jobvsl specifies the computation of the left Schur vectors:
//
//   - \c 'V': The left Schur vectors of \a (A,B) are computed and returned in \a VSL.
//   - \c 'N': The left Schur vectors of \a (A,B) are not computed.
//
// The parameter \a jobvsr specifies the computation of the right Schur vectors:
//
//   - \c 'V': The right Schur vectors of \a (A,B) are computed and returned in \a VSR.
//   - \c 'N': The right Schur vectors of \a (A,B) are not computed.
//
// Via the parameter \a selctg it is possible to select specific eigenvalues. \a selctg is
// a pointer to a function of three single precision floating point arguments returning a
// \c blas_int_t value. In case \a sort is set to \c 'N', \a selctg is not referenced, else
// if sort is set to \c 'S' \a selctg is used to select eigenvalues to sort to the top left
// of the Schur form. An eigenvalue \f$ (alphar[j]+alphai[j]*i)/beta[j] \f$ is selected if

   \code
   selctg( &alphar[j], &alphai[j], &beta[j])
   \endcode

// evaluates to true; i.e. if either one of a complex conjugate pair of eigenvalues is selected,
// then both complex eigenvalues are selected. On exit, the parameter \a sdim will contain the
// number of eigenvalues (after sorting) for which \a selctg is \a true. Note that complex
// conjugate pairs for which \a selctg is \a true for either eigenvalue count as 2.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 1,...,\a n: The QZ iteration failed. (A,B) are not in Schur form, but \a alphar(j),
//          \a alphai(j), and \a beta(j) should be correct for j=\a info,...,\a n-1.
//   - > \a n: see online reference for details.
//
// For more information on the sgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort,
                  blas_int_t (*selctg)( const float*, const float*, const float* ), blas_int_t n,
                  float* A, blas_int_t lda, float* B, blas_int_t ldb,
                  blas_int_t* sdim, float* alphar, float* alphai, float* beta,
                  float* VSL, blas_int_t ldvsl, float* VSR,
                  blas_int_t ldvsr, float* work, blas_int_t lwork,
                  blas_int_t* bwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgges_( &jobvsl, &jobvsr, &sort, selctg, &n, A, &lda, B, &ldb, sdim, alphar, alphai, beta,
           VSL, &ldvsl, VSR, &ldvsr, work, &lwork, bwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition for a pair of
//        non-symmetric double precision column-major matrices.
// \ingroup lapack_eigenvalue
//
// \param jobvsl \c 'V' to compute the left Schur vectors, \c 'N' to not compute them.
// \param jobvsr \c 'V' to compute the right Schur vectors, \c 'N' to not compute them.
// \param sort \c 'S' to order the eigenvalues on the diagonal, \c 'N' to not order them (see \a selctg).
// \param selctg Logical function of three double precision arguments.
// \param n The order of the matrices \a A, \a B, \a VSL, and \a VSR \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix \a A.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param B Pointer to the first element of the double precision column-major matrix \a B.
// \param ldb The total number of elements between two columns of the matrix B; size >= max( 1, \a n ).
// \param sdim Returns the number of eigenvalues (after sorting) for which \a selctg is \a true.
// \param alphar Pointer to the first element of the real part of the eigenvalue numerator; size == n.
// \param alphai Pointer to the first element of the imaginary part of the eigenvalue numerator; size == n.
// \param beta Pointer to the first element of the eigenvalue denominator; size == n.
// \param VSL Pointer to the first element of the column-major matrix for the left Schur vectors.
// \param ldvsl The total number of elements between two columns of the matrix VSL \f$[0..\infty)\f$.
// \param VSR Pointer to the first element of the column-major matrix for the right Schur vectors.
// \param ldvsr The total number of elements between two columns of the matrix VSR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param bwork Auxiliary array; size == n.
// \param info Return code of the function call.
//
// This function computes for a pair of \a n-by-\a n double precision non-symmetric matrices
// (\a A,\a B) the generalized eigenvalues, the generalized real Schur form (\a S, \a T), and
// optionally the left and/or right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

                      \f[ (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T ). \f]

// Optionally, it also orders the eigenvalues such that a selected cluster of eigenvalues appears
// in the leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The leading columns of \a VSL and \a VSR then form an orthonormal basis for the
// corresponding left and right eigenspaces (deflating subspaces).
//
// The parameter \a jobvsl specifies the computation of the left Schur vectors:
//
//   - \c 'V': The left Schur vectors of \a (A,B) are computed and returned in \a VSL.
//   - \c 'N': The left Schur vectors of \a (A,B) are not computed.
//
// The parameter \a jobvsr specifies the computation of the right Schur vectors:
//
//   - \c 'V': The right Schur vectors of \a (A,B) are computed and returned in \a VSR.
//   - \c 'N': The right Schur vectors of \a (A,B) are not computed.
//
// Via the parameter \a selctg it is possible to select specific eigenvalues. \a selctg is
// a pointer to a function of three double precision floating point arguments returning an
// \c blas_int_t value. In case \a sort is set to \c 'N', \a selctg is not referenced, else
// if sort is set to \c 'S' \a selctg is used to select eigenvalues to sort to the top left
// of the Schur form. An eigenvalue \f$ (alphar[j]+alphai[j]*i)/beta[j] \f$ is selected if

   \code
   selctg( &alphar[j], &alphai[j], &beta[j])
   \endcode

// evaluates to true; i.e. if either one of a complex conjugate pair of eigenvalues is selected,
// then both complex eigenvalues are selected. On exit, the parameter \a sdim will contain the
// number of eigenvalues (after sorting) for which \a selctg is \a true. Note that complex
// conjugate pairs for which \a selctg is \a true for either eigenvalue count as 2.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 1,...,\a n: The QZ iteration failed. (A,B) are not in Schur form, but \a alphar(j),
//          \a alphai(j), and \a beta(j) should be correct for j=\a info,...,\a n-1.
//   - > \a n: see online reference for details.
//
// For more information on the sgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort,
                  blas_int_t (*selctg)( const double*, const double*, const double* ), blas_int_t n,
                  double* A, blas_int_t lda, double* B, blas_int_t ldb,
                  blas_int_t* sdim, double* alphar, double* alphai, double* beta,
                  double* VSL, blas_int_t ldvsl, double* VSR,
                  blas_int_t ldvsr, double* work, blas_int_t lwork,
                  blas_int_t* bwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgges_( &jobvsl, &jobvsr, &sort, selctg, &n, A, &lda, B, &ldb, sdim, alphar, alphai, beta,
           VSL, &ldvsl, VSR, &ldvsr, work, &lwork, bwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition for a pair of
//        non-symmetric single precision complex column-major matrices.
// \ingroup lapack_eigenvalue
//
// \param jobvsl \c 'V' to compute the left Schur vectors, \c 'N' to not compute them.
// \param jobvsr \c 'V' to compute the right Schur vectors, \c 'N' to not compute them.
// \param sort \c 'S' to order the eigenvalues on the diagonal, \c 'N' to not order them (see \a selctg).
// \param selctg Logical function of two single precision complex arguments.
// \param n The order of the matrices \a A, \a B, \a VSL, and \a VSR \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix \a A.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param B Pointer to the first element of the single precision complex column-major matrix \a B.
// \param ldb The total number of elements between two columns of the matrix B; size >= max( 1, \a n ).
// \param sdim Returns the number of eigenvalues (after sorting) for which \a selctg is \a true.
// \param alpha Pointer to the first element of the eigenvalue numerator; size == n.
// \param beta Pointer to the first element of the eigenvalue denominator; size == n.
// \param VSL Pointer to the first element of the column-major matrix for the left Schur vectors.
// \param ldvsl The total number of elements between two columns of the matrix VSL \f$[0..\infty)\f$.
// \param VSR Pointer to the first element of the column-major matrix for the right Schur vectors.
// \param ldvsr The total number of elements between two columns of the matrix VSR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size == 8*n.
// \param bwork Auxiliary array; size == n.
// \param info Return code of the function call.
//
// This function computes for a pair of \a n-by-\a n single precision complex non-symmetric
// matrices (\a A,\a B) the generalized eigenvalues, the generalized real Schur form (\a S, \a T),
// and optionally the left and/or right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

                      \f[ (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T ). \f]

// Optionally, it also orders the eigenvalues such that a selected cluster of eigenvalues appears
// in the leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The leading columns of \a VSL and \a VSR then form an orthonormal basis for the
// corresponding left and right eigenspaces (deflating subspaces).
//
// The parameter \a jobvsl specifies the computation of the left Schur vectors:
//
//   - \c 'V': The left Schur vectors of \a (A,B) are computed and returned in \a VSL.
//   - \c 'N': The left Schur vectors of \a (A,B) are not computed.
//
// The parameter \a jobvsr specifies the computation of the right Schur vectors:
//
//   - \c 'V': The right Schur vectors of \a (A,B) are computed and returned in \a VSR.
//   - \c 'N': The right Schur vectors of \a (A,B) are not computed.
//
// Via the parameter \a selctg it is possible to select specific eigenvalues. \a selctg is a
// pointer to a function of two single precision complex arguments returning an \c blas_int_t
// value. In case \a sort is set to \c 'N', \a selctg is not referenced, else if sort is set
// to \c 'S' \a selctg is used to select eigenvalues to sort to the top left of the Schur form.
// An eigenvalue \f$ alpha[j]/beta[j] \f$ is selected if

   \code
   selctg( &alpha[j], &beta[j])
   \endcode

// evaluates to true; i.e. if either one of a complex conjugate pair of eigenvalues is selected,
// then both complex eigenvalues are selected. On exit, the parameter \a sdim will contain the
// number of eigenvalues (after sorting) for which \a selctg is \a true. Note that complex
// conjugate pairs for which \a selctg is \a true for either eigenvalue count as 2.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 1,...,\a n: The QZ iteration failed. (A,B) are not in Schur form, but \a alphar(j),
//          \a alphai(j), and \a beta(j) should be correct for j=\a info,...,\a n-1.
//   - > \a n: see online reference for details.
//
// For more information on the sgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort,
                  blas_int_t (*selctg)( const complex<float>*, const complex<float>* ), blas_int_t n,
                  complex<float>* A, blas_int_t lda, complex<float>* B, blas_int_t ldb,
                  blas_int_t* sdim, complex<float>* alpha, complex<float>* beta,
                  complex<float>* VSL, blas_int_t ldvsl, complex<float>* VSR,
                  blas_int_t ldvsr, complex<float>* work, blas_int_t lwork,
                  float* rwork, blas_int_t* bwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
   using Selctg = MKL_INT (*)( const ET*, const ET* );
#else
   using ET = float;
   using Selctg = blas_int_t (*)( const ET*, const ET* );
#endif

   cgges_( &jobvsl, &jobvsr, &sort, reinterpret_cast<Selctg>( selctg ), &n,
           reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( B ), &ldb, sdim,
           reinterpret_cast<ET*>( alpha ), reinterpret_cast<ET*>( beta ),
           reinterpret_cast<ET*>( VSL ), &ldvsl, reinterpret_cast<ET*>( VSR ), &ldvsr,
           reinterpret_cast<ET*>( work ), &lwork, rwork, bwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the generalized Schur decomposition for a pair of
//        non-symmetric double precision complex column-major matrices.
// \ingroup lapack_eigenvalue
//
// \param jobvsl \c 'V' to compute the left Schur vectors, \c 'N' to not compute them.
// \param jobvsr \c 'V' to compute the right Schur vectors, \c 'N' to not compute them.
// \param sort \c 'S' to order the eigenvalues on the diagonal, \c 'N' to not order them (see \a selctg).
// \param selctg Logical function of two double precision complex arguments.
// \param n The order of the matrices \a A, \a B, \a VSL, and \a VSR \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix \a A.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param B Pointer to the first element of the double precision complex column-major matrix \a B.
// \param ldb The total number of elements between two columns of the matrix B; size >= max( 1, \a n ).
// \param sdim Returns the number of eigenvalues (after sorting) for which \a selctg is \a true.
// \param alpha Pointer to the first element of the eigenvalue numerator; size == n.
// \param beta Pointer to the first element of the eigenvalue denominator; size == n.
// \param VSL Pointer to the first element of the column-major matrix for the left Schur vectors.
// \param ldvsl The total number of elements between two columns of the matrix VSL \f$[0..\infty)\f$.
// \param VSR Pointer to the first element of the column-major matrix for the right Schur vectors.
// \param ldvsr The total number of elements between two columns of the matrix VSR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size == 8*n.
// \param bwork Auxiliary array; size == n.
// \param info Return code of the function call.
//
// This function computes for a pair of \a n-by-\a n double precision complex non-symmetric
// matrices (\a A,\a B) the generalized eigenvalues, the generalized real Schur form (\a S, \a T),
// and optionally the left and/or right matrices of Schur vectors (\a VSL and \a VSR). This gives
// the generalized Schur factorization

                      \f[ (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T ). \f]

// Optionally, it also orders the eigenvalues such that a selected cluster of eigenvalues appears
// in the leading diagonal blocks of the upper quasi-triangular matrix \a S and the upper triangular
// matrix \a T. The leading columns of \a VSL and \a VSR then form an orthonormal basis for the
// corresponding left and right eigenspaces (deflating subspaces).
//
// The parameter \a jobvsl specifies the computation of the left Schur vectors:
//
//   - \c 'V': The left Schur vectors of \a (A,B) are computed and returned in \a VSL.
//   - \c 'N': The left Schur vectors of \a (A,B) are not computed.
//
// The parameter \a jobvsr specifies the computation of the right Schur vectors:
//
//   - \c 'V': The right Schur vectors of \a (A,B) are computed and returned in \a VSR.
//   - \c 'N': The right Schur vectors of \a (A,B) are not computed.
//
// Via the parameter \a selctg it is possible to select specific eigenvalues. \a selctg is a
// pointer to a function of two double precision complex arguments returning an \c blas_int_t
// value. In case \a sort is set to \c 'N', \a selctg is not referenced, else if sort is set
// to \c 'S' \a selctg is used to select eigenvalues to sort to the top left of the Schur form.
// An eigenvalue \f$ alpha[j]/beta[j] \f$ is selected if

   \code
   selctg( &alpha[j], &beta[j])
   \endcode

// evaluates to true; i.e. if either one of a complex conjugate pair of eigenvalues is selected,
// then both complex eigenvalues are selected. On exit, the parameter \a sdim will contain the
// number of eigenvalues (after sorting) for which \a selctg is \a true. Note that complex
// conjugate pairs for which \a selctg is \a true for either eigenvalue count as 2.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 1,...,\a n: The QZ iteration failed. (A,B) are not in Schur form, but \a alphar(j),
//          \a alphai(j), and \a beta(j) should be correct for j=\a info,...,\a n-1.
//   - > \a n: see online reference for details.
//
// For more information on the sgges() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void gges( char jobvsl, char jobvsr, char sort,
                  blas_int_t (*selctg)( const complex<double>*, const complex<double>* ), blas_int_t n,
                  complex<double>* A, blas_int_t lda, complex<double>* B, blas_int_t ldb,
                  blas_int_t* sdim, complex<double>* alpha, complex<double>* beta,
                  complex<double>* VSL, blas_int_t ldvsl, complex<double>* VSR,
                  blas_int_t ldvsr, complex<double>* work, blas_int_t lwork,
                  double* rwork, blas_int_t* bwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
   using Selctg = MKL_INT (*)( const ET*, const ET* );
#else
   using ET = double;
   using Selctg = blas_int_t (*)( const ET*, const ET* );
#endif

   zgges_( &jobvsl, &jobvsr, &sort, reinterpret_cast<Selctg>( selctg ), &n,
           reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( B ), &ldb, sdim,
           reinterpret_cast<ET*>( alpha ), reinterpret_cast<ET*>( beta ),
           reinterpret_cast<ET*>( VSL ), &ldvsl, reinterpret_cast<ET*>( VSR ), &ldvsr,
           reinterpret_cast<ET*>( work ), &lwork, rwork, bwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************

} // namespace blaze

#endif
