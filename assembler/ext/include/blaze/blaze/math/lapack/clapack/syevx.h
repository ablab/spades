//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/syevx.h
//  \brief Header file for the CLAPACK syevx wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_SYEVX_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_SYEVX_H_


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

void ssyevx_( char* jobz, char* range, char* uplo, blaze::blas_int_t* n, float* A,
              blaze::blas_int_t* lda, float* vl, float* vu, blaze::blas_int_t* il,
              blaze::blas_int_t* iu, float* abstol, blaze::blas_int_t* m, float* w, float* Z,
              blaze::blas_int_t* ldz, float* work, blaze::blas_int_t* lwork,
              blaze::blas_int_t* iwork, blaze::blas_int_t* ifail, blaze::blas_int_t* info,
              blaze::fortran_charlen_t njobz, blaze::fortran_charlen_t nrange,
              blaze::fortran_charlen_t nuplo );
void dsyevx_( char* jobz, char* range, char* uplo, blaze::blas_int_t* n, double* A,
              blaze::blas_int_t* lda, double* vl, double* vu, blaze::blas_int_t* il,
              blaze::blas_int_t* iu, double* abstol, blaze::blas_int_t* m, double* w, double* Z,
              blaze::blas_int_t* ldz, double* work, blaze::blas_int_t* lwork,
              blaze::blas_int_t* iwork, blaze::blas_int_t* ifail, blaze::blas_int_t* info,
              blaze::fortran_charlen_t njobz, blaze::fortran_charlen_t nrange,
              blaze::fortran_charlen_t nuplo );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK SYMMETRIC MATRIX EIGENVALUE FUNCTIONS (SYEVX)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK symmetric matrix eigenvalue functions (syevx) */
//@{
void syevx( char jobz, char range, char uplo, blas_int_t n, float* A,
            blas_int_t lda, float vl, float vu, blas_int_t il,
            blas_int_t iu, float abstol, blas_int_t* m, float* w,
            float* Z, blas_int_t ldz, float* work, blas_int_t lwork,
            blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );

void syevx( char jobz, char range, char uplo, blas_int_t n, double* A,
            blas_int_t lda, double vl, double vu, blas_int_t il,
            blas_int_t iu, double abstol, blas_int_t* m, double* w,
            double* Z, blas_int_t ldz, double* work, blas_int_t lwork,
            blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense symmetric single
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobz \c 'V' to compute the eigenvectors of \a A, \c 'N' to only compute the eigenvalues.
// \param range Specifies the range of eigenvalues to find (\c 'A', \c 'V', or \c 'I').
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A \f$[0..\infty)\f$.
// \param vl The lower bound of the interval to be searched for eigenvalues (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for eigenvalues (\a vl < \a vu).
// \param il The index of the smallest eigenvalue to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest eigenvalue to be returned (0 <= \a il <= \a iu).
// \param abstol The absolute error tolerance for the computation of the eigenvalues.
// \param m The total number of eigenvalues found (0 <= \a m <= \a n).
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param Z Pointer to the first element of the column-major matrix for the eigenvectors.
// \param ldz The total number of elements between two columns of the matrix Z \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param iwork Auxiliary array; size >= 5*\a n.
// \param ifail Index array of eigenvectors that failed to converge.
// \param info Return code of the function call.
// \return void
//
// This function computes a specified number of eigenvalues of a symmetric \a n-by-\a n single
// precision column-major matrix based on the LAPACK ssyevx() function. Optionally, it computes
// a specified number of eigenvectors. The selected real eigenvalues are returned in ascending
// order in the given \a n-dimensional vector \a w.
//
// The parameter \a jobz specifies the computation of the orthonormal eigenvectors of \a A:
//
//   - \c 'V': The eigenvectors of \a A are computed and returned within the matrix \a Z.
//   - \c 'N': The eigenvectors of \a A are not computed.
//
// The parameter \a range specifies the amount of eigenvalues/eigenvectors to be found:
//
//   - \c 'A': All eigenvalues will be found.
//   - \c 'V': All eigenvalues in the half-open interval \f$(vl..vu]\f$ will be found.
//   - \c 'I': The \a il-th through \a iu-th eigenvalues will be found.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the algorithm failed to converge; i values did not converge to zero.
//
// For more information on the ssyevx() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void syevx( char jobz, char range, char uplo, blas_int_t n, float* A,
                   blas_int_t lda, float vl, float vu, blas_int_t il,
                   blas_int_t iu, float abstol, blas_int_t* m, float* w,
                   float* Z, blas_int_t ldz, float* work, blas_int_t lwork,
                   blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   ++il;
   ++iu;

   ssyevx_( &jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu,
            &abstol, m, w, Z, &ldz, work, &lwork, iwork, ifail, info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense symmetric double
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobz \c 'V' to compute the eigenvectors of \a A, \c 'N' to only compute the eigenvalues.
// \param range Specifies the range of eigenvalues to find (\c 'A', \c 'V', or \c 'I').
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A \f$[0..\infty)\f$.
// \param vl The lower bound of the interval to be searched for eigenvalues (\a vl < \a vu).
// \param vu The upper bound of the interval to be searched for eigenvalues (\a vl < \a vu).
// \param il The index of the smallest eigenvalue to be returned (0 <= \a il <= \a iu).
// \param iu The index of the largest eigenvalue to be returned (0 <= \a il <= \a iu).
// \param abstol The absolute error tolerance for the computation of the eigenvalues.
// \param m The total number of eigenvalues found (0 <= \a m <= \a n).
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param Z Pointer to the first element of the column-major matrix for the eigenvectors.
// \param ldz The total number of elements between two columns of the matrix Z \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param iwork Auxiliary array; size >= 5*\a n.
// \param ifail Index array of eigenvectors that failed to converge.
// \param info Return code of the function call.
// \return void
//
// This function computes a specified number of eigenvalues of a symmetric \a n-by-\a n double
// precision column-major matrix based on the LAPACK dsyevx() function. Optionally, it computes
// a specified number of eigenvectors. The selected real eigenvalues are returned in ascending
// order in the given \a n-dimensional vector \a w.
//
// The parameter \a jobz specifies the computation of the orthonormal eigenvectors of \a A:
//
//   - \c 'V': The eigenvectors of \a A are computed and returned within the matrix \a Z.
//   - \c 'N': The eigenvectors of \a A are not computed.
//
// The parameter \a range specifies the amount of eigenvalues/eigenvectors to be found:
//
//   - \c 'A': All eigenvalues will be found.
//   - \c 'V': All eigenvalues in the half-open interval \f$(vl..vu]\f$ will be found.
//   - \c 'I': The \a il-th through \a iu-th eigenvalues will be found.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the algorithm failed to converge; i values did not converge to zero.
//
// For more information on the dsyevx() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void syevx( char jobz, char range, char uplo, blas_int_t n, double* A,
                   blas_int_t lda, double vl, double vu, blas_int_t il,
                   blas_int_t iu, double abstol, blas_int_t* m, double* w,
                   double* Z, blas_int_t ldz, double* work, blas_int_t lwork,
                   blas_int_t* iwork, blas_int_t* ifail, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   ++il;
   ++iu;

   dsyevx_( &jobz, &range, &uplo, &n, A, &lda, &vl, &vu, &il, &iu,
            &abstol, m, w, Z, &ldz, work, &lwork, iwork, ifail, info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************

} // namespace blaze

#endif
