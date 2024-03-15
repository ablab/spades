//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/heev.h
//  \brief Header file for the CLAPACK heev wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_HEEV_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_HEEV_H_


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

void cheev_( char* jobz, char* uplo, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
             float* w, float* work, blaze::blas_int_t* lwork, float* rwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobz, blaze::fortran_charlen_t nuplo );
void zheev_( char* jobz, char* uplo, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
             double* w, double* work, blaze::blas_int_t* lwork, double* rwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobz, blaze::fortran_charlen_t nuplo );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK HERMITIAN MATRIX EIGENVALUE FUNCTIONS (HEEV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK Hermitian matrix eigenvalue functions (heev) */
//@{
void heev( char jobz, char uplo, blas_int_t n, complex<float>* A,
           blas_int_t lda, float* w, complex<float>* work,
           blas_int_t lwork, float* rwork, blas_int_t* info );

void heev( char jobz, char uplo, blas_int_t n, complex<double>* A,
           blas_int_t lda, double* w, complex<double>* work,
           blas_int_t lwork, double* rwork, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense Hermitian single
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobz \c 'V' to compute the eigenvectors of \a A, \c 'N' to only compute the eigenvalues.
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A \f$[0..\infty)\f$.
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size >= max(1,3*\a n-2).
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of an Hermitian \a n-by-\a n single precision complex
// column-major matrix based on the LAPACK cheev() function. Optionally, it computes the left
// and right eigenvectors. The real eigenvalues are returned in ascending order in the given
// \a n-dimensional vector \a w.
//
// The parameter \a jobz specifies the computation of the orthonormal eigenvectors of \a A:
//
//   - \c 'V': The eigenvectors of \a A are computed and returned within the matrix \a A.
//   - \c 'N': The eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the algorithm failed to converge; i values did not converge to zero.
//
// For more information on the cheev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void heev( char jobz, char uplo, blas_int_t n, complex<float>* A,
                  blas_int_t lda, float* w, complex<float>* work,
                  blas_int_t lwork, float* rwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cheev_( &jobz, &uplo, &n, reinterpret_cast<ET*>( A ), &lda, w,
           reinterpret_cast<ET*>( work ), &lwork, rwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense Hermitian double
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobz \c 'V' to compute the eigenvectors of \a A, \c 'N' to only compute the eigenvalues.
// \param uplo \c 'L' to use the lower part of the matrix, \c 'U' to use the upper part.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A \f$[0..\infty)\f$.
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size >= max(1,3*\a n-2).
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of an Hermitian \a n-by-\a n double precision complex
// column-major matrix based on the LAPACK zheev() function. Optionally, it computes the left
// and right eigenvectors. The real eigenvalues are returned in ascending order in the given
// \a n-dimensional vector \a w.
//
// The parameter \a jobz specifies the computation of the orthonormal eigenvectors of \a A:
//
//   - \c 'V': The eigenvectors of \a A are computed and returned within the matrix \a A.
//   - \c 'N': The eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the algorithm failed to converge; i values did not converge to zero.
//
// For more information on the zheev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void heev( char jobz, char uplo, blas_int_t n, complex<double>* A,
                  blas_int_t lda, double* w, complex<double>* work,
                  blas_int_t lwork, double* rwork, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zheev_( &jobz, &uplo, &n, reinterpret_cast<ET*>( A ), &lda, w,
           reinterpret_cast<ET*>( work ), &lwork, rwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************

} // namespace blaze

#endif
