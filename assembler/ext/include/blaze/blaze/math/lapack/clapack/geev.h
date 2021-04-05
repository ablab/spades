//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/geev.h
//  \brief Header file for the CLAPACK geev wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_GEEV_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_GEEV_H_


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

void sgeev_( char* jobvl, char* jobvr, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
             float* wr, float* wi, float* VL, blaze::blas_int_t* ldvl, float* VR, blaze::blas_int_t* ldvr,
             float* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvl, blaze::fortran_charlen_t njobvr );
void dgeev_( char* jobvl, char* jobvr, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
             double* wr, double* wi, double* VL, blaze::blas_int_t* ldvl, double* VR, blaze::blas_int_t* ldvr,
             double* work, blaze::blas_int_t* lwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvl, blaze::fortran_charlen_t njobvr );
void cgeev_( char* jobvl, char* jobvr, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
             float* w, float* VL, blaze::blas_int_t* ldvl, float* VR, blaze::blas_int_t* ldvr,
             float* work, blaze::blas_int_t* lwork, float* rwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvl, blaze::fortran_charlen_t njobvr );
void zgeev_( char* jobvl, char* jobvr, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
             double* w, double* VL, blaze::blas_int_t* ldvl, double* VR, blaze::blas_int_t* ldvr,
             double* work, blaze::blas_int_t* lwork, double* rwork, blaze::blas_int_t* info,
             blaze::fortran_charlen_t njobvl, blaze::fortran_charlen_t njobvr );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK GENERAL MATRIX EIGENVALUE FUNCTIONS (GEEV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK general matrix eigenvalue functions (geev) */
//@{
void geev( char jobvl, char jobvr, blas_int_t n, float* A, blas_int_t lda,
           float* wr, float* wi, float* VL, blas_int_t ldvl, float* VR, blas_int_t ldvr,
           float* work, blas_int_t lwork, blas_int_t* info );

void geev( char jobvl, char jobvr, blas_int_t n, double* A, blas_int_t lda,
           double* wr, double* wi, double* VL, blas_int_t ldvl, double* VR, blas_int_t ldvr,
           double* work, blas_int_t lwork, blas_int_t* info );

void geev( char jobvl, char jobvr, blas_int_t n, complex<float>* A, blas_int_t lda,
           complex<float>* w, complex<float>* VL, blas_int_t ldvl, complex<float>* VR,
           blas_int_t ldvr, complex<float>* work, blas_int_t lwork, float* rwork,
           blas_int_t* info );

void geev( char jobvl, char jobvr, blas_int_t n, complex<double>* A, blas_int_t lda,
           complex<double>* w, complex<double>* VL, blas_int_t ldvl, complex<double>* VR,
           blas_int_t ldvr, complex<double>* work, blas_int_t lwork, double* rwork,
           blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general single
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobvl \c 'V' to compute the left eigenvectors of \a A, \c 'N' to not compute them.
// \param jobvr \c 'V' to compute the right eigenvectors of \a A, \c 'N' to not compute them.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param wr Pointer to the first element of the vector for the real part of the eigenvalues.
// \param wi Pointer to the first element of the vector for the imaginary part of the eigenvalues.
// \param VL Pointer to the first element of the column-major matrix for the left eigenvectors.
// \param ldvl The total number of elements between two columns of the matrix VL \f$[0..\infty)\f$.
// \param VR Pointer to the first element of the column-major matrix for the right eigenvectors.
// \param ldvr The total number of elements between two columns of the matrix VR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of a non-symmetric \a n-by-\a n single precision
// column-major matrix based on the LAPACK sgeev() function. Optionally, it computes the left
// and/or right eigenvectors. The right eigenvector \f$v[j]\f$ of \a satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the
// positive imaginary part first. The computed eigenvectors are normalized to have Euclidean
// norm equal to 1 and largest component real.
//
// The parameter \a jobvl specifies the computation of the left eigenvectors:
//
//   - \c 'V': The left eigenvectors of \a A are computed and returned in \a VL.
//   - \c 'N': The left eigenvectors of \a A are not computed.
//
// The parameter \a jobvr specifies the computation of the right eigenvectors:
//
//   - \c 'V': The right eigenvectors of \a A are computed and returned in \a VR.
//   - \c 'N': The right eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the QR algorithm failed to compute all the eigenvalues and no eigenvectors
//          have been computed; the elements with index larger than \a i have converged.
//
// For more information on the sgeev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geev( char jobvl, char jobvr, blas_int_t n, float* A, blas_int_t lda,
                  float* wr, float* wi, float* VL, blas_int_t ldvl, float* VR, blas_int_t ldvr,
                  float* work, blas_int_t lwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   sgeev_( &jobvl, &jobvr, &n, A, &lda, wr, wi, VL, &ldvl, VR, &ldvr, work, &lwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general double
//        precision column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobvl \c 'V' to compute the left eigenvectors of \a A, \c 'N' to not compute them.
// \param jobvr \c 'V' to compute the right eigenvectors of \a A, \c 'N' to not compute them.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param wr Pointer to the first element of the vector for the real part of the eigenvalues.
// \param wi Pointer to the first element of the vector for the imaginary part of the eigenvalues.
// \param VL Pointer to the first element of the column-major matrix for the left eigenvectors.
// \param ldvl The total number of elements between two columns of the matrix VL \f$[0..\infty)\f$.
// \param VR Pointer to the first element of the column-major matrix for the right eigenvectors.
// \param ldvr The total number of elements between two columns of the matrix VR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of a non-symmetric \a n-by-\a n double precision
// column-major matrix based on the LAPACK dgeev() function. Optionally, it computes the left
// and/or right eigenvectors. The right eigenvector \f$v[j]\f$ of \a satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the
// positive imaginary part first. The computed eigenvectors are normalized to have Euclidean
// norm equal to 1 and largest component real.
//
// The parameter \a jobvl specifies the computation of the left eigenvectors:
//
//   - \c 'V': The left eigenvectors of \a A are computed and returned in \a VL.
//   - \c 'N': The left eigenvectors of \a A are not computed.
//
// The parameter \a jobvr specifies the computation of the right eigenvectors:
//
//   - \c 'V': The right eigenvectors of \a A are computed and returned in \a VR.
//   - \c 'N': The right eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the QR algorithm failed to compute all the eigenvalues and no eigenvectors
//          have been computed; the elements with index larger than \a i have converged.
//
// For more information on the dgeev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geev( char jobvl, char jobvr, blas_int_t n, double* A, blas_int_t lda,
                  double* wr, double* wi, double* VL, blas_int_t ldvl, double* VR, blas_int_t ldvr,
                  double* work, blas_int_t lwork, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dgeev_( &jobvl, &jobvr, &n, A, &lda, wr, wi, VL, &ldvl, VR, &ldvr, work, &lwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general single
//        precision complex column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobvl \c 'V' to compute the left eigenvectors of \a A, \c 'N' to not compute them.
// \param jobvr \c 'V' to compute the right eigenvectors of \a A, \c 'N' to not compute them.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param VL Pointer to the first element of the column-major matrix for the left eigenvectors.
// \param ldvl The total number of elements between two columns of the matrix VL \f$[0..\infty)\f$.
// \param VR Pointer to the first element of the column-major matrix for the right eigenvectors.
// \param ldvr The total number of elements between two columns of the matrix VR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size >= 2*\a n.
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of a non-symmetric \a n-by-\a n single precision
// complex column-major matrix based on the LAPACK cgeev() function. Optionally, it computes
// the left and/or right eigenvectors. The right eigenvector \f$v[j]\f$ of \a satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the
// positive imaginary part first. The computed eigenvectors are normalized to have Euclidean
// norm equal to 1 and largest component real.
//
// The parameter \a jobvl specifies the computation of the left eigenvectors:
//
//   - \c 'V': The left eigenvectors of \a A are computed and returned in \a VL.
//   - \c 'N': The left eigenvectors of \a A are not computed.
//
// The parameter \a jobvr specifies the computation of the right eigenvectors:
//
//   - \c 'V': The right eigenvectors of \a A are computed and returned in \a VR.
//   - \c 'N': The right eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the QR algorithm failed to compute all the eigenvalues and no eigenvectors
//          have been computed; the elements with index larger than \a i have converged.
//
// For more information on the cgeev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geev( char jobvl, char jobvr, blas_int_t n, complex<float>* A, blas_int_t lda,
                  complex<float>* w, complex<float>* VL, blas_int_t ldvl, complex<float>* VR,
                  blas_int_t ldvr, complex<float>* work, blas_int_t lwork, float* rwork,
                  blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   cgeev_( &jobvl, &jobvr, &n, reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( w ),
           reinterpret_cast<ET*>( VL ), &ldvl, reinterpret_cast<ET*>( VR ), &ldvr,
           reinterpret_cast<ET*>( work ), &lwork, rwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for computing the eigenvalues of the given dense general double
//        precision complex column-major matrix.
// \ingroup lapack_eigenvalue
//
// \param jobvl \c 'V' to compute the left eigenvectors of \a A, \c 'N' to not compute them.
// \param jobvr \c 'V' to compute the right eigenvectors of \a A, \c 'N' to not compute them.
// \param n The number of rows and columns of the given matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix A; size >= max( 1, \a n ).
// \param w Pointer to the first element of the vector for the eigenvalues.
// \param VL Pointer to the first element of the column-major matrix for the left eigenvectors.
// \param ldvl The total number of elements between two columns of the matrix VL \f$[0..\infty)\f$.
// \param VR Pointer to the first element of the column-major matrix for the right eigenvectors.
// \param ldvr The total number of elements between two columns of the matrix VR \f$[0..\infty)\f$.
// \param work Auxiliary array; size >= max( 1, \a lwork ).
// \param lwork The dimension of the array \a work; see online reference for details.
// \param rwork Auxiliary array; size >= 2*\a n.
// \param info Return code of the function call.
// \return void
//
// This function computes the eigenvalues of a non-symmetric \a n-by-\a n double precision
// complex column-major matrix based on the LAPACK zgeev() function. Optionally, it computes
// the left and/or right eigenvectors. The right eigenvector \f$v[j]\f$ of \a satisfies

                          \f[ A * v[j] = lambda[j] * v[j], \f]

// where \f$lambda[j]\f$ is its eigenvalue. The left eigenvector \f$u[j]\f$ of \a A satisfies

                       \f[ u[j]^{H} * A = lambda[j] * u[j]^{H}, \f]

// where \f$u[j]^{H}\f$ denotes the conjugate transpose of \f$u[j]\f$.
//
// Complex conjugate pairs of eigenvalues appear consecutively with the eigenvalue having the
// positive imaginary part first. The computed eigenvectors are normalized to have Euclidean
// norm equal to 1 and largest component real.
//
// The parameter \a jobvl specifies the computation of the left eigenvectors:
//
//   - \c 'V': The left eigenvectors of \a A are computed and returned in \a VL.
//   - \c 'N': The left eigenvectors of \a A are not computed.
//
// The parameter \a jobvr specifies the computation of the right eigenvectors:
//
//   - \c 'V': The right eigenvectors of \a A are computed and returned in \a VR.
//   - \c 'N': The right eigenvectors of \a A are not computed.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The computation finished successfully.
//   - < 0: If info = -i, the i-th argument had an illegal value.
//   - > 0: If info = i, the QR algorithm failed to compute all the eigenvalues and no eigenvectors
//          have been computed; the elements with index larger than \a i have converged.
//
// For more information on the zgeev() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void geev( char jobvl, char jobvr, blas_int_t n, complex<double>* A, blas_int_t lda,
                  complex<double>* w, complex<double>* VL, blas_int_t ldvl, complex<double>* VR,
                  blas_int_t ldvr, complex<double>* work, blas_int_t lwork, double* rwork,
                  blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zgeev_( &jobvl, &jobvr, &n, reinterpret_cast<ET*>( A ), &lda, reinterpret_cast<ET*>( w ),
           reinterpret_cast<ET*>( VL ), &ldvl, reinterpret_cast<ET*>( VR ), &ldvr,
           reinterpret_cast<ET*>( work ), &lwork, rwork, info
#if !defined(INTEL_MKL_VERSION)
         , blaze::fortran_charlen_t(1), blaze::fortran_charlen_t(1)
#endif
         );
}
//*************************************************************************************************

} // namespace blaze

#endif
