//=================================================================================================
/*!
//  \file blaze/math/lapack/clapack/sytri.h
//  \brief Header file for the CLAPACK sytri wrapper functions
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

#ifndef _BLAZE_MATH_LAPACK_CLAPACK_SYTRI_H_
#define _BLAZE_MATH_LAPACK_CLAPACK_SYTRI_H_


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

void ssytri_( char* uplo, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, float* work, blaze::blas_int_t* info,
              blaze::fortran_charlen_t nuplo );
void dsytri_( char* uplo, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, double* work, blaze::blas_int_t* info,
              blaze::fortran_charlen_t nuplo );
void csytri_( char* uplo, blaze::blas_int_t* n, float* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, float* work, blaze::blas_int_t* info,
              blaze::fortran_charlen_t nuplo );
void zsytri_( char* uplo, blaze::blas_int_t* n, double* A, blaze::blas_int_t* lda,
              blaze::blas_int_t* ipiv, double* work, blaze::blas_int_t* info,
              blaze::fortran_charlen_t nuplo );

}
#endif
/*! \endcond */
//*************************************************************************************************




namespace blaze {

//=================================================================================================
//
//  LAPACK LDLT-BASED INVERSION FUNCTIONS (SYTRI)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK LDLT-based inversion functions (sytri) */
//@{
void sytri( char uplo, blas_int_t n, float* A, blas_int_t lda,
            const blas_int_t* ipiv, float* work, blas_int_t* info );

void sytri( char uplo, blas_int_t n, double* A, blas_int_t lda,
            const blas_int_t* ipiv, double* work, blas_int_t* info );

void sytri( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda,
            const blas_int_t* ipiv, complex<float>* work, blas_int_t* info );

void sytri( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda,
            const blas_int_t* ipiv, complex<double>* work, blas_int_t* info );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense symmetric indefinite single precision
//        column-major square matrix.
// \ingroup lapack_inversion
//
// \param uplo \c 'L' in case of a lower matrix, \c 'U' in case of an upper matrix.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array of size \a n.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK ssytri() function for
// symmetric indefinite single precision column-major matrices that have already been factorized
// by the ssytrf() function.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element D(i,i) is exactly zero and the inverse could not be computed.
//
// For more information on the ssytri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void sytri( char uplo, blas_int_t n, float* A, blas_int_t lda,
                   const blas_int_t* ipiv, float* work, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   ssytri_( &uplo, &n, A, &lda, const_cast<blas_int_t*>( ipiv ), work, info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense symmetric indefinite double precision
//        column-major square matrix.
// \ingroup lapack_inversion
//
// \param uplo \c 'L' in case of a lower matrix, \c 'U' in case of an upper matrix.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array of size \a n.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK dsytri() function for
// symmetric indefinite double precision column-major matrices that have already been factorized
// by the dsytrf() function.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element D(i,i) is exactly zero and the inverse could not be computed.
//
// For more information on the dsytri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void sytri( char uplo, blas_int_t n, double* A, blas_int_t lda,
                   const blas_int_t* ipiv, double* work, blas_int_t* info )
{
#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
#endif

   dsytri_( &uplo, &n, A, &lda, const_cast<blas_int_t*>( ipiv ), work, info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense symmetric indefinite single precision
//        complex column-major square matrix.
// \ingroup lapack_inversion
//
// \param uplo \c 'L' in case of a lower matrix, \c 'U' in case of an upper matrix.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the single precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array of size \a n.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK csytri() function for
// symmetric indefinite single precision complex column-major matrices that have already been
// factorized by the csytrf() function.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element D(i,i) is exactly zero and the inverse could not be computed.
//
// For more information on the csytri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void sytri( char uplo, blas_int_t n, complex<float>* A, blas_int_t lda,
                   const blas_int_t* ipiv, complex<float>* work, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<float> ) == 2UL*sizeof( float ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex8 ) == sizeof( complex<float> ) );
   using ET = MKL_Complex8;
#else
   using ET = float;
#endif

   csytri_( &uplo, &n, reinterpret_cast<ET*>( A ), &lda,
            const_cast<blas_int_t*>( ipiv ), reinterpret_cast<ET*>( work ), info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the inversion of the given dense symmetric indefinite double precision
//        complex column-major square matrix.
// \ingroup lapack_inversion
//
// \param uplo \c 'L' in case of a lower matrix, \c 'U' in case of an upper matrix.
// \param n The number of rows/columns of the symmetric matrix \f$[0..\infty)\f$.
// \param A Pointer to the first element of the double precision complex column-major matrix.
// \param lda The total number of elements between two columns of the matrix \f$[0..\infty)\f$.
// \param ipiv Auxiliary array of size \a n for the pivot indices.
// \param work Auxiliary array of size \a n.
// \param info Return code of the function call.
// \return void
//
// This function performs the dense matrix inversion based on the LAPACK zsytri() function for
// symmetric indefinite double precision complex column-major matrices that have already been
// factorized by the zsytrf() function.
//
// The \a info argument provides feedback on the success of the function call:
//
//   - = 0: The inversion finished successfully.
//   - < 0: If \a info = -i, the i-th argument had an illegal value.
//   - > 0: If \a info = i, element D(i,i) is exactly zero and the inverse could not be computed.
//
// For more information on the zsytri() function, see the LAPACK online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
inline void sytri( char uplo, blas_int_t n, complex<double>* A, blas_int_t lda,
                   const blas_int_t* ipiv, complex<double>* work, blas_int_t* info )
{
   BLAZE_STATIC_ASSERT( sizeof( complex<double> ) == 2UL*sizeof( double ) );

#if defined(INTEL_MKL_VERSION)
   BLAZE_STATIC_ASSERT( sizeof( MKL_INT ) == sizeof( blas_int_t ) );
   BLAZE_STATIC_ASSERT( sizeof( MKL_Complex16 ) == sizeof( complex<double> ) );
   using ET = MKL_Complex16;
#else
   using ET = double;
#endif

   zsytri_( &uplo, &n, reinterpret_cast<ET*>( A ), &lda,
            const_cast<blas_int_t*>( ipiv ), reinterpret_cast<ET*>( work ), info
#if !defined(INTEL_MKL_VERSION)
          , blaze::fortran_charlen_t(1)
#endif
          );
}
//*************************************************************************************************

} // namespace blaze

#endif
