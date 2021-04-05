//=================================================================================================
/*!
//  \file blaze/math/blas/trmv.h
//  \brief Header file for BLAS triangular matrix/vector multiplication functions (trmv)
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

#ifndef _BLAZE_MATH_BLAS_TRMV_H_
#define _BLAZE_MATH_BLAS_TRMV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/blas/cblas/trmv.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/Assert.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  BLAS TRIANGULAR MATRIX/VECTOR MULTIPLICATION FUNCTIONS (TRMV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS triangular matrix/vector multiplication functions (trmv) */
//@{
#if BLAZE_BLAS_MODE

template< typename VT, typename MT, bool SO >
void trmv( DenseVector<VT,false>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo );

template< typename VT, typename MT, bool SO >
void trmv( DenseVector<VT,true>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense vector multiplication
//        (\f$ \vec{x}=A*\vec{x} \f$).
// \ingroup blas
//
// \param x The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a triangular matrix by a vector based on the BLAS
// trmv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
inline void trmv( DenseVector<VT,false>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).rows() )    );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );

   trmv( ( IsRowMajorMatrix_v<MT> )?( CblasRowMajor ):( CblasColMajor ),
         uplo, CblasNoTrans, CblasNonUnit, n, (*A).data(), lda, (*x).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/triangular dense matrix multiplication
//        (\f$ \vec{x}^T=\vec{x}^T*A \f$).
// \ingroup blas
//
// \param x The target left-hand side dense vector.
// \param A The dense matrix operand.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \return void
//
// This function performs the multiplication of a vector and a triangular matrix based on the BLAS
// trmv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename VT  // Type of the target vector
        , typename MT  // Type of the matrix operand
        , bool SO >    // Storage order of the matrix operand
inline void trmv( DenseVector<VT,true>& x, const DenseMatrix<MT,SO>& A, CBLAS_UPLO uplo )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).rows() )    );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );

   trmv( ( IsRowMajorMatrix_v<MT> )?( CblasRowMajor ):( CblasColMajor ),
         uplo, CblasTrans, CblasNonUnit, n, (*A).data(), lda, (*x).data(), 1 );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
