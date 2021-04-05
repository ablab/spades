//=================================================================================================
/*!
//  \file blaze/math/blas/gemm.h
//  \brief Header file for BLAS general matrix/matrix multiplication functions (gemm)
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

#ifndef _BLAZE_MATH_BLAS_GEMM_H_
#define _BLAZE_MATH_BLAS_GEMM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/blas/cblas/gemm.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/NumericCast.h>


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

template< typename MT1, bool SO1, typename MT2, bool SO2, typename MT3, bool SO3, typename ST >
void gemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
           const DenseMatrix<MT3,SO3>& B, ST alpha, ST beta );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense matrix multiplication (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup blas
//
// \param C The target left-hand side dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function performs the dense matrix/dense matrix multiplication based on the BLAS
// gemm() functions. Note that the function only works for matrices with \c float, \c double,
// \c complex<float>, and \c complex<double> element type. The attempt to call the function
// with matrices of any other element type results in a compile time error.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1   // Type of the left-hand side target matrix
        , bool SO1       // Storage order of the left-hand side target matrix
        , typename MT2   // Type of the left-hand side matrix operand
        , bool SO2       // Storage order of the left-hand side matrix operand
        , typename MT3   // Type of the right-hand side matrix operand
        , bool SO3       // Storage order of the right-hand side matrix operand
        , typename ST >  // Type of the scalar factors
inline void gemm( DenseMatrix<MT1,SO1>& C, const DenseMatrix<MT2,SO2>& A,
                  const DenseMatrix<MT3,SO3>& B, ST alpha, ST beta )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   const blas_int_t m  ( numeric_cast<blas_int_t>( (*A).rows() )    );
   const blas_int_t n  ( numeric_cast<blas_int_t>( (*B).columns() ) );
   const blas_int_t k  ( numeric_cast<blas_int_t>( (*A).columns() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb( numeric_cast<blas_int_t>( (*B).spacing() ) );
   const blas_int_t ldc( numeric_cast<blas_int_t>( (*C).spacing() ) );

   gemm( ( IsRowMajorMatrix_v<MT1> )?( CblasRowMajor ):( CblasColMajor ),
         ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
         ( SO1 == SO3 )?( CblasNoTrans ):( CblasTrans ),
         m, n, k, alpha, (*A).data(), lda, (*B).data(), ldb, beta, (*C).data(), ldc );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
