//=================================================================================================
/*!
//  \file blaze/math/blas/trmm.h
//  \brief Header file for BLAS triangular matrix/matrix multiplication functions (trmm)
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

#ifndef _BLAZE_MATH_BLAS_TRMM_H_
#define _BLAZE_MATH_BLAS_TRMM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/blas/cblas/trmm.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/BLAS.h>
#include <blaze/util/Assert.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  BLAS TRIANGULAR MATRIX MULTIPLICATION FUNCTIONS (TRMM)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS triangular matrix multiplication functions (trmm) */
//@{
#if BLAZE_BLAS_MODE

template< typename MT1, bool SO1, typename MT2, bool SO2, typename ST >
void trmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
           CBLAS_SIDE side, CBLAS_UPLO uplo, ST alpha );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a triangular dense matrix/dense matrix multiplication
//        (\f$ B=\alpha*A*B \f$ or \f$ B=\alpha*B*A \f$).
// \ingroup blas
//
// \param B The target dense matrix.
// \param A The dense matrix multiplication operand.
// \param side \a CblasLeft to compute \f$ B=\alpha*A*B \f$, \a CblasRight to compute \f$ B=\alpha*B*A \f$.
// \param uplo \a CblasLower to use the lower triangle from \a A, \a CblasUpper to use the upper triangle.
// \param alpha The scaling factor for \f$ A*B \f$ or \f$ B*A \f$.
// \return void
//
// This function performs the scaling and multiplication of a triangular matrix by a matrix
// based on the BLAS trmm() functions. Note that the function only works for matrices with
// \c float, \c double, \c complex<float>, and \c complex<double> element type. The attempt to
// call the function with matrices of any other element type results in a compile time error.
// Also note that matrix \a A is expected to be a square matrix.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT1   // Type of the target matrix
        , bool SO1       // Storage order of the target matrix
        , typename MT2   // Type of the left-hand side matrix operand
        , bool SO2       // Storage order of the left-hand side matrix operand
        , typename ST >  // Type of the scalar factor
inline void trmm( DenseMatrix<MT1,SO1>& B, const DenseMatrix<MT2,SO2>& A,
                  CBLAS_SIDE side, CBLAS_UPLO uplo, ST alpha )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*A).columns(), "Non-square triangular matrix detected" );
   BLAZE_INTERNAL_ASSERT( side == CblasLeft  || side == CblasRight, "Invalid side argument detected" );
   BLAZE_INTERNAL_ASSERT( uplo == CblasLower || uplo == CblasUpper, "Invalid uplo argument detected" );

   const blas_int_t m  ( numeric_cast<blas_int_t>( (*B).rows() )    );
   const blas_int_t n  ( numeric_cast<blas_int_t>( (*B).columns() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );
   const blas_int_t ldb( numeric_cast<blas_int_t>( (*B).spacing() ) );

   trmm( ( IsRowMajorMatrix_v<MT1> )?( CblasRowMajor ):( CblasColMajor ),
         side,
         ( SO1 == SO2 )?( uplo ):( ( uplo == CblasLower )?( CblasUpper ):( CblasLower ) ),
         ( SO1 == SO2 )?( CblasNoTrans ):( CblasTrans ),
         CblasNonUnit,
         m, n, alpha, (*A).data(), lda, (*B).data(), ldb );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
