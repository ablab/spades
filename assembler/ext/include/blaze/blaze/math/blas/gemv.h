//=================================================================================================
/*!
//  \file blaze/math/blas/gemv.h
//  \brief Header file for BLAS general matrix/vector multiplication functions (gemv)
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

#ifndef _BLAZE_MATH_BLAS_GEMV_H_
#define _BLAZE_MATH_BLAS_GEMV_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/blas/cblas/gemv.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/ConstDataAccess.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  BLAS GENERAL MATRIX/VECTOR MULTIPLICATION FUNCTIONS (GEMV)
//
//=================================================================================================

//*************************************************************************************************
/*!\name BLAS general matrix/vector multiplication functions (gemv) */
//@{
#if BLAZE_BLAS_MODE

template< typename VT1, typename MT1, bool SO, typename VT2, typename ST >
void gemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
           const DenseVector<VT2,false>& x, ST alpha, ST beta );

template< typename VT1, typename VT2, typename MT1, bool SO, typename ST >
void gemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
           const DenseMatrix<MT1,SO>& A, ST alpha, ST beta );

#endif
//@}
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a dense matrix/dense vector multiplication
//        (\f$ \vec{y}=\alpha*A*\vec{x}+\beta*\vec{y} \f$).
// \ingroup blas
//
// \param y The target left-hand side dense vector.
// \param A The left-hand side dense matrix operand.
// \param x The right-hand side dense vector operand.
// \param alpha The scaling factor for \f$ A*\vec{x} \f$.
// \param beta The scaling factor for \f$ \vec{y} \f$.
// \return void
//
// This function performs the dense matrix/dense vector multiplication based on the BLAS gemv()
// functions. Note that the function only works for vectors and matrices with \c float, \c double,
// \c complex<float>, or \c complex<double> element type. The attempt to call the function with
// vectors and matrices of any other element type results in a compile time error.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename VT1   // Type of the left-hand side target vector
        , typename MT1   // Type of the left-hand side matrix operand
        , bool SO        // Storage order of the left-hand side matrix operand
        , typename VT2   // Type of the right-hand side vector operand
        , typename ST >  // Type of the scalar factors
inline void gemv( DenseVector<VT1,false>& y, const DenseMatrix<MT1,SO>& A,
                  const DenseVector<VT2,false>& x, ST alpha, ST beta )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   const blas_int_t m  ( numeric_cast<blas_int_t>( (*A).rows() )    );
   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).columns() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );

   gemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasNoTrans, m, n, alpha,
         (*A).data(), lda, (*x).data(), 1, beta, (*y).data(), 1 );
}
#endif
//*************************************************************************************************


//*************************************************************************************************
#if BLAZE_BLAS_MODE
/*!\brief BLAS kernel for a transpose dense vector/dense matrix multiplication
//        (\f$ \vec{y}^T=\alpha*\vec{x}^T*A+\beta*\vec{y}^T \f$).
// \ingroup blas
//
// \param y The target left-hand side dense vector.
// \param x The left-hand side dense vector operand.
// \param A The right-hand side dense matrix operand.
// \param alpha The scaling factor for \f$ \vec{x}^T*A \f$.
// \param beta The scaling factor for \f$ \vec{y}^T \f$.
// \return void
//
// This function performs the transpose dense vector/dense matrix multiplication based on the
// BLAS gemv() functions. Note that the function only works for vectors and matrices with \c float,
// \c double, \c complex<float>, or \c complex<double> element type. The attempt to call the
// function with vectors and matrices of any other element type results in a compile time error.
//
// \note This function can only be used if a fitting BLAS library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename VT1   // Type of the left-hand side target vector
        , typename VT2   // Type of the left-hand side vector operand
        , typename MT1   // Type of the right-hand side matrix operand
        , bool SO        // Storage order of the right-hand side matrix operand
        , typename ST >  // Type of the scalar factors
inline void gemv( DenseVector<VT1,true>& y, const DenseVector<VT2,true>& x,
                  const DenseMatrix<MT1,SO>& A, ST alpha, ST beta )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT2 );

   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( VT1 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( VT2 );
   BLAZE_CONSTRAINT_MUST_HAVE_CONST_DATA_ACCESS  ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   const blas_int_t m  ( numeric_cast<blas_int_t>( (*A).rows() )    );
   const blas_int_t n  ( numeric_cast<blas_int_t>( (*A).columns() ) );
   const blas_int_t lda( numeric_cast<blas_int_t>( (*A).spacing() ) );

   gemv( ( SO )?( CblasColMajor ):( CblasRowMajor ), CblasTrans, m, n, alpha,
         (*A).data(), lda, (*x).data(), 1, beta, (*y).data(), 1 );
}
#endif
//*************************************************************************************************

} // namespace blaze

#endif
