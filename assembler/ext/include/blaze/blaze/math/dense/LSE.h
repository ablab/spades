//=================================================================================================
/*!
//  \file blaze/math/dense/LSE.h
//  \brief Header file for the dense matrix linear system solver kernels
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

#ifndef _BLAZE_MATH_DENSE_LSE_H_
#define _BLAZE_MATH_DENSE_LSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/constraints/Uniform.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/lapack/gesv.h>
#include <blaze/math/lapack/hesv.h>
#include <blaze/math/lapack/sysv.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsGeneral.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 0x0 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 0 \times 0 \f$ linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 0 \times 0 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 0 \times 0 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT1      // Type of the system matrix
        , bool SO1          // Storage order of the system matrix
        , typename MT2      // Type of the solution matrix
        , bool SO2          // Storage order of the solution matrix
        , typename MT3      // Type of the right-hand side matrix
        , bool SO3 >        // Storage order of the right-hand side matrix
void solve0x0( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 0UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 0UL, "Invalid number of rows detected"    );

   MAYBE_UNUSED( A );

   resize( *X, 0UL, (*B).columns() );

   BLAZE_INTERNAL_ASSERT( isIntact( *X ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 1x1 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 1 \times 1 \f$ linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 1 \times 1 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given system matrix, \a x is the solution vector, and
// \a b is the given right-hand side vector. The function fails if the given matrix is not a
// \f$ 1 \times 1 \f$ square matrix or singular. In this case a \a std::runtime_error exception
// is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT       // Type of the system matrix
        , bool SO           // Storage order of the system matrix
        , typename VT1      // Type of the solution vector
        , bool TF1          // Transpose flag of the solution vector
        , typename VT2      // Type of the right-hand side vector
        , bool TF2 >        // Transpose flag of the right-hand side vector
void solve1x1( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 1UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 1UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 1UL, "Invalid vector size detected"       );

   CompositeType_t<MT> Atmp( *A );

   resize( *x, 1UL );
   smpAssign( *x, *b );

   if( !isDivisor( Atmp(0,0) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   (*x)[0] /= Atmp(0,0);

   BLAZE_INTERNAL_ASSERT( isIntact( *x ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 1 \times 1 \f$ linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 1 \times 1 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 1 \times 1 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT1      // Type of the system matrix
        , bool SO1          // Storage order of the system matrix
        , typename MT2      // Type of the solution matrix
        , bool SO2          // Storage order of the solution matrix
        , typename MT3      // Type of the right-hand side matrix
        , bool SO3 >        // Storage order of the right-hand side matrix
void solve1x1( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 1UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 1UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 1UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> Atmp( *A );

   if( !isDivisor( Atmp(0,0) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( *X, 1UL, (*B).columns() );
   smpAssign( *X, *B );

   const ElementType_t<MT1> invD( inv( Atmp(0,0) ) );

   for( size_t j=0UL; j<(*B).columns(); ++j ) {
      (*X)(0,j) *= invD;
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *X ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 2x2 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution vector,
// and \a b is the given right-hand side vector. The function fails if the given matrix is not a
// \f$ 2 \times 2 \f$ square matrix or singular. In this case a \a std::runtime_error exception
// is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> || ( IsHermitian_v<MT> && !IsSymmetric_v<MT> ) >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET D( A_(0,0)*A_(1,1) - A_(0,1)*A_(1,0) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   x_[0] = ( A_(1,1)*b_[0] - A_(0,1)*b_[1] ) * invD;
   x_[1] = ( A_(0,0)*b_[1] - A_(1,0)*b_[0] ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET D( A_(0,0)*A_(1,1) - A_(0,1)*A_(0,1) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   x_[0] = ( A_(1,1)*b_[0] - A_(0,1)*b_[1] ) * invD;
   x_[1] = ( A_(0,0)*b_[1] - A_(0,1)*b_[0] ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = ( b_[0] ) / A_(0,0);
   x_[1] = ( b_[1] - A_(1,0)*x_[0] ) / A_(1,1);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[0] = ( b_[0] );
   x_[1] = ( b_[1] - A_(1,0)*x_[0] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[1] = ( b_[1] ) / A_(1,1);
   x_[0] = ( b_[0] - A_(0,1)*x_[1] ) / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[1] = ( b_[1] );
   x_[0] = ( b_[0] - A_(0,1)*x_[1] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solve2x2( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 2UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[1] = b_[1] / A_(1,1);
   x_[0] = b_[0] / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< ( IsGeneral_v<MT1> || IsSymmetric_v<MT1> || IsHermitian_v<MT1> ) && !IsDiagonal_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   resize( *X, (*B).rows(), (*B).columns() );
   const ResultType_t<MT1> invA( inv( *A ) );
   smpAssign( *X, invA * (*B) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = ( B_(0,j) ) * invD0;
      X_(1,j) = ( B_(1,j) - A_(1,0)*X_(0,j) ) * invD1;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j);
      X_(1,j) = B_(1,j) - A_(1,0)*X_(0,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(1,j) = ( B_(1,j) ) * invD1;
      X_(0,j) = ( B_(0,j) - A_(0,1)*X_(1,j) ) * invD0;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(1,j) = B_(1,j);
      X_(0,j) = B_(0,j) - A_(0,1)*X_(1,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 2 \times 2 \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 2 \times 2 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 2 \times 2 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solve2x2( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 2UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 2UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j) * invD0;
      X_(1,j) = B_(1,j) * invD1;
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 3x3 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> || ( IsHermitian_v<MT> && !IsSymmetric_v<MT> ) >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1( A_(0,0)*A_(1,1) - A_(1,0)*A_(0,1) );
   const ET tmp2( A_(1,0)*A_(0,2) - A_(0,0)*A_(1,2) );
   const ET tmp3( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );

   const ET D( tmp1*A_(2,2) + tmp2*A_(2,1) + tmp3*A_(2,0) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp4( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp5( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp6( A_(0,0)*b_[1] - A_(1,0)*b_[0] );

   x_[0] = ( tmp4*A_(2,2) + tmp5*A_(2,1) + tmp3*b_[2] ) * invD;
   x_[1] = ( tmp6*A_(2,2) + tmp2*b_[2] - tmp5*A_(2,0) ) * invD;
   x_[2] = ( tmp1*b_[2] - tmp4*A_(2,0) - tmp6*A_(2,1) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1( A_(0,0)*A_(1,1) - A_(0,1)*A_(0,1) );
   const ET tmp2( A_(0,1)*A_(0,2) - A_(0,0)*A_(1,2) );
   const ET tmp3( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );

   const ET D( tmp1*A_(2,2) + tmp2*A_(1,2) + tmp3*A_(0,2) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp4( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp5( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp6( A_(0,0)*b_[1] - A_(0,1)*b_[0] );

   x_[0] = ( tmp4*A_(2,2) + tmp5*A_(1,2) + tmp3*b_[2] ) * invD;
   x_[1] = ( tmp6*A_(2,2) + tmp2*b_[2] - tmp5*A_(0,2) ) * invD;
   x_[2] = ( tmp1*b_[2] - tmp4*A_(0,2) - tmp6*A_(1,2) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = ( b_[0] ) / A_(0,0);
   x_[1] = ( b_[1] - A_(1,0)*x_[0] ) / A_(1,1);
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] ) / A_(2,2);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[0] = ( b_[0] );
   x_[1] = ( b_[1] - A_(1,0)*x_[0] );
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[2] = ( b_[2] ) / A_(2,2);
   x_[1] = ( b_[1] - A_(1,2)*x_[2] ) / A_(1,1);
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] ) / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[2] = ( b_[2] );
   x_[1] = ( b_[1] - A_(1,2)*x_[2] );
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solve3x3( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 3UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = b_[0] / A_(0,0);
   x_[1] = b_[1] / A_(1,1);
   x_[2] = b_[2] / A_(2,2);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< ( IsGeneral_v<MT1> || IsSymmetric_v<MT1> || IsHermitian_v<MT1> ) && !IsDiagonal_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   resize( *X, (*B).rows(), (*B).columns() );
   const ResultType_t<MT1> invA( inv( *A ) );
   smpAssign( *X, invA * (*B) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = ( B_(0,j) ) * invD0;
      X_(1,j) = ( B_(1,j) - A_(1,0)*X_(0,j) ) * invD1;
      X_(2,j) = ( B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j) ) * invD2;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j);
      X_(1,j) = B_(1,j) - A_(1,0)*X_(0,j);
      X_(2,j) = B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(2,j) = ( B_(2,j) ) * invD2;
      X_(1,j) = ( B_(1,j) - A_(1,2)*X_(2,j) ) * invD1;
      X_(0,j) = ( B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) ) * invD0;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(2,j) = B_(2,j);
      X_(1,j) = B_(1,j) - A_(1,2)*X_(2,j);
      X_(0,j) = B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 3 \times 3 \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 3 \times 3 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 3 \times 3 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solve3x3( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 3UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 3UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j) * invD0;
      X_(1,j) = B_(1,j) * invD1;
      X_(2,j) = B_(2,j) * invD2;
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 4x4 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> || ( IsHermitian_v<MT> && !IsSymmetric_v<MT> ) >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(1,0)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(1,0)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(1,0)*A_(0,3) );
   const ET tmp4 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp5 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp6 ( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );

   const ET tmp7 ( tmp1*A_(2,2) - tmp2*A_(2,1) + tmp4*A_(2,0) );
   const ET tmp8 ( tmp1*A_(2,3) - tmp3*A_(2,1) + tmp5*A_(2,0) );
   const ET tmp9 ( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp6*A_(2,0) );
   const ET tmp10( tmp4*A_(2,3) - tmp5*A_(2,2) + tmp6*A_(2,1) );

   const ET D( tmp7*A_(3,3) - tmp8*A_(3,2) + tmp9*A_(3,1) - tmp10*A_(3,0) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp11( A_(0,0)*b_[1] - A_(1,0)*b_[0] );
   const ET tmp12( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp13( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp14( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp15( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp16( A_(1,3)*b_[0] - A_(0,3)*b_[1] );

   const ET tmp17( tmp14*A_(2,2) - tmp15*A_(2,1) + tmp4*b_[2] );
   const ET tmp18( tmp14*A_(2,3) - tmp16*A_(2,1) + tmp5*b_[2] );
   const ET tmp19( tmp15*A_(2,3) - tmp16*A_(2,2) + tmp6*b_[2] );

   const ET tmp20( tmp11*A_(2,2) - tmp2*b_[2] + tmp15*A_(2,0) );
   const ET tmp21( tmp11*A_(2,3) - tmp3*b_[2] + tmp16*A_(2,0) );
   const ET tmp22( tmp12*A_(2,3) - tmp5*b_[2] + tmp16*A_(2,1) );

   const ET tmp23( tmp1*b_[2] - tmp11*A_(2,1) + tmp12*A_(2,0) );
   const ET tmp24( tmp2*b_[2] - tmp11*A_(2,2) + tmp13*A_(2,0) );
   const ET tmp25( tmp4*b_[2] - tmp12*A_(2,2) + tmp13*A_(2,1) );

   x_[0] = ( tmp17*A_(3,3) - tmp18*A_(3,2) + tmp19*A_(3,1) - tmp10*b_[3] ) * invD;
   x_[1] = ( tmp20*A_(3,3) - tmp21*A_(3,2) + tmp9*b_[3] - tmp19*A_(3,0) ) * invD;
   x_[2] = ( tmp23*A_(3,3) - tmp8*b_[3] + tmp21*A_(3,1) - tmp22*A_(3,0) ) * invD;
   x_[3] = ( tmp7*b_[3] - tmp23*A_(3,2) + tmp24*A_(3,1) - tmp25*A_(3,0) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(0,1)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(0,1)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(0,1)*A_(0,3) );
   const ET tmp4 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp5 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp6 ( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );

   const ET tmp7 ( tmp1*A_(2,2) - tmp2*A_(1,2) + tmp4*A_(0,2) );
   const ET tmp8 ( tmp1*A_(2,3) - tmp3*A_(1,2) + tmp5*A_(0,2) );
   const ET tmp9 ( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp6*A_(0,2) );
   const ET tmp10( tmp4*A_(2,3) - tmp5*A_(2,2) + tmp6*A_(1,2) );

   const ET D( tmp7*A_(3,3) - tmp8*A_(2,3) + tmp9*A_(1,3) - tmp10*A_(0,3) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp11( A_(0,0)*b_[1] - A_(0,1)*b_[0] );
   const ET tmp12( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp13( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp14( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp15( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp16( A_(1,3)*b_[0] - A_(0,3)*b_[1] );

   const ET tmp17( tmp14*A_(2,2) - tmp15*A_(1,2) + tmp4*b_[2] );
   const ET tmp18( tmp14*A_(2,3) - tmp16*A_(1,2) + tmp5*b_[2] );
   const ET tmp19( tmp15*A_(2,3) - tmp16*A_(2,2) + tmp6*b_[2] );

   const ET tmp20( tmp11*A_(2,2) - tmp2*b_[2] + tmp15*A_(0,2) );
   const ET tmp21( tmp11*A_(2,3) - tmp3*b_[2] + tmp16*A_(0,2) );
   const ET tmp22( tmp12*A_(2,3) - tmp5*b_[2] + tmp16*A_(1,2) );

   const ET tmp23( tmp1*b_[2] - tmp11*A_(1,2) + tmp12*A_(0,2) );
   const ET tmp24( tmp2*b_[2] - tmp11*A_(2,2) + tmp13*A_(0,2) );
   const ET tmp25( tmp4*b_[2] - tmp12*A_(2,2) + tmp13*A_(1,2) );

   x_[0] = ( tmp17*A_(3,3) - tmp18*A_(2,3) + tmp19*A_(1,3) - tmp10*b_[3] ) * invD;
   x_[1] = ( tmp20*A_(3,3) - tmp21*A_(2,3) + tmp9*b_[3] - tmp19*A_(0,3) ) * invD;
   x_[2] = ( tmp23*A_(3,3) - tmp8*b_[3] + tmp21*A_(1,3) - tmp22*A_(0,3) ) * invD;
   x_[3] = ( tmp7*b_[3] - tmp23*A_(2,3) + tmp24*A_(1,3) - tmp25*A_(0,3) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = ( b_[0] ) / A_(0,0);
   x_[1] = ( b_[1] - A_(1,0)*x_[0] ) / A_(1,1);
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] ) / A_(2,2);
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] ) / A_(3,3);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[0] = ( b_[0] );
   x_[1] = ( b_[1] - A_(1,0)*x_[0] );
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] );
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[3] = ( b_[3] ) / A_(3,3);
   x_[2] = ( b_[2] - A_(2,3)*x_[3] ) / A_(2,2);
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] ) / A_(1,1);
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] ) / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[3] = ( b_[3] );
   x_[2] = ( b_[2] - A_(2,3)*x_[3] );
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] );
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solve4x4( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 4UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = b_[0] / A_(0,0);
   x_[1] = b_[1] / A_(1,1);
   x_[2] = b_[2] / A_(2,2);
   x_[3] = b_[3] / A_(3,3);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< ( IsGeneral_v<MT1> || IsSymmetric_v<MT1> || IsHermitian_v<MT1> ) && !IsDiagonal_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   resize( *X, (*B).rows(), (*B).columns() );
   const ResultType_t<MT1> invA( inv( *A ) );
   smpAssign( *X, invA * (*B) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = ( B_(0,j) ) * invD0;
      X_(1,j) = ( B_(1,j) - A_(1,0)*X_(0,j) ) * invD1;
      X_(2,j) = ( B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j) ) * invD2;
      X_(3,j) = ( B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j) ) * invD3;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j);
      X_(1,j) = B_(1,j) - A_(1,0)*X_(0,j);
      X_(2,j) = B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j);
      X_(3,j) = B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(3,j) = ( B_(3,j) ) * invD3;
      X_(2,j) = ( B_(2,j) - A_(2,3)*X_(3,j) ) * invD2;
      X_(1,j) = ( B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j) ) * invD1;
      X_(0,j) = ( B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j) ) * invD0;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(3,j) = B_(3,j);
      X_(2,j) = B_(2,j) - A_(2,3)*X_(3,j);
      X_(1,j) = B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j);
      X_(0,j) = B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 4 \times 4 \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 4 \times 4 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 4 \times 4 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solve4x4( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 4UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 4UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j) * invD0;
      X_(1,j) = B_(1,j) * invD1;
      X_(2,j) = B_(2,j) * invD2;
      X_(3,j) = B_(3,j) * invD3;
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 5x5 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> || ( IsHermitian_v<MT> && !IsSymmetric_v<MT> ) >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(1,0)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(1,0)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(1,0)*A_(0,3) );
   const ET tmp4 ( A_(0,0)*A_(1,4) - A_(1,0)*A_(0,4) );
   const ET tmp5 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp6 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp7 ( A_(0,1)*A_(1,4) - A_(1,1)*A_(0,4) );
   const ET tmp8 ( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );
   const ET tmp9 ( A_(0,2)*A_(1,4) - A_(1,2)*A_(0,4) );
   const ET tmp10( A_(0,3)*A_(1,4) - A_(1,3)*A_(0,4) );

   const ET tmp11( A_(0,0)*b_[1] - A_(1,0)*b_[0] );
   const ET tmp12( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp13( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp14( A_(0,3)*b_[1] - A_(1,3)*b_[0] );

   const ET tmp15( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp16( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp17( A_(1,3)*b_[0] - A_(0,3)*b_[1] );
   const ET tmp18( A_(1,4)*b_[0] - A_(0,4)*b_[1] );

   const ET tmp19( tmp1*A_(2,2) - tmp2*A_(2,1) + tmp5*A_(2,0) );
   const ET tmp20( tmp1*A_(2,3) - tmp3*A_(2,1) + tmp6*A_(2,0) );
   const ET tmp21( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp8*A_(2,0) );
   const ET tmp22( tmp5*A_(2,3) - tmp6*A_(2,2) + tmp8*A_(2,1) );
   const ET tmp23( tmp1*A_(2,4) - tmp4*A_(2,1) + tmp7*A_(2,0) );
   const ET tmp24( tmp2*A_(2,4) - tmp4*A_(2,2) + tmp9*A_(2,0) );
   const ET tmp25( tmp5*A_(2,4) - tmp7*A_(2,2) + tmp9*A_(2,1) );
   const ET tmp26( tmp3*A_(2,4) - tmp4*A_(2,3) + tmp10*A_(2,0) );
   const ET tmp27( tmp6*A_(2,4) - tmp7*A_(2,3) + tmp10*A_(2,1) );
   const ET tmp28( tmp8*A_(2,4) - tmp9*A_(2,3) + tmp10*A_(2,2) );

   const ET tmp29( tmp19*A_(3,3) - tmp20*A_(3,2) + tmp21*A_(3,1) - tmp22*A_(3,0) );
   const ET tmp30( tmp19*A_(3,4) - tmp23*A_(3,2) + tmp24*A_(3,1) - tmp25*A_(3,0) );
   const ET tmp31( tmp20*A_(3,4) - tmp23*A_(3,3) + tmp26*A_(3,1) - tmp27*A_(3,0) );
   const ET tmp32( tmp21*A_(3,4) - tmp24*A_(3,3) + tmp26*A_(3,2) - tmp28*A_(3,0) );
   const ET tmp33( tmp22*A_(3,4) - tmp25*A_(3,3) + tmp27*A_(3,2) - tmp28*A_(3,1) );

   const ET D( tmp29*A_(4,4) - tmp30*A_(4,3) + tmp31*A_(4,2) - tmp32*A_(4,1) + tmp33*A_(4,0) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp34( tmp15*A_(2,2) - tmp16*A_(2,1) + tmp5*b_[2] );
   const ET tmp35( tmp15*A_(2,3) - tmp17*A_(2,1) + tmp6*b_[2] );
   const ET tmp36( tmp16*A_(2,3) - tmp17*A_(2,2) + tmp8*b_[2] );
   const ET tmp37( tmp15*A_(2,4) - tmp18*A_(2,1) + tmp7*b_[2] );
   const ET tmp38( tmp16*A_(2,4) - tmp18*A_(2,2) + tmp9*b_[2] );
   const ET tmp39( tmp17*A_(2,4) - tmp18*A_(2,3) + tmp10*b_[2] );

   const ET tmp40( tmp11*A_(2,2) - tmp2*b_[2] + tmp16*A_(2,0) );
   const ET tmp41( tmp11*A_(2,3) - tmp3*b_[2] + tmp17*A_(2,0) );
   const ET tmp42( tmp12*A_(2,3) - tmp6*b_[2] + tmp17*A_(2,1) );
   const ET tmp43( tmp11*A_(2,4) - tmp4*b_[2] + tmp18*A_(2,0) );
   const ET tmp44( tmp12*A_(2,4) - tmp7*b_[2] + tmp18*A_(2,1) );
   const ET tmp45( tmp13*A_(2,4) - tmp9*b_[2] + tmp18*A_(2,2) );

   const ET tmp46( tmp1*b_[2] - tmp11*A_(2,1) + tmp12*A_(2,0) );
   const ET tmp47( tmp2*b_[2] - tmp11*A_(2,2) + tmp13*A_(2,0) );
   const ET tmp48( tmp5*b_[2] - tmp12*A_(2,2) + tmp13*A_(2,1) );
   const ET tmp49( tmp3*b_[2] - tmp11*A_(2,3) + tmp14*A_(2,0) );
   const ET tmp50( tmp6*b_[2] - tmp12*A_(2,3) + tmp14*A_(2,1) );
   const ET tmp51( tmp8*b_[2] - tmp13*A_(2,3) + tmp14*A_(2,2) );

   const ET tmp52( tmp34*A_(3,3) - tmp35*A_(3,2) + tmp36*A_(3,1) - tmp22*b_[3] );
   const ET tmp53( tmp34*A_(3,4) - tmp37*A_(3,2) + tmp38*A_(3,1) - tmp25*b_[3] );
   const ET tmp54( tmp35*A_(3,4) - tmp37*A_(3,3) + tmp39*A_(3,1) - tmp27*b_[3] );
   const ET tmp55( tmp36*A_(3,4) - tmp38*A_(3,3) + tmp39*A_(3,2) - tmp28*b_[3] );

   const ET tmp56( tmp40*A_(3,3) - tmp41*A_(3,2) + tmp21*b_[3] - tmp36*A_(3,0) );
   const ET tmp57( tmp40*A_(3,4) - tmp43*A_(3,2) + tmp24*b_[3] - tmp38*A_(3,0) );
   const ET tmp58( tmp41*A_(3,4) - tmp43*A_(3,3) + tmp26*b_[3] - tmp39*A_(3,0) );
   const ET tmp59( tmp42*A_(3,4) - tmp44*A_(3,3) + tmp27*b_[3] - tmp39*A_(3,1) );

   const ET tmp60( tmp46*A_(3,3) - tmp20*b_[3] + tmp41*A_(3,1) - tmp42*A_(3,0) );
   const ET tmp61( tmp46*A_(3,4) - tmp23*b_[3] + tmp43*A_(3,1) - tmp44*A_(3,0) );
   const ET tmp62( tmp47*A_(3,4) - tmp24*b_[3] + tmp43*A_(3,2) - tmp45*A_(3,0) );
   const ET tmp63( tmp48*A_(3,4) - tmp25*b_[3] + tmp44*A_(3,2) - tmp45*A_(3,1) );

   const ET tmp64( tmp19*b_[3] - tmp46*A_(3,2) + tmp47*A_(3,1) - tmp48*A_(3,0) );
   const ET tmp65( tmp20*b_[3] - tmp46*A_(3,3) + tmp49*A_(3,1) - tmp50*A_(3,0) );
   const ET tmp66( tmp21*b_[3] - tmp47*A_(3,3) + tmp49*A_(3,2) - tmp51*A_(3,0) );
   const ET tmp67( tmp22*b_[3] - tmp48*A_(3,3) + tmp50*A_(3,2) - tmp51*A_(3,1) );

   x_[0] = ( tmp52*A_(4,4) - tmp53*A_(4,3) + tmp54*A_(4,2) - tmp55*A_(4,1) + tmp33*b_[4] ) * invD;
   x_[1] = ( tmp56*A_(4,4) - tmp57*A_(4,3) + tmp58*A_(4,2) - tmp32*b_[4] + tmp55*A_(4,0) ) * invD;
   x_[2] = ( tmp60*A_(4,4) - tmp61*A_(4,3) + tmp31*b_[4] - tmp58*A_(4,1) + tmp59*A_(4,0) ) * invD;
   x_[3] = ( tmp64*A_(4,4) - tmp30*b_[4] + tmp61*A_(4,2) - tmp62*A_(4,1) + tmp63*A_(4,0) ) * invD;
   x_[4] = ( tmp29*b_[4] - tmp64*A_(4,3) + tmp65*A_(4,2) - tmp66*A_(4,1) + tmp67*A_(4,0) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(0,1)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(0,1)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(0,1)*A_(0,3) );
   const ET tmp4 ( A_(0,0)*A_(1,4) - A_(0,1)*A_(0,4) );
   const ET tmp5 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp6 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp7 ( A_(0,1)*A_(1,4) - A_(1,1)*A_(0,4) );
   const ET tmp8 ( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );
   const ET tmp9 ( A_(0,2)*A_(1,4) - A_(1,2)*A_(0,4) );
   const ET tmp10( A_(0,3)*A_(1,4) - A_(1,3)*A_(0,4) );

   const ET tmp11( A_(0,0)*b_[1] - A_(0,1)*b_[0] );
   const ET tmp12( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp13( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp14( A_(0,3)*b_[1] - A_(1,3)*b_[0] );

   const ET tmp15( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp16( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp17( A_(1,3)*b_[0] - A_(0,3)*b_[1] );
   const ET tmp18( A_(1,4)*b_[0] - A_(0,4)*b_[1] );

   const ET tmp19( tmp1*A_(2,2) - tmp2*A_(1,2) + tmp5*A_(0,2) );
   const ET tmp20( tmp1*A_(2,3) - tmp3*A_(1,2) + tmp6*A_(0,2) );
   const ET tmp21( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp8*A_(0,2) );
   const ET tmp22( tmp5*A_(2,3) - tmp6*A_(2,2) + tmp8*A_(1,2) );
   const ET tmp23( tmp1*A_(2,4) - tmp4*A_(1,2) + tmp7*A_(0,2) );
   const ET tmp24( tmp2*A_(2,4) - tmp4*A_(2,2) + tmp9*A_(0,2) );
   const ET tmp25( tmp5*A_(2,4) - tmp7*A_(2,2) + tmp9*A_(1,2) );
   const ET tmp26( tmp3*A_(2,4) - tmp4*A_(2,3) + tmp10*A_(0,2) );
   const ET tmp27( tmp6*A_(2,4) - tmp7*A_(2,3) + tmp10*A_(1,2) );
   const ET tmp28( tmp8*A_(2,4) - tmp9*A_(2,3) + tmp10*A_(2,2) );

   const ET tmp29( tmp19*A_(3,3) - tmp20*A_(2,3) + tmp21*A_(1,3) - tmp22*A_(0,3) );
   const ET tmp30( tmp19*A_(3,4) - tmp23*A_(2,3) + tmp24*A_(1,3) - tmp25*A_(0,3) );
   const ET tmp31( tmp20*A_(3,4) - tmp23*A_(3,3) + tmp26*A_(1,3) - tmp27*A_(0,3) );
   const ET tmp32( tmp21*A_(3,4) - tmp24*A_(3,3) + tmp26*A_(2,3) - tmp28*A_(0,3) );
   const ET tmp33( tmp22*A_(3,4) - tmp25*A_(3,3) + tmp27*A_(2,3) - tmp28*A_(1,3) );

   const ET D( tmp29*A_(4,4) - tmp30*A_(3,4) + tmp31*A_(2,4) - tmp32*A_(1,4) + tmp33*A_(0,4) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp34( tmp15*A_(2,2) - tmp16*A_(1,2) + tmp5*b_[2] );
   const ET tmp35( tmp15*A_(2,3) - tmp17*A_(1,2) + tmp6*b_[2] );
   const ET tmp36( tmp16*A_(2,3) - tmp17*A_(2,2) + tmp8*b_[2] );
   const ET tmp37( tmp15*A_(2,4) - tmp18*A_(1,2) + tmp7*b_[2] );
   const ET tmp38( tmp16*A_(2,4) - tmp18*A_(2,2) + tmp9*b_[2] );
   const ET tmp39( tmp17*A_(2,4) - tmp18*A_(2,3) + tmp10*b_[2] );

   const ET tmp40( tmp11*A_(2,2) - tmp2*b_[2] + tmp16*A_(0,2) );
   const ET tmp41( tmp11*A_(2,3) - tmp3*b_[2] + tmp17*A_(0,2) );
   const ET tmp42( tmp12*A_(2,3) - tmp6*b_[2] + tmp17*A_(1,2) );
   const ET tmp43( tmp11*A_(2,4) - tmp4*b_[2] + tmp18*A_(0,2) );
   const ET tmp44( tmp12*A_(2,4) - tmp7*b_[2] + tmp18*A_(1,2) );
   const ET tmp45( tmp13*A_(2,4) - tmp9*b_[2] + tmp18*A_(2,2) );

   const ET tmp46( tmp1*b_[2] - tmp11*A_(1,2) + tmp12*A_(0,2) );
   const ET tmp47( tmp2*b_[2] - tmp11*A_(2,2) + tmp13*A_(0,2) );
   const ET tmp48( tmp5*b_[2] - tmp12*A_(2,2) + tmp13*A_(1,2) );
   const ET tmp49( tmp3*b_[2] - tmp11*A_(2,3) + tmp14*A_(0,2) );
   const ET tmp50( tmp6*b_[2] - tmp12*A_(2,3) + tmp14*A_(1,2) );
   const ET tmp51( tmp8*b_[2] - tmp13*A_(2,3) + tmp14*A_(2,2) );

   const ET tmp52( tmp34*A_(3,3) - tmp35*A_(2,3) + tmp36*A_(1,3) - tmp22*b_[3] );
   const ET tmp53( tmp34*A_(3,4) - tmp37*A_(2,3) + tmp38*A_(1,3) - tmp25*b_[3] );
   const ET tmp54( tmp35*A_(3,4) - tmp37*A_(3,3) + tmp39*A_(1,3) - tmp27*b_[3] );
   const ET tmp55( tmp36*A_(3,4) - tmp38*A_(3,3) + tmp39*A_(2,3) - tmp28*b_[3] );

   const ET tmp56( tmp40*A_(3,3) - tmp41*A_(2,3) + tmp21*b_[3] - tmp36*A_(0,3) );
   const ET tmp57( tmp40*A_(3,4) - tmp43*A_(2,3) + tmp24*b_[3] - tmp38*A_(0,3) );
   const ET tmp58( tmp41*A_(3,4) - tmp43*A_(3,3) + tmp26*b_[3] - tmp39*A_(0,3) );
   const ET tmp59( tmp42*A_(3,4) - tmp44*A_(3,3) + tmp27*b_[3] - tmp39*A_(1,3) );

   const ET tmp60( tmp46*A_(3,3) - tmp20*b_[3] + tmp41*A_(1,3) - tmp42*A_(0,3) );
   const ET tmp61( tmp46*A_(3,4) - tmp23*b_[3] + tmp43*A_(1,3) - tmp44*A_(0,3) );
   const ET tmp62( tmp47*A_(3,4) - tmp24*b_[3] + tmp43*A_(2,3) - tmp45*A_(0,3) );
   const ET tmp63( tmp48*A_(3,4) - tmp25*b_[3] + tmp44*A_(2,3) - tmp45*A_(1,3) );

   const ET tmp64( tmp19*b_[3] - tmp46*A_(2,3) + tmp47*A_(1,3) - tmp48*A_(0,3) );
   const ET tmp65( tmp20*b_[3] - tmp46*A_(3,3) + tmp49*A_(1,3) - tmp50*A_(0,3) );
   const ET tmp66( tmp21*b_[3] - tmp47*A_(3,3) + tmp49*A_(2,3) - tmp51*A_(0,3) );
   const ET tmp67( tmp22*b_[3] - tmp48*A_(3,3) + tmp50*A_(2,3) - tmp51*A_(1,3) );

   x_[0] = ( tmp52*A_(4,4) - tmp53*A_(3,4) + tmp54*A_(2,4) - tmp55*A_(1,4) + tmp33*b_[4] ) * invD;
   x_[1] = ( tmp56*A_(4,4) - tmp57*A_(3,4) + tmp58*A_(2,4) - tmp32*b_[4] + tmp55*A_(0,4) ) * invD;
   x_[2] = ( tmp60*A_(4,4) - tmp61*A_(3,4) + tmp31*b_[4] - tmp58*A_(1,4) + tmp59*A_(0,4) ) * invD;
   x_[3] = ( tmp64*A_(4,4) - tmp30*b_[4] + tmp61*A_(2,4) - tmp62*A_(1,4) + tmp63*A_(0,4) ) * invD;
   x_[4] = ( tmp29*b_[4] - tmp64*A_(3,4) + tmp65*A_(2,4) - tmp66*A_(1,4) + tmp67*A_(0,4) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = ( b_[0] ) / A_(0,0);
   x_[1] = ( b_[1] - A_(1,0)*x_[0] ) / A_(1,1);
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] ) / A_(2,2);
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] ) / A_(3,3);
   x_[4] = ( b_[4] - A_(4,0)*x_[0] - A_(4,1)*x_[1] - A_(4,2)*x_[2] - A_(4,3)*x_[3] ) / A_(4,4);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[0] = ( b_[0] );
   x_[1] = ( b_[1] - A_(1,0)*x_[0] );
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] );
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] );
   x_[4] = ( b_[4] - A_(4,0)*x_[0] - A_(4,1)*x_[1] - A_(4,2)*x_[2] - A_(4,3)*x_[3] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[4] = ( b_[4] ) / A_(4,4);
   x_[3] = ( b_[3] - A_(3,4)*x_[4] ) / A_(3,3);
   x_[2] = ( b_[2] - A_(2,3)*x_[3] - A_(2,4)*x_[4] ) / A_(2,2);
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] - A_(1,4)*x_[4] ) / A_(1,1);
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] - A_(0,4)*x_[4] ) / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[4] = ( b_[4] );
   x_[3] = ( b_[3] - A_(3,4)*x_[4] );
   x_[2] = ( b_[2] - A_(2,3)*x_[3] - A_(2,4)*x_[4] );
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] - A_(1,4)*x_[4] );
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] - A_(0,4)*x_[4] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solve5x5( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 5UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = b_[0] / A_(0,0);
   x_[1] = b_[1] / A_(1,1);
   x_[2] = b_[2] / A_(2,2);
   x_[3] = b_[3] / A_(3,3);
   x_[4] = b_[4] / A_(4,4);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< ( IsGeneral_v<MT1> || IsSymmetric_v<MT1> || IsHermitian_v<MT1> ) && !IsDiagonal_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   resize( *X, (*B).rows(), (*B).columns() );
   const ResultType_t<MT1> invA( inv( *A ) );
   smpAssign( *X, invA * (*B) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = ( B_(0,j) ) * invD0;
      X_(1,j) = ( B_(1,j) - A_(1,0)*X_(0,j) ) * invD1;
      X_(2,j) = ( B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j) ) * invD2;
      X_(3,j) = ( B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j) ) * invD3;
      X_(4,j) = ( B_(4,j) - A_(4,0)*X_(0,j) - A_(4,1)*X_(1,j) - A_(4,2)*X_(2,j) - A_(4,3)*X_(3,j) ) * invD4;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j);
      X_(1,j) = B_(1,j) - A_(1,0)*X_(0,j);
      X_(2,j) = B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j);
      X_(3,j) = B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j);
      X_(4,j) = B_(4,j) - A_(4,0)*X_(0,j) - A_(4,1)*X_(1,j) - A_(4,2)*X_(2,j) - A_(4,3)*X_(3,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(4,j) = ( B_(4,j) ) * invD4;
      X_(3,j) = ( B_(3,j) - A_(3,4)*X_(4,j) ) * invD3;
      X_(2,j) = ( B_(2,j) - A_(2,3)*X_(3,j) - A_(2,4)*X_(4,j) ) * invD2;
      X_(1,j) = ( B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j) - A_(1,4)*X_(4,j) ) * invD1;
      X_(0,j) = ( B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j) - A_(0,4)*X_(4,j) ) * invD0;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(4,j) = B_(4,j);
      X_(3,j) = B_(3,j) - A_(3,4)*X_(4,j);
      X_(2,j) = B_(2,j) - A_(2,3)*X_(3,j) - A_(2,4)*X_(4,j);
      X_(1,j) = B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j) - A_(1,4)*X_(4,j);
      X_(0,j) = B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j) - A_(0,4)*X_(4,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 5 \times 5 \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 5 \times 5 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 5 \times 5 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solve5x5( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 5UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 5UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j) * invD0;
      X_(1,j) = B_(1,j) * invD1;
      X_(2,j) = B_(2,j) * invD2;
      X_(3,j) = B_(3,j) * invD3;
      X_(4,j) = B_(4,j) * invD4;
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING 6x6 LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> || ( IsHermitian_v<MT> && !IsSymmetric_v<MT> ) >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(1,0)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(1,0)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(1,0)*A_(0,3) );
   const ET tmp4 ( A_(0,0)*A_(1,4) - A_(1,0)*A_(0,4) );
   const ET tmp5 ( A_(0,0)*A_(1,5) - A_(1,0)*A_(0,5) );
   const ET tmp6 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp7 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp8 ( A_(0,1)*A_(1,4) - A_(1,1)*A_(0,4) );
   const ET tmp9 ( A_(0,1)*A_(1,5) - A_(1,1)*A_(0,5) );
   const ET tmp10( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );
   const ET tmp11( A_(0,2)*A_(1,4) - A_(1,2)*A_(0,4) );
   const ET tmp12( A_(0,2)*A_(1,5) - A_(1,2)*A_(0,5) );
   const ET tmp13( A_(0,3)*A_(1,4) - A_(1,3)*A_(0,4) );
   const ET tmp14( A_(0,3)*A_(1,5) - A_(1,3)*A_(0,5) );
   const ET tmp15( A_(0,4)*A_(1,5) - A_(1,4)*A_(0,5) );

   const ET tmp16( tmp1*A_(2,2) - tmp2*A_(2,1) + tmp6*A_(2,0) );
   const ET tmp17( tmp1*A_(2,3) - tmp3*A_(2,1) + tmp7*A_(2,0) );
   const ET tmp18( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp10*A_(2,0) );
   const ET tmp19( tmp6*A_(2,3) - tmp7*A_(2,2) + tmp10*A_(2,1) );
   const ET tmp20( tmp1*A_(2,4) - tmp4*A_(2,1) + tmp8*A_(2,0) );
   const ET tmp21( tmp2*A_(2,4) - tmp4*A_(2,2) + tmp11*A_(2,0) );
   const ET tmp22( tmp6*A_(2,4) - tmp8*A_(2,2) + tmp11*A_(2,1) );
   const ET tmp23( tmp3*A_(2,4) - tmp4*A_(2,3) + tmp13*A_(2,0) );
   const ET tmp24( tmp7*A_(2,4) - tmp8*A_(2,3) + tmp13*A_(2,1) );
   const ET tmp25( tmp10*A_(2,4) - tmp11*A_(2,3) + tmp13*A_(2,2) );
   const ET tmp26( tmp1*A_(2,5) - tmp5*A_(2,1) + tmp9*A_(2,0) );
   const ET tmp27( tmp2*A_(2,5) - tmp5*A_(2,2) + tmp12*A_(2,0) );
   const ET tmp28( tmp6*A_(2,5) - tmp9*A_(2,2) + tmp12*A_(2,1) );
   const ET tmp29( tmp3*A_(2,5) - tmp5*A_(2,3) + tmp14*A_(2,0) );
   const ET tmp30( tmp7*A_(2,5) - tmp9*A_(2,3) + tmp14*A_(2,1) );
   const ET tmp31( tmp10*A_(2,5) - tmp12*A_(2,3) + tmp14*A_(2,2) );
   const ET tmp32( tmp4*A_(2,5) - tmp5*A_(2,4) + tmp15*A_(2,0) );
   const ET tmp33( tmp8*A_(2,5) - tmp9*A_(2,4) + tmp15*A_(2,1) );
   const ET tmp34( tmp11*A_(2,5) - tmp12*A_(2,4) + tmp15*A_(2,2) );
   const ET tmp35( tmp13*A_(2,5) - tmp14*A_(2,4) + tmp15*A_(2,3) );

   const ET tmp36( tmp16*A_(3,3) - tmp17*A_(3,2) + tmp18*A_(3,1) - tmp19*A_(3,0) );
   const ET tmp37( tmp16*A_(3,4) - tmp20*A_(3,2) + tmp21*A_(3,1) - tmp22*A_(3,0) );
   const ET tmp38( tmp17*A_(3,4) - tmp20*A_(3,3) + tmp23*A_(3,1) - tmp24*A_(3,0) );
   const ET tmp39( tmp18*A_(3,4) - tmp21*A_(3,3) + tmp23*A_(3,2) - tmp25*A_(3,0) );
   const ET tmp40( tmp19*A_(3,4) - tmp22*A_(3,3) + tmp24*A_(3,2) - tmp25*A_(3,1) );
   const ET tmp41( tmp16*A_(3,5) - tmp26*A_(3,2) + tmp27*A_(3,1) - tmp28*A_(3,0) );
   const ET tmp42( tmp17*A_(3,5) - tmp26*A_(3,3) + tmp29*A_(3,1) - tmp30*A_(3,0) );
   const ET tmp43( tmp18*A_(3,5) - tmp27*A_(3,3) + tmp29*A_(3,2) - tmp31*A_(3,0) );
   const ET tmp44( tmp19*A_(3,5) - tmp28*A_(3,3) + tmp30*A_(3,2) - tmp31*A_(3,1) );
   const ET tmp45( tmp20*A_(3,5) - tmp26*A_(3,4) + tmp32*A_(3,1) - tmp33*A_(3,0) );
   const ET tmp46( tmp21*A_(3,5) - tmp27*A_(3,4) + tmp32*A_(3,2) - tmp34*A_(3,0) );
   const ET tmp47( tmp22*A_(3,5) - tmp28*A_(3,4) + tmp33*A_(3,2) - tmp34*A_(3,1) );
   const ET tmp48( tmp23*A_(3,5) - tmp29*A_(3,4) + tmp32*A_(3,3) - tmp35*A_(3,0) );
   const ET tmp49( tmp24*A_(3,5) - tmp30*A_(3,4) + tmp33*A_(3,3) - tmp35*A_(3,1) );
   const ET tmp50( tmp25*A_(3,5) - tmp31*A_(3,4) + tmp34*A_(3,3) - tmp35*A_(3,2) );

   const ET tmp51( tmp36*A_(4,4) - tmp37*A_(4,3) + tmp38*A_(4,2) - tmp39*A_(4,1) + tmp40*A_(4,0) );
   const ET tmp52( tmp36*A_(4,5) - tmp41*A_(4,3) + tmp42*A_(4,2) - tmp43*A_(4,1) + tmp44*A_(4,0) );
   const ET tmp53( tmp37*A_(4,5) - tmp41*A_(4,4) + tmp45*A_(4,2) - tmp46*A_(4,1) + tmp47*A_(4,0) );
   const ET tmp54( tmp38*A_(4,5) - tmp42*A_(4,4) + tmp45*A_(4,3) - tmp48*A_(4,1) + tmp49*A_(4,0) );
   const ET tmp55( tmp39*A_(4,5) - tmp43*A_(4,4) + tmp46*A_(4,3) - tmp48*A_(4,2) + tmp50*A_(4,0) );
   const ET tmp56( tmp40*A_(4,5) - tmp44*A_(4,4) + tmp47*A_(4,3) - tmp49*A_(4,2) + tmp50*A_(4,1) );

   const ET D( tmp51*A_(5,5) - tmp52*A_(5,4) + tmp53*A_(5,3) - tmp54*A_(5,2) + tmp55*A_(5,1) - tmp56*A_(5,0) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp57( A_(0,0)*b_[1] - A_(1,0)*b_[0] );
   const ET tmp58( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp59( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp60( A_(0,3)*b_[1] - A_(1,3)*b_[0] );
   const ET tmp61( A_(0,4)*b_[1] - A_(1,4)*b_[0] );

   const ET tmp62( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp63( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp64( A_(1,3)*b_[0] - A_(0,3)*b_[1] );
   const ET tmp65( A_(1,4)*b_[0] - A_(0,4)*b_[1] );
   const ET tmp66( A_(1,5)*b_[0] - A_(0,5)*b_[1] );

   const ET tmp67( tmp62*A_(2,2) - tmp63*A_(2,1) + tmp6*b_[2] );
   const ET tmp68( tmp62*A_(2,3) - tmp64*A_(2,1) + tmp7*b_[2] );
   const ET tmp69( tmp63*A_(2,3) - tmp64*A_(2,2) + tmp10*b_[2] );
   const ET tmp70( tmp62*A_(2,4) - tmp65*A_(2,1) + tmp8*b_[2] );
   const ET tmp71( tmp63*A_(2,4) - tmp65*A_(2,2) + tmp11*b_[2] );
   const ET tmp72( tmp64*A_(2,4) - tmp65*A_(2,3) + tmp13*b_[2] );
   const ET tmp73( tmp62*A_(2,5) - tmp66*A_(2,1) + tmp9*b_[2] );
   const ET tmp74( tmp63*A_(2,5) - tmp66*A_(2,2) + tmp12*b_[2] );
   const ET tmp75( tmp64*A_(2,5) - tmp66*A_(2,3) + tmp14*b_[2] );
   const ET tmp76( tmp65*A_(2,5) - tmp66*A_(2,4) + tmp15*b_[2] );

   const ET tmp77( tmp57*A_(2,2) - tmp2*b_[2] + tmp63*A_(2,0) );
   const ET tmp78( tmp57*A_(2,3) - tmp3*b_[2] + tmp64*A_(2,0) );
   const ET tmp79( tmp58*A_(2,3) - tmp7*b_[2] + tmp64*A_(2,1) );
   const ET tmp80( tmp57*A_(2,4) - tmp4*b_[2] + tmp65*A_(2,0) );
   const ET tmp81( tmp58*A_(2,4) - tmp8*b_[2] + tmp65*A_(2,1) );
   const ET tmp82( tmp59*A_(2,4) - tmp11*b_[2] + tmp65*A_(2,2) );
   const ET tmp83( tmp57*A_(2,5) - tmp5*b_[2] + tmp66*A_(2,0) );
   const ET tmp84( tmp58*A_(2,5) - tmp9*b_[2] + tmp66*A_(2,1) );
   const ET tmp85( tmp59*A_(2,5) - tmp12*b_[2] + tmp66*A_(2,2) );
   const ET tmp86( tmp60*A_(2,5) - tmp14*b_[2] + tmp66*A_(2,3) );

   const ET tmp87( tmp1*b_[2] - tmp57*A_(2,1) + tmp58*A_(2,0) );
   const ET tmp88( tmp2*b_[2] - tmp57*A_(2,2) + tmp59*A_(2,0) );
   const ET tmp89( tmp6*b_[2] - tmp58*A_(2,2) + tmp59*A_(2,1) );
   const ET tmp90( tmp3*b_[2] - tmp57*A_(2,3) + tmp60*A_(2,0) );
   const ET tmp91( tmp7*b_[2] - tmp58*A_(2,3) + tmp60*A_(2,1) );
   const ET tmp92( tmp10*b_[2] - tmp59*A_(2,3) + tmp60*A_(2,2) );
   const ET tmp93( tmp4*b_[2] - tmp57*A_(2,4) + tmp61*A_(2,0) );
   const ET tmp94( tmp8*b_[2] - tmp58*A_(2,4) + tmp61*A_(2,1) );
   const ET tmp95( tmp11*b_[2] - tmp59*A_(2,4) + tmp61*A_(2,2) );
   const ET tmp96( tmp13*b_[2] - tmp60*A_(2,4) + tmp61*A_(2,3) );

   const ET tmp97 ( tmp67*A_(3,3) - tmp68*A_(3,2) + tmp69*A_(3,1) - tmp19*b_[3] );
   const ET tmp98 ( tmp67*A_(3,4) - tmp70*A_(3,2) + tmp71*A_(3,1) - tmp22*b_[3] );
   const ET tmp99 ( tmp68*A_(3,4) - tmp70*A_(3,3) + tmp72*A_(3,1) - tmp24*b_[3] );
   const ET tmp100( tmp69*A_(3,4) - tmp71*A_(3,3) + tmp72*A_(3,2) - tmp25*b_[3] );
   const ET tmp101( tmp67*A_(3,5) - tmp73*A_(3,2) + tmp74*A_(3,1) - tmp28*b_[3] );
   const ET tmp102( tmp68*A_(3,5) - tmp73*A_(3,3) + tmp75*A_(3,1) - tmp30*b_[3] );
   const ET tmp103( tmp69*A_(3,5) - tmp74*A_(3,3) + tmp75*A_(3,2) - tmp31*b_[3] );
   const ET tmp104( tmp70*A_(3,5) - tmp73*A_(3,4) + tmp76*A_(3,1) - tmp33*b_[3] );
   const ET tmp105( tmp71*A_(3,5) - tmp74*A_(3,4) + tmp76*A_(3,2) - tmp34*b_[3] );
   const ET tmp106( tmp72*A_(3,5) - tmp75*A_(3,4) + tmp76*A_(3,3) - tmp35*b_[3] );

   const ET tmp107( tmp77*A_(3,3) - tmp78*A_(3,2) + tmp18*b_[3] - tmp69*A_(3,0) );
   const ET tmp108( tmp77*A_(3,4) - tmp80*A_(3,2) + tmp21*b_[3] - tmp71*A_(3,0) );
   const ET tmp109( tmp78*A_(3,4) - tmp80*A_(3,3) + tmp23*b_[3] - tmp72*A_(3,0) );
   const ET tmp110( tmp79*A_(3,4) - tmp81*A_(3,3) + tmp24*b_[3] - tmp72*A_(3,1) );
   const ET tmp111( tmp77*A_(3,5) - tmp83*A_(3,2) + tmp27*b_[3] - tmp74*A_(3,0) );
   const ET tmp112( tmp78*A_(3,5) - tmp83*A_(3,3) + tmp29*b_[3] - tmp75*A_(3,0) );
   const ET tmp113( tmp79*A_(3,5) - tmp84*A_(3,3) + tmp30*b_[3] - tmp75*A_(3,1) );
   const ET tmp114( tmp80*A_(3,5) - tmp83*A_(3,4) + tmp32*b_[3] - tmp76*A_(3,0) );
   const ET tmp115( tmp81*A_(3,5) - tmp84*A_(3,4) + tmp33*b_[3] - tmp76*A_(3,1) );
   const ET tmp116( tmp82*A_(3,5) - tmp85*A_(3,4) + tmp34*b_[3] - tmp76*A_(3,2) );

   const ET tmp117( tmp87*A_(3,3) - tmp17*b_[3] + tmp78*A_(3,1) - tmp79*A_(3,0) );
   const ET tmp118( tmp87*A_(3,4) - tmp20*b_[3] + tmp80*A_(3,1) - tmp81*A_(3,0) );
   const ET tmp119( tmp88*A_(3,4) - tmp21*b_[3] + tmp80*A_(3,2) - tmp82*A_(3,0) );
   const ET tmp120( tmp89*A_(3,4) - tmp22*b_[3] + tmp81*A_(3,2) - tmp82*A_(3,1) );
   const ET tmp121( tmp87*A_(3,5) - tmp26*b_[3] + tmp83*A_(3,1) - tmp84*A_(3,0) );
   const ET tmp122( tmp88*A_(3,5) - tmp27*b_[3] + tmp83*A_(3,2) - tmp85*A_(3,0) );
   const ET tmp123( tmp89*A_(3,5) - tmp28*b_[3] + tmp84*A_(3,2) - tmp85*A_(3,1) );
   const ET tmp124( tmp90*A_(3,5) - tmp29*b_[3] + tmp83*A_(3,3) - tmp86*A_(3,0) );
   const ET tmp125( tmp91*A_(3,5) - tmp30*b_[3] + tmp84*A_(3,3) - tmp86*A_(3,1) );
   const ET tmp126( tmp92*A_(3,5) - tmp31*b_[3] + tmp85*A_(3,3) - tmp86*A_(3,2) );

   const ET tmp127( tmp16*b_[3] - tmp87*A_(3,2) + tmp88*A_(3,1) - tmp89*A_(3,0) );
   const ET tmp128( tmp17*b_[3] - tmp87*A_(3,3) + tmp90*A_(3,1) - tmp91*A_(3,0) );
   const ET tmp129( tmp18*b_[3] - tmp88*A_(3,3) + tmp90*A_(3,2) - tmp92*A_(3,0) );
   const ET tmp130( tmp19*b_[3] - tmp89*A_(3,3) + tmp91*A_(3,2) - tmp92*A_(3,1) );
   const ET tmp131( tmp20*b_[3] - tmp87*A_(3,4) + tmp93*A_(3,1) - tmp94*A_(3,0) );
   const ET tmp132( tmp21*b_[3] - tmp88*A_(3,4) + tmp93*A_(3,2) - tmp95*A_(3,0) );
   const ET tmp133( tmp22*b_[3] - tmp89*A_(3,4) + tmp94*A_(3,2) - tmp95*A_(3,1) );
   const ET tmp134( tmp23*b_[3] - tmp90*A_(3,4) + tmp93*A_(3,3) - tmp96*A_(3,0) );
   const ET tmp135( tmp24*b_[3] - tmp91*A_(3,4) + tmp94*A_(3,3) - tmp96*A_(3,1) );
   const ET tmp136( tmp25*b_[3] - tmp92*A_(3,4) + tmp95*A_(3,3) - tmp96*A_(3,2) );

   const ET tmp137( tmp97*A_(4,4) - tmp98*A_(4,3) + tmp99*A_(4,2) - tmp100*A_(4,1) + tmp40*b_[4] );
   const ET tmp138( tmp97*A_(4,5) - tmp101*A_(4,3) + tmp102*A_(4,2) - tmp103*A_(4,1) + tmp44*b_[4] );
   const ET tmp139( tmp98*A_(4,5) - tmp101*A_(4,4) + tmp104*A_(4,2) - tmp105*A_(4,1) + tmp47*b_[4] );
   const ET tmp140( tmp99*A_(4,5) - tmp102*A_(4,4) + tmp104*A_(4,3) - tmp106*A_(4,1) + tmp49*b_[4] );
   const ET tmp141( tmp100*A_(4,5) - tmp103*A_(4,4) + tmp105*A_(4,3) - tmp106*A_(4,2) + tmp50*b_[4] );

   const ET tmp142( tmp107*A_(4,4) - tmp108*A_(4,3) + tmp109*A_(4,2) - tmp39*b_[4] + tmp100*A_(4,0) );
   const ET tmp143( tmp107*A_(4,5) - tmp111*A_(4,3) + tmp112*A_(4,2) - tmp43*b_[4] + tmp103*A_(4,0) );
   const ET tmp144( tmp108*A_(4,5) - tmp111*A_(4,4) + tmp114*A_(4,2) - tmp46*b_[4] + tmp105*A_(4,0) );
   const ET tmp145( tmp109*A_(4,5) - tmp112*A_(4,4) + tmp114*A_(4,3) - tmp48*b_[4] + tmp106*A_(4,0) );
   const ET tmp146( tmp110*A_(4,5) - tmp113*A_(4,4) + tmp115*A_(4,3) - tmp49*b_[4] + tmp106*A_(4,1) );

   const ET tmp147( tmp117*A_(4,4) - tmp118*A_(4,3) + tmp38*b_[4] - tmp109*A_(4,1) + tmp110*A_(4,0) );
   const ET tmp148( tmp117*A_(4,5) - tmp121*A_(4,3) + tmp42*b_[4] - tmp112*A_(4,1) + tmp113*A_(4,0) );
   const ET tmp149( tmp118*A_(4,5) - tmp121*A_(4,4) + tmp45*b_[4] - tmp114*A_(4,1) + tmp115*A_(4,0) );
   const ET tmp150( tmp119*A_(4,5) - tmp122*A_(4,4) + tmp46*b_[4] - tmp114*A_(4,2) + tmp116*A_(4,0) );
   const ET tmp151( tmp120*A_(4,5) - tmp123*A_(4,4) + tmp47*b_[4] - tmp115*A_(4,2) + tmp116*A_(4,1) );

   const ET tmp152( tmp127*A_(4,4) - tmp37*b_[4] + tmp118*A_(4,2) - tmp119*A_(4,1) + tmp120*A_(4,0) );
   const ET tmp153( tmp127*A_(4,5) - tmp41*b_[4] + tmp121*A_(4,2) - tmp122*A_(4,1) + tmp123*A_(4,0) );
   const ET tmp154( tmp128*A_(4,5) - tmp42*b_[4] + tmp121*A_(4,3) - tmp124*A_(4,1) + tmp125*A_(4,0) );
   const ET tmp155( tmp129*A_(4,5) - tmp43*b_[4] + tmp122*A_(4,3) - tmp124*A_(4,2) + tmp126*A_(4,0) );
   const ET tmp156( tmp130*A_(4,5) - tmp44*b_[4] + tmp123*A_(4,3) - tmp125*A_(4,2) + tmp126*A_(4,1) );

   const ET tmp157( tmp36*b_[4] - tmp127*A_(4,3) + tmp128*A_(4,2) - tmp129*A_(4,1) + tmp130*A_(4,0) );
   const ET tmp158( tmp37*b_[4] - tmp127*A_(4,4) + tmp131*A_(4,2) - tmp132*A_(4,1) + tmp133*A_(4,0) );
   const ET tmp159( tmp38*b_[4] - tmp128*A_(4,4) + tmp131*A_(4,3) - tmp134*A_(4,1) + tmp135*A_(4,0) );
   const ET tmp160( tmp39*b_[4] - tmp129*A_(4,4) + tmp132*A_(4,3) - tmp134*A_(4,2) + tmp136*A_(4,0) );
   const ET tmp161( tmp40*b_[4] - tmp130*A_(4,4) + tmp133*A_(4,3) - tmp135*A_(4,2) + tmp136*A_(4,1) );

   x_[0] = ( tmp137*A_(5,5) - tmp138*A_(5,4) + tmp139*A_(5,3) - tmp140*A_(5,2) + tmp141*A_(5,1) - tmp56*b_[5] ) * invD;
   x_[1] = ( tmp142*A_(5,5) - tmp143*A_(5,4) + tmp144*A_(5,3) - tmp145*A_(5,2) + tmp55*b_[5] - tmp141*A_(5,0) ) * invD;
   x_[2] = ( tmp147*A_(5,5) - tmp148*A_(5,4) + tmp149*A_(5,3) - tmp54*b_[5] + tmp145*A_(5,1) - tmp146*A_(5,0) ) * invD;
   x_[3] = ( tmp152*A_(5,5) - tmp153*A_(5,4) + tmp53*b_[5] - tmp149*A_(5,2) + tmp150*A_(5,1) - tmp151*A_(5,0) ) * invD;
   x_[4] = ( tmp157*A_(5,5) - tmp52*b_[5] + tmp153*A_(5,3) - tmp154*A_(5,2) + tmp155*A_(5,1) - tmp156*A_(5,0) ) * invD;
   x_[5] = ( tmp51*b_[5] - tmp157*A_(5,4) + tmp158*A_(5,3) - tmp159*A_(5,2) + tmp160*A_(5,1) - tmp161*A_(5,0) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   using ET = ElementType_t<MT>;

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   const ET tmp1 ( A_(0,0)*A_(1,1) - A_(0,1)*A_(0,1) );
   const ET tmp2 ( A_(0,0)*A_(1,2) - A_(0,1)*A_(0,2) );
   const ET tmp3 ( A_(0,0)*A_(1,3) - A_(0,1)*A_(0,3) );
   const ET tmp4 ( A_(0,0)*A_(1,4) - A_(0,1)*A_(0,4) );
   const ET tmp5 ( A_(0,0)*A_(1,5) - A_(0,1)*A_(0,5) );
   const ET tmp6 ( A_(0,1)*A_(1,2) - A_(1,1)*A_(0,2) );
   const ET tmp7 ( A_(0,1)*A_(1,3) - A_(1,1)*A_(0,3) );
   const ET tmp8 ( A_(0,1)*A_(1,4) - A_(1,1)*A_(0,4) );
   const ET tmp9 ( A_(0,1)*A_(1,5) - A_(1,1)*A_(0,5) );
   const ET tmp10( A_(0,2)*A_(1,3) - A_(1,2)*A_(0,3) );
   const ET tmp11( A_(0,2)*A_(1,4) - A_(1,2)*A_(0,4) );
   const ET tmp12( A_(0,2)*A_(1,5) - A_(1,2)*A_(0,5) );
   const ET tmp13( A_(0,3)*A_(1,4) - A_(1,3)*A_(0,4) );
   const ET tmp14( A_(0,3)*A_(1,5) - A_(1,3)*A_(0,5) );
   const ET tmp15( A_(0,4)*A_(1,5) - A_(1,4)*A_(0,5) );

   const ET tmp16( tmp1*A_(2,2) - tmp2*A_(1,2) + tmp6*A_(0,2) );
   const ET tmp17( tmp1*A_(2,3) - tmp3*A_(1,2) + tmp7*A_(0,2) );
   const ET tmp18( tmp2*A_(2,3) - tmp3*A_(2,2) + tmp10*A_(0,2) );
   const ET tmp19( tmp6*A_(2,3) - tmp7*A_(2,2) + tmp10*A_(1,2) );
   const ET tmp20( tmp1*A_(2,4) - tmp4*A_(1,2) + tmp8*A_(0,2) );
   const ET tmp21( tmp2*A_(2,4) - tmp4*A_(2,2) + tmp11*A_(0,2) );
   const ET tmp22( tmp6*A_(2,4) - tmp8*A_(2,2) + tmp11*A_(1,2) );
   const ET tmp23( tmp3*A_(2,4) - tmp4*A_(2,3) + tmp13*A_(0,2) );
   const ET tmp24( tmp7*A_(2,4) - tmp8*A_(2,3) + tmp13*A_(1,2) );
   const ET tmp25( tmp10*A_(2,4) - tmp11*A_(2,3) + tmp13*A_(2,2) );
   const ET tmp26( tmp1*A_(2,5) - tmp5*A_(1,2) + tmp9*A_(0,2) );
   const ET tmp27( tmp2*A_(2,5) - tmp5*A_(2,2) + tmp12*A_(0,2) );
   const ET tmp28( tmp6*A_(2,5) - tmp9*A_(2,2) + tmp12*A_(1,2) );
   const ET tmp29( tmp3*A_(2,5) - tmp5*A_(2,3) + tmp14*A_(0,2) );
   const ET tmp30( tmp7*A_(2,5) - tmp9*A_(2,3) + tmp14*A_(1,2) );
   const ET tmp31( tmp10*A_(2,5) - tmp12*A_(2,3) + tmp14*A_(2,2) );
   const ET tmp32( tmp4*A_(2,5) - tmp5*A_(2,4) + tmp15*A_(0,2) );
   const ET tmp33( tmp8*A_(2,5) - tmp9*A_(2,4) + tmp15*A_(1,2) );
   const ET tmp34( tmp11*A_(2,5) - tmp12*A_(2,4) + tmp15*A_(2,2) );
   const ET tmp35( tmp13*A_(2,5) - tmp14*A_(2,4) + tmp15*A_(2,3) );

   const ET tmp36( tmp16*A_(3,3) - tmp17*A_(2,3) + tmp18*A_(1,3) - tmp19*A_(0,3) );
   const ET tmp37( tmp16*A_(3,4) - tmp20*A_(2,3) + tmp21*A_(1,3) - tmp22*A_(0,3) );
   const ET tmp38( tmp17*A_(3,4) - tmp20*A_(3,3) + tmp23*A_(1,3) - tmp24*A_(0,3) );
   const ET tmp39( tmp18*A_(3,4) - tmp21*A_(3,3) + tmp23*A_(2,3) - tmp25*A_(0,3) );
   const ET tmp40( tmp19*A_(3,4) - tmp22*A_(3,3) + tmp24*A_(2,3) - tmp25*A_(1,3) );
   const ET tmp41( tmp16*A_(3,5) - tmp26*A_(2,3) + tmp27*A_(1,3) - tmp28*A_(0,3) );
   const ET tmp42( tmp17*A_(3,5) - tmp26*A_(3,3) + tmp29*A_(1,3) - tmp30*A_(0,3) );
   const ET tmp43( tmp18*A_(3,5) - tmp27*A_(3,3) + tmp29*A_(2,3) - tmp31*A_(0,3) );
   const ET tmp44( tmp19*A_(3,5) - tmp28*A_(3,3) + tmp30*A_(2,3) - tmp31*A_(1,3) );
   const ET tmp45( tmp20*A_(3,5) - tmp26*A_(3,4) + tmp32*A_(1,3) - tmp33*A_(0,3) );
   const ET tmp46( tmp21*A_(3,5) - tmp27*A_(3,4) + tmp32*A_(2,3) - tmp34*A_(0,3) );
   const ET tmp47( tmp22*A_(3,5) - tmp28*A_(3,4) + tmp33*A_(2,3) - tmp34*A_(1,3) );
   const ET tmp48( tmp23*A_(3,5) - tmp29*A_(3,4) + tmp32*A_(3,3) - tmp35*A_(0,3) );
   const ET tmp49( tmp24*A_(3,5) - tmp30*A_(3,4) + tmp33*A_(3,3) - tmp35*A_(1,3) );
   const ET tmp50( tmp25*A_(3,5) - tmp31*A_(3,4) + tmp34*A_(3,3) - tmp35*A_(2,3) );

   const ET tmp51( tmp36*A_(4,4) - tmp37*A_(3,4) + tmp38*A_(2,4) - tmp39*A_(1,4) + tmp40*A_(0,4) );
   const ET tmp52( tmp36*A_(4,5) - tmp41*A_(3,4) + tmp42*A_(2,4) - tmp43*A_(1,4) + tmp44*A_(0,4) );
   const ET tmp53( tmp37*A_(4,5) - tmp41*A_(4,4) + tmp45*A_(2,4) - tmp46*A_(1,4) + tmp47*A_(0,4) );
   const ET tmp54( tmp38*A_(4,5) - tmp42*A_(4,4) + tmp45*A_(3,4) - tmp48*A_(1,4) + tmp49*A_(0,4) );
   const ET tmp55( tmp39*A_(4,5) - tmp43*A_(4,4) + tmp46*A_(3,4) - tmp48*A_(2,4) + tmp50*A_(0,4) );
   const ET tmp56( tmp40*A_(4,5) - tmp44*A_(4,4) + tmp47*A_(3,4) - tmp49*A_(2,4) + tmp50*A_(1,4) );

   const ET D( tmp51*A_(5,5) - tmp52*A_(4,5) + tmp53*A_(3,5) - tmp54*A_(2,5) + tmp55*A_(1,5) - tmp56*A_(0,5) );

   if( !isDivisor( D ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const ET invD( inv( D ) );

   const ET tmp57( A_(0,0)*b_[1] - A_(0,1)*b_[0] );
   const ET tmp58( A_(0,1)*b_[1] - A_(1,1)*b_[0] );
   const ET tmp59( A_(0,2)*b_[1] - A_(1,2)*b_[0] );
   const ET tmp60( A_(0,3)*b_[1] - A_(1,3)*b_[0] );
   const ET tmp61( A_(0,4)*b_[1] - A_(1,4)*b_[0] );

   const ET tmp62( A_(1,1)*b_[0] - A_(0,1)*b_[1] );
   const ET tmp63( A_(1,2)*b_[0] - A_(0,2)*b_[1] );
   const ET tmp64( A_(1,3)*b_[0] - A_(0,3)*b_[1] );
   const ET tmp65( A_(1,4)*b_[0] - A_(0,4)*b_[1] );
   const ET tmp66( A_(1,5)*b_[0] - A_(0,5)*b_[1] );

   const ET tmp67( tmp62*A_(2,2) - tmp63*A_(1,2) + tmp6*b_[2] );
   const ET tmp68( tmp62*A_(2,3) - tmp64*A_(1,2) + tmp7*b_[2] );
   const ET tmp69( tmp63*A_(2,3) - tmp64*A_(2,2) + tmp10*b_[2] );
   const ET tmp70( tmp62*A_(2,4) - tmp65*A_(1,2) + tmp8*b_[2] );
   const ET tmp71( tmp63*A_(2,4) - tmp65*A_(2,2) + tmp11*b_[2] );
   const ET tmp72( tmp64*A_(2,4) - tmp65*A_(2,3) + tmp13*b_[2] );
   const ET tmp73( tmp62*A_(2,5) - tmp66*A_(1,2) + tmp9*b_[2] );
   const ET tmp74( tmp63*A_(2,5) - tmp66*A_(2,2) + tmp12*b_[2] );
   const ET tmp75( tmp64*A_(2,5) - tmp66*A_(2,3) + tmp14*b_[2] );
   const ET tmp76( tmp65*A_(2,5) - tmp66*A_(2,4) + tmp15*b_[2] );

   const ET tmp77( tmp57*A_(2,2) - tmp2*b_[2] + tmp63*A_(0,2) );
   const ET tmp78( tmp57*A_(2,3) - tmp3*b_[2] + tmp64*A_(0,2) );
   const ET tmp79( tmp58*A_(2,3) - tmp7*b_[2] + tmp64*A_(1,2) );
   const ET tmp80( tmp57*A_(2,4) - tmp4*b_[2] + tmp65*A_(0,2) );
   const ET tmp81( tmp58*A_(2,4) - tmp8*b_[2] + tmp65*A_(1,2) );
   const ET tmp82( tmp59*A_(2,4) - tmp11*b_[2] + tmp65*A_(2,2) );
   const ET tmp83( tmp57*A_(2,5) - tmp5*b_[2] + tmp66*A_(0,2) );
   const ET tmp84( tmp58*A_(2,5) - tmp9*b_[2] + tmp66*A_(1,2) );
   const ET tmp85( tmp59*A_(2,5) - tmp12*b_[2] + tmp66*A_(2,2) );
   const ET tmp86( tmp60*A_(2,5) - tmp14*b_[2] + tmp66*A_(2,3) );

   const ET tmp87( tmp1*b_[2] - tmp57*A_(1,2) + tmp58*A_(0,2) );
   const ET tmp88( tmp2*b_[2] - tmp57*A_(2,2) + tmp59*A_(0,2) );
   const ET tmp89( tmp6*b_[2] - tmp58*A_(2,2) + tmp59*A_(1,2) );
   const ET tmp90( tmp3*b_[2] - tmp57*A_(2,3) + tmp60*A_(0,2) );
   const ET tmp91( tmp7*b_[2] - tmp58*A_(2,3) + tmp60*A_(1,2) );
   const ET tmp92( tmp10*b_[2] - tmp59*A_(2,3) + tmp60*A_(2,2) );
   const ET tmp93( tmp4*b_[2] - tmp57*A_(2,4) + tmp61*A_(0,2) );
   const ET tmp94( tmp8*b_[2] - tmp58*A_(2,4) + tmp61*A_(1,2) );
   const ET tmp95( tmp11*b_[2] - tmp59*A_(2,4) + tmp61*A_(2,2) );
   const ET tmp96( tmp13*b_[2] - tmp60*A_(2,4) + tmp61*A_(2,3) );

   const ET tmp97 ( tmp67*A_(3,3) - tmp68*A_(2,3) + tmp69*A_(1,3) - tmp19*b_[3] );
   const ET tmp98 ( tmp67*A_(3,4) - tmp70*A_(2,3) + tmp71*A_(1,3) - tmp22*b_[3] );
   const ET tmp99 ( tmp68*A_(3,4) - tmp70*A_(3,3) + tmp72*A_(1,3) - tmp24*b_[3] );
   const ET tmp100( tmp69*A_(3,4) - tmp71*A_(3,3) + tmp72*A_(2,3) - tmp25*b_[3] );
   const ET tmp101( tmp67*A_(3,5) - tmp73*A_(2,3) + tmp74*A_(1,3) - tmp28*b_[3] );
   const ET tmp102( tmp68*A_(3,5) - tmp73*A_(3,3) + tmp75*A_(1,3) - tmp30*b_[3] );
   const ET tmp103( tmp69*A_(3,5) - tmp74*A_(3,3) + tmp75*A_(2,3) - tmp31*b_[3] );
   const ET tmp104( tmp70*A_(3,5) - tmp73*A_(3,4) + tmp76*A_(1,3) - tmp33*b_[3] );
   const ET tmp105( tmp71*A_(3,5) - tmp74*A_(3,4) + tmp76*A_(2,3) - tmp34*b_[3] );
   const ET tmp106( tmp72*A_(3,5) - tmp75*A_(3,4) + tmp76*A_(3,3) - tmp35*b_[3] );

   const ET tmp107( tmp77*A_(3,3) - tmp78*A_(2,3) + tmp18*b_[3] - tmp69*A_(0,3) );
   const ET tmp108( tmp77*A_(3,4) - tmp80*A_(2,3) + tmp21*b_[3] - tmp71*A_(0,3) );
   const ET tmp109( tmp78*A_(3,4) - tmp80*A_(3,3) + tmp23*b_[3] - tmp72*A_(0,3) );
   const ET tmp110( tmp79*A_(3,4) - tmp81*A_(3,3) + tmp24*b_[3] - tmp72*A_(1,3) );
   const ET tmp111( tmp77*A_(3,5) - tmp83*A_(2,3) + tmp27*b_[3] - tmp74*A_(0,3) );
   const ET tmp112( tmp78*A_(3,5) - tmp83*A_(3,3) + tmp29*b_[3] - tmp75*A_(0,3) );
   const ET tmp113( tmp79*A_(3,5) - tmp84*A_(3,3) + tmp30*b_[3] - tmp75*A_(1,3) );
   const ET tmp114( tmp80*A_(3,5) - tmp83*A_(3,4) + tmp32*b_[3] - tmp76*A_(0,3) );
   const ET tmp115( tmp81*A_(3,5) - tmp84*A_(3,4) + tmp33*b_[3] - tmp76*A_(1,3) );
   const ET tmp116( tmp82*A_(3,5) - tmp85*A_(3,4) + tmp34*b_[3] - tmp76*A_(2,3) );

   const ET tmp117( tmp87*A_(3,3) - tmp17*b_[3] + tmp78*A_(1,3) - tmp79*A_(0,3) );
   const ET tmp118( tmp87*A_(3,4) - tmp20*b_[3] + tmp80*A_(1,3) - tmp81*A_(0,3) );
   const ET tmp119( tmp88*A_(3,4) - tmp21*b_[3] + tmp80*A_(2,3) - tmp82*A_(0,3) );
   const ET tmp120( tmp89*A_(3,4) - tmp22*b_[3] + tmp81*A_(2,3) - tmp82*A_(1,3) );
   const ET tmp121( tmp87*A_(3,5) - tmp26*b_[3] + tmp83*A_(1,3) - tmp84*A_(0,3) );
   const ET tmp122( tmp88*A_(3,5) - tmp27*b_[3] + tmp83*A_(2,3) - tmp85*A_(0,3) );
   const ET tmp123( tmp89*A_(3,5) - tmp28*b_[3] + tmp84*A_(2,3) - tmp85*A_(1,3) );
   const ET tmp124( tmp90*A_(3,5) - tmp29*b_[3] + tmp83*A_(3,3) - tmp86*A_(0,3) );
   const ET tmp125( tmp91*A_(3,5) - tmp30*b_[3] + tmp84*A_(3,3) - tmp86*A_(1,3) );
   const ET tmp126( tmp92*A_(3,5) - tmp31*b_[3] + tmp85*A_(3,3) - tmp86*A_(2,3) );

   const ET tmp127( tmp16*b_[3] - tmp87*A_(2,3) + tmp88*A_(1,3) - tmp89*A_(0,3) );
   const ET tmp128( tmp17*b_[3] - tmp87*A_(3,3) + tmp90*A_(1,3) - tmp91*A_(0,3) );
   const ET tmp129( tmp18*b_[3] - tmp88*A_(3,3) + tmp90*A_(2,3) - tmp92*A_(0,3) );
   const ET tmp130( tmp19*b_[3] - tmp89*A_(3,3) + tmp91*A_(2,3) - tmp92*A_(1,3) );
   const ET tmp131( tmp20*b_[3] - tmp87*A_(3,4) + tmp93*A_(1,3) - tmp94*A_(0,3) );
   const ET tmp132( tmp21*b_[3] - tmp88*A_(3,4) + tmp93*A_(2,3) - tmp95*A_(0,3) );
   const ET tmp133( tmp22*b_[3] - tmp89*A_(3,4) + tmp94*A_(2,3) - tmp95*A_(1,3) );
   const ET tmp134( tmp23*b_[3] - tmp90*A_(3,4) + tmp93*A_(3,3) - tmp96*A_(0,3) );
   const ET tmp135( tmp24*b_[3] - tmp91*A_(3,4) + tmp94*A_(3,3) - tmp96*A_(1,3) );
   const ET tmp136( tmp25*b_[3] - tmp92*A_(3,4) + tmp95*A_(3,3) - tmp96*A_(2,3) );

   const ET tmp137( tmp97*A_(4,4) - tmp98*A_(3,4) + tmp99*A_(2,4) - tmp100*A_(1,4) + tmp40*b_[4] );
   const ET tmp138( tmp97*A_(4,5) - tmp101*A_(3,4) + tmp102*A_(2,4) - tmp103*A_(1,4) + tmp44*b_[4] );
   const ET tmp139( tmp98*A_(4,5) - tmp101*A_(4,4) + tmp104*A_(2,4) - tmp105*A_(1,4) + tmp47*b_[4] );
   const ET tmp140( tmp99*A_(4,5) - tmp102*A_(4,4) + tmp104*A_(3,4) - tmp106*A_(1,4) + tmp49*b_[4] );
   const ET tmp141( tmp100*A_(4,5) - tmp103*A_(4,4) + tmp105*A_(3,4) - tmp106*A_(2,4) + tmp50*b_[4] );

   const ET tmp142( tmp107*A_(4,4) - tmp108*A_(3,4) + tmp109*A_(2,4) - tmp39*b_[4] + tmp100*A_(0,4) );
   const ET tmp143( tmp107*A_(4,5) - tmp111*A_(3,4) + tmp112*A_(2,4) - tmp43*b_[4] + tmp103*A_(0,4) );
   const ET tmp144( tmp108*A_(4,5) - tmp111*A_(4,4) + tmp114*A_(2,4) - tmp46*b_[4] + tmp105*A_(0,4) );
   const ET tmp145( tmp109*A_(4,5) - tmp112*A_(4,4) + tmp114*A_(3,4) - tmp48*b_[4] + tmp106*A_(0,4) );
   const ET tmp146( tmp110*A_(4,5) - tmp113*A_(4,4) + tmp115*A_(3,4) - tmp49*b_[4] + tmp106*A_(1,4) );

   const ET tmp147( tmp117*A_(4,4) - tmp118*A_(3,4) + tmp38*b_[4] - tmp109*A_(1,4) + tmp110*A_(0,4) );
   const ET tmp148( tmp117*A_(4,5) - tmp121*A_(3,4) + tmp42*b_[4] - tmp112*A_(1,4) + tmp113*A_(0,4) );
   const ET tmp149( tmp118*A_(4,5) - tmp121*A_(4,4) + tmp45*b_[4] - tmp114*A_(1,4) + tmp115*A_(0,4) );
   const ET tmp150( tmp119*A_(4,5) - tmp122*A_(4,4) + tmp46*b_[4] - tmp114*A_(2,4) + tmp116*A_(0,4) );
   const ET tmp151( tmp120*A_(4,5) - tmp123*A_(4,4) + tmp47*b_[4] - tmp115*A_(2,4) + tmp116*A_(1,4) );

   const ET tmp152( tmp127*A_(4,4) - tmp37*b_[4] + tmp118*A_(2,4) - tmp119*A_(1,4) + tmp120*A_(0,4) );
   const ET tmp153( tmp127*A_(4,5) - tmp41*b_[4] + tmp121*A_(2,4) - tmp122*A_(1,4) + tmp123*A_(0,4) );
   const ET tmp154( tmp128*A_(4,5) - tmp42*b_[4] + tmp121*A_(3,4) - tmp124*A_(1,4) + tmp125*A_(0,4) );
   const ET tmp155( tmp129*A_(4,5) - tmp43*b_[4] + tmp122*A_(3,4) - tmp124*A_(2,4) + tmp126*A_(0,4) );
   const ET tmp156( tmp130*A_(4,5) - tmp44*b_[4] + tmp123*A_(3,4) - tmp125*A_(2,4) + tmp126*A_(1,4) );

   const ET tmp157( tmp36*b_[4] - tmp127*A_(3,4) + tmp128*A_(2,4) - tmp129*A_(1,4) + tmp130*A_(0,4) );
   const ET tmp158( tmp37*b_[4] - tmp127*A_(4,4) + tmp131*A_(2,4) - tmp132*A_(1,4) + tmp133*A_(0,4) );
   const ET tmp159( tmp38*b_[4] - tmp128*A_(4,4) + tmp131*A_(3,4) - tmp134*A_(1,4) + tmp135*A_(0,4) );
   const ET tmp160( tmp39*b_[4] - tmp129*A_(4,4) + tmp132*A_(3,4) - tmp134*A_(2,4) + tmp136*A_(0,4) );
   const ET tmp161( tmp40*b_[4] - tmp130*A_(4,4) + tmp133*A_(3,4) - tmp135*A_(2,4) + tmp136*A_(1,4) );

   x_[0] = ( tmp137*A_(5,5) - tmp138*A_(4,5) + tmp139*A_(3,5) - tmp140*A_(2,5) + tmp141*A_(1,5) - tmp56*b_[5] ) * invD;
   x_[1] = ( tmp142*A_(5,5) - tmp143*A_(4,5) + tmp144*A_(3,5) - tmp145*A_(2,5) + tmp55*b_[5] - tmp141*A_(0,5) ) * invD;
   x_[2] = ( tmp147*A_(5,5) - tmp148*A_(4,5) + tmp149*A_(3,5) - tmp54*b_[5] + tmp145*A_(1,5) - tmp146*A_(0,5) ) * invD;
   x_[3] = ( tmp152*A_(5,5) - tmp153*A_(4,5) + tmp53*b_[5] - tmp149*A_(2,5) + tmp150*A_(1,5) - tmp151*A_(0,5) ) * invD;
   x_[4] = ( tmp157*A_(5,5) - tmp52*b_[5] + tmp153*A_(3,5) - tmp154*A_(2,5) + tmp155*A_(1,5) - tmp156*A_(0,5) ) * invD;
   x_[5] = ( tmp51*b_[5] - tmp157*A_(4,5) + tmp158*A_(3,5) - tmp159*A_(2,5) + tmp160*A_(1,5) - tmp161*A_(0,5) ) * invD;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = ( b_[0] ) / A_(0,0);
   x_[1] = ( b_[1] - A_(1,0)*x_[0] ) / A_(1,1);
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] ) / A_(2,2);
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] ) / A_(3,3);
   x_[4] = ( b_[4] - A_(4,0)*x_[0] - A_(4,1)*x_[1] - A_(4,2)*x_[2] - A_(4,3)*x_[3] ) / A_(4,4);
   x_[5] = ( b_[5] - A_(5,0)*x_[0] - A_(5,1)*x_[1] - A_(5,2)*x_[2] - A_(5,3)*x_[3] - A_(5,4)*x_[4] ) / A_(5,5);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[0] = ( b_[0] );
   x_[1] = ( b_[1] - A_(1,0)*x_[0] );
   x_[2] = ( b_[2] - A_(2,0)*x_[0] - A_(2,1)*x_[1] );
   x_[3] = ( b_[3] - A_(3,0)*x_[0] - A_(3,1)*x_[1] - A_(3,2)*x_[2] );
   x_[4] = ( b_[4] - A_(4,0)*x_[0] - A_(4,1)*x_[1] - A_(4,2)*x_[2] - A_(4,3)*x_[3] );
   x_[5] = ( b_[5] - A_(5,0)*x_[0] - A_(5,1)*x_[1] - A_(5,2)*x_[2] - A_(5,3)*x_[3] - A_(5,4)*x_[4] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[5] = ( b_[5] ) / A_(5,5);
   x_[4] = ( b_[4] - A_(4,5)*x_[5] ) / A_(4,4);
   x_[3] = ( b_[3] - A_(3,4)*x_[4] - A_(3,5)*x_[5] ) / A_(3,3);
   x_[2] = ( b_[2] - A_(2,3)*x_[3] - A_(2,4)*x_[4] - A_(2,5)*x_[5] ) / A_(2,2);
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] - A_(1,4)*x_[4] - A_(1,5)*x_[5] ) / A_(1,1);
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] - A_(0,4)*x_[4] - A_(0,5)*x_[5] ) / A_(0,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   resize( x_, b_.size() );

   x_[5] = ( b_[5] );
   x_[4] = ( b_[4] - A_(4,5)*x_[5] );
   x_[3] = ( b_[3] - A_(3,4)*x_[4] - A_(3,5)*x_[5] );
   x_[2] = ( b_[2] - A_(2,3)*x_[3] - A_(2,4)*x_[4] - A_(2,5)*x_[5] );
   x_[1] = ( b_[1] - A_(1,2)*x_[2] - A_(1,3)*x_[3] - A_(1,4)*x_[4] - A_(1,5)*x_[5] );
   x_[0] = ( b_[0] - A_(0,1)*x_[1] - A_(0,2)*x_[2] - A_(0,3)*x_[3] - A_(0,4)*x_[4] - A_(0,5)*x_[5] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solve6x6( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*b).size()    == 6UL, "Invalid vector size detected"       );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, b_.size() );

   x_[0] = b_[0] / A_(0,0);
   x_[1] = b_[1] / A_(1,1);
   x_[2] = b_[2] / A_(2,2);
   x_[3] = b_[3] / A_(3,3);
   x_[4] = b_[4] / A_(4,4);
   x_[5] = b_[5] / A_(5,5);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< ( IsGeneral_v<MT1> || IsSymmetric_v<MT1> || IsHermitian_v<MT1> ) && !IsDiagonal_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   resize( *X, (*B).rows(), (*B).columns() );
   const ResultType_t<MT1> invA( inv( *A ) );
   smpAssign( *X, invA * (*B) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );
   const ET invD5( inv( A_(5,5) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = ( B_(0,j) ) * invD0;
      X_(1,j) = ( B_(1,j) - A_(1,0)*X_(0,j) ) * invD1;
      X_(2,j) = ( B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j) ) * invD2;
      X_(3,j) = ( B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j) ) * invD3;
      X_(4,j) = ( B_(4,j) - A_(4,0)*X_(0,j) - A_(4,1)*X_(1,j) - A_(4,2)*X_(2,j) - A_(4,3)*X_(3,j) ) * invD4;
      X_(5,j) = ( B_(5,j) - A_(5,0)*X_(0,j) - A_(5,1)*X_(1,j) - A_(5,2)*X_(2,j) - A_(5,3)*X_(3,j) - A_(5,4)*X_(4,j) ) * invD5;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j);
      X_(1,j) = B_(1,j) - A_(1,0)*X_(0,j);
      X_(2,j) = B_(2,j) - A_(2,0)*X_(0,j) - A_(2,1)*X_(1,j);
      X_(3,j) = B_(3,j) - A_(3,0)*X_(0,j) - A_(3,1)*X_(1,j) - A_(3,2)*X_(2,j);
      X_(4,j) = B_(4,j) - A_(4,0)*X_(0,j) - A_(4,1)*X_(1,j) - A_(4,2)*X_(2,j) - A_(4,3)*X_(3,j);
      X_(5,j) = B_(5,j) - A_(5,0)*X_(0,j) - A_(5,1)*X_(1,j) - A_(5,2)*X_(2,j) - A_(5,3)*X_(3,j) - A_(5,4)*X_(4,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );
   const ET invD5( inv( A_(5,5) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(5,j) = ( B_(5,j) ) * invD5;
      X_(4,j) = ( B_(4,j) - A_(4,5)*X_(5,j) ) * invD4;
      X_(3,j) = ( B_(3,j) - A_(3,4)*X_(4,j) - A_(3,5)*X_(5,j) ) * invD3;
      X_(2,j) = ( B_(2,j) - A_(2,3)*X_(3,j) - A_(2,4)*X_(4,j) - A_(2,5)*X_(5,j) ) * invD2;
      X_(1,j) = ( B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j) - A_(1,4)*X_(4,j) - A_(1,5)*X_(5,j) ) * invD1;
      X_(0,j) = ( B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j) - A_(0,4)*X_(4,j) - A_(0,5)*X_(5,j) ) * invD0;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t j=0UL; j<N; ++j ) {
      X_(5,j) = B_(5,j);
      X_(4,j) = B_(4,j) - A_(4,5)*X_(5,j);
      X_(3,j) = B_(3,j) - A_(3,4)*X_(4,j) - A_(3,5)*X_(5,j);
      X_(2,j) = B_(2,j) - A_(2,3)*X_(3,j) - A_(2,4)*X_(4,j) - A_(2,5)*X_(5,j);
      X_(1,j) = B_(1,j) - A_(1,2)*X_(2,j) - A_(1,3)*X_(3,j) - A_(1,4)*X_(4,j) - A_(1,5)*X_(5,j);
      X_(0,j) = B_(0,j) - A_(0,1)*X_(1,j) - A_(0,2)*X_(2,j) - A_(0,3)*X_(3,j) - A_(0,4)*X_(4,j) - A_(0,5)*X_(5,j);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ 6 \times 6 \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for a \f$ 6 \times 6 \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not a \f$ 6 \times 6 \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solve6x6( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( (*A).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*A).columns() == 6UL, "Invalid number of columns detected" );
   BLAZE_INTERNAL_ASSERT( (*B).rows()    == 6UL, "Invalid number of rows detected"    );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( A_(0,0)*A_(1,1)*A_(2,2)*A_(3,3)*A_(4,4)*A_(5,5) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   const ET invD0( inv( A_(0,0) ) );
   const ET invD1( inv( A_(1,1) ) );
   const ET invD2( inv( A_(2,2) ) );
   const ET invD3( inv( A_(3,3) ) );
   const ET invD4( inv( A_(4,4) ) );
   const ET invD5( inv( A_(5,5) ) );

   for( size_t j=0UL; j<N; ++j ) {
      X_(0,j) = B_(0,j) * invD0;
      X_(1,j) = B_(1,j) * invD1;
      X_(2,j) = B_(2,j) * invD2;
      X_(3,j) = B_(3,j) * invD3;
      X_(4,j) = B_(4,j) * invD4;
      X_(5,j) = B_(5,j) * invD5;
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING NxN LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ general linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given general system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not an \f$ N \times N \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsGeneral_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   using RT = ResultType_t<MT>;
   using OT = OppositeType_t<MT>;

   const size_t N( (*b).size() );

   RemoveAdaptor_t< If_t<SO,RT,OT> > Atmp( A );

   resize( *x, N );
   smpAssign( *x, *b );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   gesv( Atmp, *x, ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ symmetric linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given symmetric system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not an \f$ N \times N \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   using RT = ResultType_t<MT>;
   using OT = OppositeType_t<MT>;

   const size_t N( (*b).size() );

   RemoveAdaptor_t< If_t<SO,RT,OT> > Atmp( A );

   resize( *x, N );
   smpAssign( *x, *b );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   sysv( Atmp, *x, 'L', ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ Hermitian linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given Hermitian system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not an \f$ N \times N \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsHermitian_v<MT> && !IsSymmetric_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   using RT = ResultType_t<MT>;
   using OT = OppositeType_t<MT>;

   const size_t N( (*b).size() );

   RemoveAdaptor_t< If_t<SO,RT,OT> > Atmp( A );

   resize( *x, N );
   smpAssign( *x, *b );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   hesv( Atmp, *x, 'L', ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ lower triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsLower_v<MT> && !IsUniLower_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t N( (*b).size() );

   resize( x_, N );

   for( size_t i=0UL; i<N; ++i ) {
      x_[i] = b_[i];
      for( size_t j=0UL; j<i; ++j )
         x_[i] -= A_(i,j)*x_[j];
      x_[i] *= inv( A_(i,i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ lower unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given lower unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniLower_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t N( (*b).size() );

   resize( x_, N );

   for( size_t i=0UL; i<N; ++i ) {
      x_[i] = b_[i];
      for( size_t j=0UL; j<i; ++j )
         x_[i] -= A_(i,j)*x_[j];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ upper triangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper triangular system matrix, \a x is the given
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUpper_v<MT> && !IsUniUpper_v<MT> && !IsDiagonal_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t N( (*b).size() );

   resize( x_, N );

   for( size_t i=N-1UL; i<N; --i ) {
      x_[i] = b_[i];
      for( size_t j=i+1UL; j<N; ++j )
         x_[i] -= A_(i,j)*x_[j];
      x_[i] *= inv( A_(i,i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ upper unitriangular linear system of equations
//        (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given upper unitriangular system matrix, \a x is the
// solution vector, and \a b is the given right-hand side vector. The function fails if the
// given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t N( (*b).size() );

   resize( x_, N );

   for( size_t i=N-1UL; i<N; --i ) {
      x_[i] = b_[i];
      for( size_t j=i+1UL; j<N; ++j )
         x_[i] -= A_(i,j)*x_[j];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ diagonal linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param x The dense solution vector.
// \param b The dense right-hand side vector.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * x = b \f$, where \a A is the given diagonal system matrix, \a x is the solution
// vector, and \a b is the given right-hand side vector. The function fails if the given matrix
// is not an \f$ N \times N \f$ square matrix or singular. In this case a \a std::runtime_error
// exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2      // Transpose flag of the right-hand side vector
        , EnableIf_t< IsDiagonal_v<MT> >* = nullptr >
void solveNxN( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*b).size(), "Invalid matrix and vector sizes" );

   const size_t N( (*b).size() );

   CompositeType_t<MT> A_( *A );
   VT1& x_( *x );
   const VT2& b_( *b );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   resize( x_, N );

   for( size_t i=0UL; i<N; ++i ) {
      x_[i] = b_[i] / A_(i,i);
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ general linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense general system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsGeneral_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using MT4 = RemoveAdaptor_t< If_t< SO1, ResultType_t<MT1>, OppositeType_t<MT1> > >;
   using MT5 = RemoveAdaptor_t< If_t< SO3, ResultType_t<MT3>, OppositeType_t<MT3> > >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT5 );

   const size_t N( (*A).rows() );

   MT4 Atmp( A );
   MT5 Xtmp( B );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   gesv( Atmp, Xtmp, ipiv.get() );

   resize( *X, Xtmp.rows(), Xtmp.columns() );
   smpAssign( *X, Xtmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ symmetric linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense symmetric system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsSymmetric_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using MT4 = RemoveAdaptor_t< If_t< SO1, ResultType_t<MT1>, OppositeType_t<MT1> > >;
   using MT5 = RemoveAdaptor_t< If_t< SO3, ResultType_t<MT3>, OppositeType_t<MT3> > >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT5 );

   const size_t N( (*A).rows() );

   MT4 Atmp( A );
   MT5 Xtmp( B );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   sysv( Atmp, Xtmp, 'L', ipiv.get() );

   resize( *X, Xtmp.rows(), Xtmp.columns() );
   smpAssign( *X, Xtmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ Hermitian linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense Hermitian system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsHermitian_v<MT1> && !IsSymmetric_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using MT4 = RemoveAdaptor_t< If_t< SO1, ResultType_t<MT1>, OppositeType_t<MT1> > >;
   using MT5 = RemoveAdaptor_t< If_t< SO3, ResultType_t<MT3>, OppositeType_t<MT3> > >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT4 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT5 );

   const size_t N( (*A).rows() );

   MT4 Atmp( A );
   MT5 Xtmp( B );

   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[N] );

   hesv( Atmp, Xtmp, 'L', ipiv.get() );

   resize( *X, Xtmp.rows(), Xtmp.columns() );
   smpAssign( *X, Xtmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ lower triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsLower_v<MT1> && !IsUniLower_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t i=0UL; i<M; ++i ) {
      const ET invD( inv( A_(i,i) ) );
      for( size_t j=0UL; j<N; ++j ) {
         X_(i,j) = B_(i,j);
         for( size_t k=0UL; k<i; ++k )
            X_(i,j) -= A_(i,k)*X_(k,j);
         X_(i,j) *= invD;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ lower unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense lower unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniLower_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t i=0UL; i<M; ++i ) {
      for( size_t j=0UL; j<N; ++j ) {
         X_(i,j) = B_(i,j);
         for( size_t k=0UL; k<i; ++k )
            X_(i,j) -= A_(i,k)*X_(k,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ upper triangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper triangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUpper_v<MT1> && !IsUniUpper_v<MT1> && !IsDiagonal_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t i=M-1UL; i<M; --i ) {
      const ET invD( inv( A_(i,i) ) );
      for( size_t j=0; j<N; ++j ) {
         X_(i,j) = B_(i,j);
         for( size_t k=i+1UL; k<M; ++k )
            X_(i,j) -= A_(i,k)*X_(k,j);
         X_(i,j) *= invD;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ upper unitriangular linear system of equations
//        (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense upper unitriangular system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsUniUpper_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t i=M-1UL; i<M; --i ) {
      for( size_t j=0; j<N; ++j ) {
         X_(i,j) = B_(i,j);
         for( size_t k=i+1UL; k<M; ++k )
            X_(i,j) -= A_(i,k)*X_(k,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Solving the given \f$ N \times N \f$ diagonal linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The dense diagonal system matrix.
// \param X The dense solution matrix.
// \param B The dense right-hand side matrix.
// \return void
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for an \f$ N \times N \f$ linear system of equations (LSE)
// \f$ A * X = B \f$, where \a A is the given system matrix, the columns of \a X are the solution
// vectors, and the columns of \a B are the given right-hand side vectors. The function fails if
// the given matrix is not an \f$ N \times N \f$ square matrix or singular. In this case a
// \a std::runtime_error exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3      // Storage order of the right-hand side matrix
        , EnableIf_t< IsDiagonal_v<MT1> >* = nullptr >
void solveNxN( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   BLAZE_INTERNAL_ASSERT( isSquare( *A ), "Non-square matrix detected" );
   BLAZE_INTERNAL_ASSERT( (*A).rows() == (*B).rows(), "Invalid matrix and vector sizes" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).rows() == (*X).rows(), "Invalid number of rows detected" );
   BLAZE_INTERNAL_ASSERT( IsResizable_v<MT2> || (*B).columns() == (*X).columns(), "Invalid number of columns detected" );

   using ET = ElementType_t<MT1>;

   CompositeType_t<MT1> A_( *A );
   MT2& X_( *X );
   const MT3& B_( *B );

   if( !isDivisor( trace( A_ ) ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Solving LSE with singular system matrix failed" );
   }

   const size_t M( B_.rows()    );
   const size_t N( B_.columns() );

   resize( X_, M, N );

   for( size_t i=0UL; i<M; ++i ) {
      const ET invD( inv( A_(i,i) ) );
      for( size_t j=0UL; j<N; ++j ) {
         X_(i,j) = B_(i,j) * invD;
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  FUNCTIONS FOR SOLVING LINEAR SYSTEMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Solving the given \f$ N \times N \f$ linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_matrix
//
// \param A The NxN dense system matrix.
// \param x The dense solution vector.
// \param b The N-dimensional dense right-hand side vector.
// \return void
// \exception std::invalid_argument Invalid non-square system matrix provided.
// \exception std::invalid_argument Invalid right-hand side vector provided.
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for the given linear system of equations \f$ A*x=b \f$, where
// \a A is the given system matrix, \a x is the solution vector, and \a b is the given right-hand
// side vector.

   \code
   blaze::DynamicMatrix<double> A;  // The square general system matrix
   blaze::DynamicVector<double> b;  // The right-hand side vector
   // ... Resizing and initialization

   blaze::DynamicVector<double> x;  // The solution vector
   solve( A, x, b );
   \endcode

// The \c solve() function will automatically select the most suited direct solver algorithm
// depending on the size and type of the given system matrix. For small matrices of up to 6x6,
// both functions use manually optimized kernels for maximum performance. For matrices larger
// than 6x6 the computation is performed by means of the most suited LAPACK solver method.
//
// In case the type of the matrix does not provide additional compile time information about
// its structure (symmetric, lower, upper, diagonal, ...), the information can be provided
// manually by means of declaration operations when calling the \c solve() function:

   \code
   blaze::DynamicMatrix<double> A;  // The square lower system matrix
   blaze::DynamicVector<double> b;  // The right-hand side vector
   // ... Resizing and initialization

   blaze::DynamicVector<double> x;  // The solution vector

   solve( declsym( A ), x, b );     // Solving the LSE with a symmetric system matrix
   solve( declherm( A ), x, b );    // Solving the LSE with an Hermitian system matrix
   solve( decllow( A ), x, b );     // Solving the LSE with a lower system matrix
   solve( declunilow( A ), x, b );  // Solving the LSE with an unilower system matrix
   solve( declupp( A ), x, b );     // Solving the LSE with an upper system matrix
   solve( decluniupp( A ), x, b );  // Solving the LSE with an uniupper system matrix
   solve( decldiag( A ), x, b );    // Solving the LSE with a diagonal system matrix
   \endcode

// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the size of the right-hand side vector doesn't match the dimensions of the system matrix;
//  - ... the given system matrix is singular.
//
// In all failure cases an exception is thrown.
//
// \note This function can only be used for dense matrices and vector with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices and vectors of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a x may already have been modified.
*/
template< typename MT   // Type of the system matrix
        , bool SO       // Storage order of the system matrix
        , typename VT1  // Type of the solution vector
        , bool TF1      // Transpose flag of the solution vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
void solve( const DenseMatrix<MT,SO>& A, DenseVector<VT1,TF1>& x, const DenseVector<VT2,TF2>& b )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<VT2> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT>, ElementType_t<VT2> );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square system matrix provided" );
   }
   else if( (*A).rows() != (*b).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side vector provided" );
   }

   switch( (*A).rows() ) {
      case 0UL:                         break;
      case 1UL: solve1x1( *A, *x, *b ); break;
      case 2UL: solve2x2( *A, *x, *b ); break;
      case 3UL: solve3x3( *A, *x, *b ); break;
      case 4UL: solve4x4( *A, *x, *b ); break;
      case 5UL: solve5x5( *A, *x, *b ); break;
      case 6UL: solve6x6( *A, *x, *b ); break;
      default : solveNxN( *A, *x, *b ); break;
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *x ), "Broken invariant detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Solving the given \f$ N \times N \f$ linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The NxN dense system matrix.
// \param X The dense solution matrix.
// \param B The N-dimensional dense right-hand side matrix.
// \return void
// \exception std::invalid_argument Invalid non-square system matrix provided.
// \exception std::invalid_argument Invalid right-hand side matrix provided.
// \exception std::runtime_error Solving LSE with singular system matrix failed.
//
// This function computes a solution for the given linear system of equations \f$ A*X=B \f$, where
// \a A is the given system matrix, the columns of \a X are the solution vectors, and the columns
// of \a B are the given right-hand side vectors.

   \code
   blaze::DynamicMatrix<double> A;  // The square general system matrix
   blaze::DynamicMatrix<double> B;  // The right-hand side matrix
   // ... Resizing and initialization

   blaze::DynamicMatrix<double> X;  // The solution matrix
   solve( A, X, B );
   \endcode

// The \c solve() function will automatically select the most suited direct solver algorithm
// depending on the size and type of the given system matrix. For small matrices of up to 6x6,
// both functions use manually optimized kernels for maximum performance. For matrices larger
// than 6x6 the computation is performed by means of the most suited LAPACK solver method.
//
// In case the type of the matrix does not provide additional compile time information about
// its structure (symmetric, lower, upper, diagonal, ...), the information can be provided
// manually by means of declaration operations when calling the \c solve() function:

   \code
   blaze::DynamicMatrix<double> A;  // The square lower system matrix
   blaze::DynamicMatrix<double> B;  // The right-hand side matrix
   // ... Resizing and initialization

   blaze::DynamicMatrix<double> X;  // The solution matrix

   solve( declsym( A ), X, B );     // Solving the LSE with a symmetric system matrix
   solve( declherm( A ), X, B );    // Solving the LSE with an Hermitian system matrix
   solve( decllow( A ), X, B );     // Solving the LSE with a lower system matrix
   solve( declunilow( A ), X, B );  // Solving the LSE with an unilower system matrix
   solve( declupp( A ), X, B );     // Solving the LSE with an upper system matrix
   solve( decluniupp( A ), X, B );  // Solving the LSE with an uniupper system matrix
   solve( decldiag( A ), X, B );    // Solving the LSE with a diagonal system matrix
   \endcode

// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the number of rows of the right-hand side matrix doesn't match the dimensions of the system matrix;
//  - ... the given system matrix is singular.
//
// In all failure cases an exception is thrown.
//
// \note The \c solve() function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a X may already have been modified.
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the solution matrix
        , bool SO2      // Storage order of the solution matrix
        , typename MT3  // Type of the right-hand side matrix
        , bool SO3 >    // Storage order of the right-hand side matrix
void solve( const DenseMatrix<MT1,SO1>& A, DenseMatrix<MT2,SO2>& X, const DenseMatrix<MT3,SO3>& B )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT1>  );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT3> );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT2> );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ElementType_t<MT1>, ElementType_t<MT3> );

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square system matrix provided" );
   }
   else if( (*A).rows() != (*B).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side matrix provided" );
   }

   switch( (*A).rows() ) {
      case 0UL: solve0x0( *A, *X, *B ); break;
      case 1UL: solve1x1( *A, *X, *B ); break;
      case 2UL: solve2x2( *A, *X, *B ); break;
      case 3UL: solve3x3( *A, *X, *B ); break;
      case 4UL: solve4x4( *A, *X, *B ); break;
      case 5UL: solve5x5( *A, *X, *B ); break;
      case 6UL: solve6x6( *A, *X, *B ); break;
      default : solveNxN( *A, *X, *B ); break;
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *X ), "Broken invariant detected" );
}
//*************************************************************************************************

} // namespace blaze

#endif
