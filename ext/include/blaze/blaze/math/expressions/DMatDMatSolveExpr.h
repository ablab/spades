//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDMatSolveExpr.h
//  \brief Header file for the dense matrix/dense matrix solver expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDMATSOLVEEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDMATSOLVEEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatSolveExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatSolveExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATDMATSOLVEEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-dense matrix solvers.
// \ingroup dense_matrix_expression
//
// The DMatDMatSolveExpr class represents the compile time expression for a dense linear system
// of equations (LSE) with multiple right-hand side vectors.
*/
template< typename MT1  // Type of the dense system matrix
        , typename MT2  // Type of the dense right-hand side matrix
        , bool SO >     // Storage order
class DMatDMatSolveExpr
   : public MatMatSolveExpr< DenseMatrix< DMatDMatSolveExpr<MT1,MT2,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;  //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;  //!< Result type of the right-hand side dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatDMatSolveExpr instance.
   using This = DMatDMatSolveExpr<MT1,MT2,SO>;

   //! Base type of this DMatDMatSolveExpr instance.
   using BaseType = MatMatSolveExpr< DenseMatrix<This,SO> >;

   using ResultType    = SolveTrait_t<RT1,RT2>;        //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<ResultType>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side dense matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatDMatSolveExpr class.
   //
   // \param lhs The left-hand side matrix operand of the solver expression.
   // \param rhs The right-hand side matrix operand of the solver expression.
   */
   inline DMatDMatSolveExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the solver expression
      , rhs_( rhs )  // Right-hand side dense matrix of the solver expression
   {}
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return rhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return rhs_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense matrix operand.
   //
   // \return The right-hand side dense matrix operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an aliasing effect is possible, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the given alias is contained in this expression, \a false if not.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the solver expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the solver expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix/dense matrix solver expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side solver expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix/dense
   // matrix solver expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO2>& lhs, const DMatDMatSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      solve( rhs.leftOperand(), *lhs, rhs.rightOperand() );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix/dense matrix solver expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side solver expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix/dense
   // matrix solver expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO2 >   // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO2>& lhs, const DMatDMatSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OppositeType, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix/dense matrix solver expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side solver expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix/dense matrix solver expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      addAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix/dense matrix solver expression to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side solver expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix/dense matrix solver expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      subAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix-dense matrix solver expression to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side solver expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix/dense matrix solver expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO2>& lhs, const DMatDMatSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to sparse matrices*************************************************
   // No special implementation for the Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT2, SO );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATSOLVEEXPR( MT1, MT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Solving the given \f$ N \times N \f$ linear system of equations (\f$ A*X=B \f$).
// \ingroup dense_matrix
//
// \param A The NxN dense system matrix.
// \param B The NxRHS-dimensional dense right-hand side matrix.
// \return void
// \exception std::invalid_argument Invalid non-square system matrix provided.
// \exception std::invalid_argument Invalid right-hand side matrix provided.
//
// This function returns an expression representing the solution of the given dense linear system
// of equations (LSE):

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::DynamicMatrix<double> B, X;
   // ... Resizing and initialization
   X = solve( A, B );
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

   X = solve( declsym( A ), B );     // Solving the LSE with a symmetric system matrix
   X = solve( declherm( A ), B );    // Solving the LSE with an Hermitian system matrix
   X = solve( decllow( A ), B );     // Solving the LSE with a lower system matrix
   X = solve( declunilow( A ), B );  // Solving the LSE with an unilower system matrix
   X = solve( declupp( A ), B );     // Solving the LSE with an upper system matrix
   X = solve( decluniupp( A ), B );  // Solving the LSE with an uniupper system matrix
   X = solve( decldiag( A ), B );    // Solving the LSE with a diagonal system matrix
   \endcode

// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the number of rows of the right-hand side matrix doesn't match the dimensions of the system matrix.
//
// In all failure cases an exception is thrown.

// \note The solve() function can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note It is not possible to use any kind of view on the expression object returned by the
// \c solve() function. Also, it is not possible to access individual elements via the function
// call operator on the expression object:

   \code
   rows( solve( A, B ), { 2UL, 4UL } );  // Compilation error: Views cannot be used on an solve() expression!
   solve( A, B )(1,2);                   // Compilation error: It is not possible to access individual elements!
   \endcode
*/
template< typename MT1  // Type of the system matrix
        , bool SO1      // Storage order of the system matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline decltype(auto) solve( const DenseMatrix<MT1,SO1>& A, const DenseMatrix<MT2,SO2>& B )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square system matrix provided" );
   }
   else if( (*A).rows() != (*B).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side matrix provided" );
   }

   using ReturnType = const DMatDMatSolveExpr<MT1,MT2,SO2>;
   return ReturnType( *A, *B );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense matrix/dense matrix solver expression.
// \ingroup dense_matrix
//
// \param dm The given dense matrix/dense matrix solver expression.
// \return The determinant of the given expression.
//
// This function computes the determinant of the given dense matrix/dense matrix solver expression.
//
// \note The computation of the determinant is numerically unreliable since especially for large
// matrices the value can overflow during the computation. Please note that this function does
// not guarantee that it is possible to compute the determinant with the given matrix!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT1  // Type of the dense system matrix
        , typename MT2  // Type of the dense right-hand side matrix
        , bool SO >     // Storage order
inline decltype(auto) det( const DMatDMatSolveExpr<MT1,MT2,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return det( evaluate( *dm ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
