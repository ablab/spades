//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDVecSolveExpr.h
//  \brief Header file for the dense matrix/dense vector solver expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDVECSOLVEEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDVECSOLVEEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/MatVecSolveExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatVecSolveExpr.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATDVECSOLVEEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-dense vector solvers.
// \ingroup dense_vector_expression
//
// The DMatDVecSolveExpr class represents the compile time expression for a dense linear system
// of equations (LSE) with a single right-hand side vector.
*/
template< typename MT  // Type of the dense system matrix
        , typename VT  // Type of the dense right-hand side vector
        , bool TF >    // Transpose flag
class DMatDVecSolveExpr
   : public MatVecSolveExpr< DenseVector< DMatDVecSolveExpr<MT,VT,TF>, TF > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using MRT = ResultType_t<MT>;  //!< Result type of the left-hand side dense matrix expression.
   using VRT = ResultType_t<VT>;  //!< Result type of the right-hand side dense vector expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatDVecSolveExpr instance.
   using This = DMatDVecSolveExpr<MT,VT,TF>;

   //! Base type of this DMatDVecSolveExpr instance.
   using BaseType = MatVecSolveExpr< DenseVector<This,TF> >;

   using ResultType    = SolveTrait_t<MRT,VRT>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<ResultType>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Composite type of the right-hand side dense vector expression.
   using RightOperand = If_t< IsExpression_v<VT>, const VT, const VT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatDVecSolveExpr class.
   //
   // \param mat The left-hand side matrix operand of the solver expression.
   // \param vec The right-hand side vector operand of the solver expression.
   */
   inline DMatDVecSolveExpr( const MT& mat, const VT& vec ) noexcept
      : mat_( mat )  // Left-hand side dense matrix of the solver expression
      , vec_( vec )  // Right-hand side dense vector of the solver expression
   {}
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return mat_.rows();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return mat_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return vec_;
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side dense matrix of the solver expression.
   RightOperand vec_;  //!< Right-hand side dense vector of the solver expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix/dense vector solver expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side solver expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix/dense
   // vector solver expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      solve( rhs.leftOperand(), *lhs, rhs.rightOperand() );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix/dense vector solver expression to a sparse vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side solver expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix/dense
   // vector solver expression to a sparse vector.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix/dense vector solver expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side solver expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix/dense vector solver expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      addAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix/dense vector solver expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side solver expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix/dense vector solver expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      subAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense matrix/dense vector solver expression to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side solver expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // matrix/dense vector solver expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a dense matrix/dense vector solver expression to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side solver expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense
   // matrix/dense vector solver expression to a dense vector.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT2,TF>& lhs, const DMatDVecSolveExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType, TF );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      divAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATVECSOLVEEXPR( MT, VT );
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
/*!\brief Solving the given \f$ N \times N \f$ linear system of equations (\f$ A*x=b \f$).
// \ingroup dense_vector
//
// \param A The NxN dense system matrix.
// \param b The N-dimensional dense right-hand side vector.
// \return void
// \exception std::invalid_argument Invalid non-square system matrix provided.
// \exception std::invalid_argument Invalid right-hand side vector provided.
//
// This function returns an expression representing the solution of the given dense linear system
// of equations (LSE):

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::DynamicVector<double> b, x;
   // ... Resizing and initialization
   x = solve( A, b );
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

   x = solve( declsym( A ), b );     // Solving the LSE with a symmetric system matrix
   x = solve( declherm( A ), b );    // Solving the LSE with an Hermitian system matrix
   x = solve( decllow( A ), b );     // Solving the LSE with a lower system matrix
   x = solve( declunilow( A ), b );  // Solving the LSE with an unilower system matrix
   x = solve( declupp( A ), b );     // Solving the LSE with an upper system matrix
   x = solve( decluniupp( A ), b );  // Solving the LSE with an uniupper system matrix
   x = solve( decldiag( A ), b );    // Solving the LSE with a diagonal system matrix
   \endcode

// The function fails if ...
//
//  - ... the given system matrix is not a square matrix;
//  - ... the size of the right-hand side vector doesn't match the dimensions of the system matrix.
//
// In all failure cases an exception is thrown.
//
// \note The solve() function can only be used for dense matrices and vectors with \c float,
// \c double, \c complex<float> or \c complex<double> element type. The attempt to call the
// function with matrices and vectors of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note It is not possible to use any kind of view on the expression object returned by the
// \c solve() function. Also, it is not possible to access individual elements via the function
// call operator on the expression object:

   \code
   row( solve( A, b ), 2UL );  // Compilation error: Views cannot be used on an solve() expression!
   solve( A, b )[2];           // Compilation error: It is not possible to access individual elements!
   \endcode
*/
template< typename MT  // Type of the system matrix
        , bool SO      // Storage order of the system matrix
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline decltype(auto) solve( const DenseMatrix<MT,SO>& A, const DenseVector<VT,TF>& b )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *A ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square system matrix provided" );
   }
   else if( (*A).rows() != (*b).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-hand side vector provided" );
   }

   using ReturnType = const DMatDVecSolveExpr<MT,VT,TF>;
   return ReturnType( *A, *b );
}
//*************************************************************************************************

} // namespace blaze

#endif
