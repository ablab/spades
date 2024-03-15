//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatEvalExpr.h
//  \brief Header file for the sparse matrix evaluation expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATEVALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATEVALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATEVALEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the forced evaluation of sparse matrices.
// \ingroup sparse_matrix_expression
//
// The SMatEvalExpr class represents the compile time expression for the forced evaluation
// of a sparse matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
class SMatEvalExpr
   : public MatEvalExpr< SparseMatrix< SMatEvalExpr<MT,SO>, SO > >
   , private Computation
{
 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatEvalExpr instance.
   using This = SMatEvalExpr<MT,SO>;

   //! Base type of this SMatEvalExpr instance.
   using BaseType = MatEvalExpr< SparseMatrix<This,SO> >;

   using ResultType    = ResultType_t<MT>;     //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<MT>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<MT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite data type of the sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatEvalExpr class.
   //
   // \param sm The sparse matrix operand of the evaluation expression.
   */
   explicit inline SMatEvalExpr( const MT& sm ) noexcept
      : sm_( sm )  // Sparse matrix of the evaluation expression
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const {
      BLAZE_INTERNAL_ASSERT( i < sm_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < sm_.columns(), "Invalid column access index" );
      return sm_(i,j);
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid matrix access index.
   */
   inline ReturnType at( size_t i, size_t j ) const {
      if( i >= sm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= sm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return sm_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return sm_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return sm_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return sm_.nonZeros(i);
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the sparse matrix operand.
   //
   // \return The sparse matrix operand.
   */
   inline Operand operand() const noexcept {
      return sm_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return sm_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return sm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return sm_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sm_;  //!< Sparse matrix of the evaluation expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix evaluation
   // expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix evaluation
   // expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix
   // evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix
   // evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void addAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix
   // evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void subAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      schurAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Schur product assignment to sparse matrices*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void schurAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      schurAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void multAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      multAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void multAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      multAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix
   // evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void smpAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix
   // evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void smpAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void smpAddAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix evaluation expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix evaluation expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void smpSubAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix evaluation expression to a dense
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSchurAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix evaluation expression to a sparse
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void smpSchurAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSchurAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a sparse matrix evaluation expression to a dense
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a sparse
   // matrix evaluation expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void smpMultAssign( DenseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpMultAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a sparse matrix evaluation expression to a sparse
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side evaluation expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a sparse
   // matrix evaluation expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void smpMultAssign( SparseMatrix<MT2,SO2>& lhs, const SMatEvalExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpMultAssign( *lhs, rhs.sm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Forces the evaluation of the given sparse matrix expression \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The evaluated sparse matrix.
//
// The \a eval function forces the evaluation of the given sparse matrix expression \a sm.
// The function returns an expression representing the operation.\n
// The following example demonstrates the use of the \a eval function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = eval( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) eval( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const SMatEvalExpr<MT,SO>;
   return ReturnType( *sm );
}
//*************************************************************************************************

} // namespace blaze

#endif
