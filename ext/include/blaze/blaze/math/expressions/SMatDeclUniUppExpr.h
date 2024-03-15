//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatDeclUniUppExpr.h
//  \brief Header file for the sparse matrix uniupper declaration expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATDECLUNIUPPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATDECLUNIUPPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/adaptors/uniuppermatrix/BaseTemplate.h>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/UniUpper.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Declaration.h>
#include <blaze/math/expressions/DeclUniUppExpr.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/traits/DeclUniUppTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/GetMemberType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATDECLUNIUPPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the explicit uniupper declaration of sparse matrices.
// \ingroup sparse_matrix_expression
//
// The SMatDeclUniUppExpr class represents the compile time expression for the explicit uniupper
// declaration of a sparse matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
class SMatDeclUniUppExpr
   : public DeclUniUppExpr< SparseMatrix< SMatDeclUniUppExpr<MT,SO>, SO > >
   , public Declaration<MT>
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;  //!< Result type of the sparse matrix expression.

   //! Definition of the GetConstIterator type trait.
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetConstIterator, ConstIterator, INVALID_TYPE );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the uniupper declaration expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the uniupper declaration expression. In case the given
       sparse matrix expression of type \a MT requires an intermediate evaluation, \a useAssign
       will be set to 1 and the uniupper declaration expression will be evaluated via the
       \a assign function family. Otherwise \a useAssign will be set to 0 and the expression
       will be evaluated via the subscript operator. */
   static constexpr bool useAssign = RequiresEvaluation_v<MT>;

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**Parallel evaluation strategy****************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case the target matrix is SMP assignable and the sparse matrix operand requires an
       intermediate evaluation, the variable is set to 1 and the expression specific evaluation
       strategy is selected. Otherwise the variable is set to 0 and the default strategy is
       chosen. */
   template< typename MT2 >
   static constexpr bool UseSMPAssign_v = ( MT2::smpAssignable && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatDeclUniUppExpr instance.
   using This = SMatDeclUniUppExpr<MT,SO>;

   //! Base type of this SMatDeclUniUppExpr instance.
   using BaseType = DeclUniUppExpr< SparseMatrix<This,SO> >;

   using ResultType    = DeclUniUppTrait_t<RT>;        //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< RequiresEvaluation_v<MT>, const ResultType, const SMatDeclUniUppExpr& >;

   //! Iterator over the elements of the dense matrix.
   using ConstIterator = GetConstIterator_t<MT>;

   //! Composite data type of the sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatDeclUniUppExpr class.
   //
   // \param sm The sparse matrix operand of the decluniupp expression.
   */
   explicit inline SMatDeclUniUppExpr( const MT& sm ) noexcept
      : sm_( sm )  // Sparse matrix of the decluniupp expression
   {
      BLAZE_INTERNAL_ASSERT( isSquare( *sm ), "Non-square matrix detected" );
   }
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

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator to the first non-zero element of row \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( sm_.begin(i) );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator just past the last non-zero element of row \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( sm_.end(i) );
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
   Operand sm_;  //!< Sparse matrix of the decluniupp expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix decluniupp
   // expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix decluniupp
   // expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Addition assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix
   // decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Addition assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix
   // decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto addAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Subtraction assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Subtraction assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto subAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Schur product assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Schur product assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto schurAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Multiplication assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto multAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief Multiplication assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto multAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
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
   /*!\brief SMP assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix
   // decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix
   // decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP addition assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP addition assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpAddAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP subtraction assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP subtraction assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpSubAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP Schur product assignment of a sparse matrix decluniupp expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP Schur product assignment of a sparse matrix decluniupp expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpSchurAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP multiplication assignment of a sparse matrix decluniupp expression to a dense
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side decluniupp expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // sparse matrix decluniupp expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpMultAssign( DenseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   /*!\brief SMP multiplication assignment of a sparse matrix decluniupp expression to a sparse
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side decluniupp expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // sparse matrix decluniupp expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpMultAssign( SparseMatrix<MT2,SO2>& lhs, const SMatDeclUniUppExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
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
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE( MT );
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
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given sparse matrix expression \a sm as uniupper.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The redeclared sparse matrix.
//
// This function declares the given sparse matrix expression \a sm as uniupper. The function
// returns an expression representing the operation.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , DisableIf_t< IsUniUpper_v<MT> || IsLower_v<MT> >* = nullptr >
inline const SMatDeclUniUppExpr<MT,SO> decluniupp_backend( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isSquare( *sm ), "Non-square matrix detected" );

   return SMatDeclUniUppExpr<MT,SO>( *sm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given lower sparse matrix expression \a sm as uniupper.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The redeclared sparse matrix.
//
// This function declares the given lower sparse matrix expression \a sm as uniupper. The
// function returns an identity matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , EnableIf_t< !IsUniUpper_v<MT> && IsLower_v<MT> >* = nullptr >
inline const IdentityMatrix<ElementType_t<MT>,SO> decluniupp_backend( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isSquare( *sm ), "Non-square matrix detected" );

   return IdentityMatrix<ElementType_t<MT>,SO>( (*sm).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Redeclares the given uniupper sparse matrix expression \a sm as uniupper.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The redeclared sparse matrix.
//
// This function redeclares the given uniupper sparse matrix expression \a sm as uniupper. The
// function returns a reference to the already uniupper matrix expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , EnableIf_t< IsUniUpper_v<MT> >* = nullptr >
inline const MT& decluniupp_backend( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isSquare( *sm ), "Non-square matrix detected" );

   return *sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Declares the given sparse matrix expression \a sm as uniupper.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The redeclared sparse matrix.
// \exception std::invalid_argument Invalid uniupper matrix specification.
//
// The \a decluniupp function declares the given sparse matrix expression \a sm as uniupper. In
// case the given matrix is not a square matrix, a \a std::invalid_argument exception is thrown.\n
// The following example demonstrates the use of the \a decluniupp function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = decluniupp( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) decluniupp( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );

   if( !isSquare( *sm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uniupper matrix specification" );
   }

   return decluniupp_backend( *sm );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsUniUpper< SMatDeclUniUppExpr<MT,SO> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
