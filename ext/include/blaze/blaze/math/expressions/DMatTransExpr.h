//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatTransExpr.h
//  \brief Header file for the dense matrix transpose expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATTRANSEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATTRANSEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DMatTransposer.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/Transformation.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/system/Inline.h>
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
//  CLASS DMATTRANSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix transpositions.
// \ingroup dense_matrix_expression
//
// The DMatTransExpr class represents the compile time expression for transpositions of
// dense matrices.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
class DMatTransExpr
   : public MatTransExpr< DenseMatrix< DMatTransExpr<MT,SO>, SO > >
   , private If_t< IsComputation_v<MT>, Computation, Transformation >
{
 private:
   //**Type definitions****************************************************************************
   //! Definition of the GetConstIterator type trait.
   BLAZE_CREATE_GET_TYPE_MEMBER_TYPE_TRAIT( GetConstIterator, ConstIterator, INVALID_TYPE );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the transposition expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the transposition expression. In case the given dense
       matrix expression of type \a MT requires an intermediate evaluation, \a useAssign will
       be set to 1 and the transposition expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to 0 and the expression will be evaluated
       via the subscript operator. */
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
       In case the target matrix is SMP assignable and the dense matrix operand requires an
       intermediate evaluation, the variable is set to 1 and the expression specific evaluation
       strategy is selected. Otherwise the variable is set to 0 and the default strategy is
       chosen. */
   template< typename MT2 >
   static constexpr bool UseSMPAssign_v = ( MT2::smpAssignable && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatTransExpr instance.
   using This = DMatTransExpr<MT,SO>;

   //! Base type of this DMatTransExpr instance.
   using BaseType = MatTransExpr< DenseMatrix<This,SO> >;

   using ResultType    = TransposeType_t<MT>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;  //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = ResultType_t<MT>;            //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;           //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;            //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DMatTransExpr& >;

   //! Iterator over the elements of the dense matrix.
   using ConstIterator = GetConstIterator_t<MT>;

   //! Composite data type of the dense matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatTransExpr class.
   //
   // \param dm The dense matrix operand of the transposition expression.
   */
   explicit inline DMatTransExpr( const MT& dm ) noexcept
      : dm_( dm )  // Dense matrix of the transposition expression
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
      BLAZE_INTERNAL_ASSERT( i < dm_.columns(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.rows()   , "Invalid column access index" );
      return dm_(j,i);
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
      if( i >= dm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= dm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Load function*******************************************************************************
   /*!\brief Access to the SIMD elements of the matrix.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return Reference to the accessed values.
   */
   BLAZE_ALWAYS_INLINE auto load( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < dm_.columns(), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < dm_.rows()   , "Invalid column access index" );
      BLAZE_INTERNAL_ASSERT( !SO || ( i % SIMDSIZE == 0UL ), "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( SO  || ( j % SIMDSIZE == 0UL ), "Invalid column access index" );
      return dm_.load(j,i);
   }
   //**********************************************************************************************

   //**Low-level data access***********************************************************************
   /*!\brief Low-level data access to the matrix elements.
   //
   // \return Pointer to the internal element storage.
   */
   inline const ElementType* data() const noexcept {
      return dm_.data();
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( dm_.begin(i) );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( dm_.end(i) );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return dm_.columns();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.rows();
   }
   //**********************************************************************************************

   //**Spacing function****************************************************************************
   /*!\brief Returns the spacing between the beginning of two rows/columns.
   //
   // \return The spacing between the beginning of two rows/columns.
   */
   inline size_t spacing() const noexcept {
      return dm_.spacing();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the dense matrix.
   //
   // \return The number of non-zero elements in the dense matrix.
   */
   inline size_t nonZeros() const {
      return dm_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row/column.
   //
   // \param i The index of the row/column.
   // \return The number of non-zero elements of row/column \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return dm_.nonZeros( i );
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense matrix operand.
   //
   // \return The dense matrix operand.
   */
   inline Operand operand() const noexcept {
      return dm_;
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
      return dm_.isAliased( alias );
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
      return dm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dm_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return dm_.canSMPAssign();
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand dm_;  //!< Dense matrix of the transposition expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix
   // transposition expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      assign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix transposition expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix
   // transposition expression to a sparse matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OppositeType, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
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
   /*!\brief Addition assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      addAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      subAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      schurAssign( tmp, rhs.dm_ );
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

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix
   // transposition expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      smpAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix transposition expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side transposition expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix
   // transposition expression to a sparse matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO == SO2, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( OppositeType, !SO );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT2, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const TmpType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      smpAddAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix transposition expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      smpSubAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a dense matrix transposition expression to a dense
   //        matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side transposition expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a dense
   // matrix transposition expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatTransExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      DMatTransposer<MT2,!SO2> tmp( *lhs );
      smpSchurAssign( tmp, rhs.dm_ );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   // No special implementation for the SMP Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, !SO );
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
/*!\brief Calculation of the transpose of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be transposed.
// \return The transpose of the matrix.
//
// This function returns an expression representing the transpose of the given dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::DynamicMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   B = trans( A );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) trans( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatTransExpr<MT,!SO>;
   return ReturnType( *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conditional calculation of the transpose of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be transposed.
// \return The transpose of the matrix.
//
// This function does not calculate the transpose of the given dense matrix and returns a
// reference to the matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) transIf( FalseType, const DenseMatrix<MT,SO>& dm )
{
   return *dm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conditional calculation of the transpose of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be transposed.
// \return The transpose of the matrix.
//
// This function calculates the transpose of the given dense matrix and returns an expression
// representing the transpose of the given dense matrix.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) transIf( TrueType, const DenseMatrix<MT,SO>& dm )
{
   return trans( *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conditional calculation of the transpose of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be transposed.
// \return The transpose of the matrix.
//
// In case the given compile time condition evaluates to \a true, this function returns an
// expression representing the transpose of the given dense matrix. Otherwise the function
// returns a reference to the given matrix.
*/
template< bool B       // Compile time condition
        , typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) transIf( const DenseMatrix<MT,SO>& dm )
{
   return transIf( BoolConstant<B>(), *dm );
}
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct HasConstDataAccess< DMatTransExpr<MT,SO> >
   : public HasConstDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsAligned< DMatTransExpr<MT,SO> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO >
struct IsPadded< DMatTransExpr<MT,SO> >
   : public IsPadded<MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
