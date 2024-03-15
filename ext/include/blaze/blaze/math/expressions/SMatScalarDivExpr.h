//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatScalarDivExpr.h
//  \brief Header file for the sparse matrix/scalar division expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATSCALARDIVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATSCALARDIVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsMultExpr.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/FloatingPoint.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATSCALARDIVEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-scalar divisions.
// \ingroup sparse_matrix_expression
//
// The SMatScalarMult class represents the compile time expression for divisions between
// a sparse matrix and a scalar value.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , typename ST  // Type of the right-hand side scalar value
        , bool SO >    // Storage order
class SMatScalarDivExpr
   : public MatScalarDivExpr< SparseMatrix< SMatScalarDivExpr<MT,ST,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;     //!< Result type of the sparse matrix expression.
   using RN = ReturnType_t<MT>;     //!< Return type of the sparse matrix expression.
   using CT = CompositeType_t<MT>;  //!< Composite type of the sparse matrix expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If the matrix operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = !IsTemporary_v<RN>;

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN>() / std::declval<ST>() );
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the division expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the division expression. In case the given sparse
       matrix expression of type \a MT requires an intermediate evaluation, \a useAssign will
       be set to 1 and the division expression will be evaluated via the \a assign function
       family. Otherwise Otherwise \a useAssign will be set to 0 and the expression will be
       evaluated via the function call operator. */
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
       In case either the target matrix or the sparse matrix operand is not SMP assignable and
       the matrix operand requires an intermediate evaluation, the variable is set to 1 and the
       expression specific evaluation strategy is selected. Otherwise the variable is set to 0
       and the default strategy is chosen. */
   template< typename MT2 >
   static constexpr bool UseSMPAssign_v =
      ( ( !MT2::smpAssignable || !MT::smpAssignable ) && useAssign );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatScalarDivExpr instance.
   using This = SMatScalarDivExpr<MT,ST,SO>;

   //! Base type of this SMatScalarDivExpr instance.
   using BaseType = MatScalarDivExpr< SparseMatrix<This,SO> >;

   using ResultType    = MultTrait_t<RT,ST>;           //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SMatScalarDivExpr& >;

   //! Composite data type of the sparse matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Composite type of the right-hand side scalar value.
   using RightOperand = ST;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the sparse matrix/scalar division expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse matrix expression.
      using Element = ValueIndexPair<ElementType>;

      //! Iterator type of the sparse matrix expression.
      using IteratorType = ConstIterator_t< RemoveReference_t<LeftOperand> >;

      using IteratorCategory = std::forward_iterator_tag;  //!< The iterator category.
      using ValueType        = Element;                    //!< Type of the underlying pointers.
      using PointerType      = ValueType*;                 //!< Pointer return type.
      using ReferenceType    = ValueType&;                 //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                  //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying pointers.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      */
      inline ConstIterator( IteratorType matrix, RightOperand scalar )
         : matrix_( matrix )  // Iterator over the elements of the left-hand side sparse matrix expression
         , scalar_( scalar )  // Right-hand side scalar of the division expression
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented expression iterator.
      */
      inline ConstIterator& operator++() {
         ++matrix_;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The element at the current iterator position.
      */
      inline const Element operator*() const {
         return Element( matrix_->value() / scalar_, matrix_->index() );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return Reference to the sparse matrix element at the current iterator position.
      */
      inline const ConstIterator* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse element.
      //
      // \return The current value of the sparse element.
      */
      inline ReturnType value() const {
         return matrix_->value() / scalar_;
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return matrix_->index();
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return matrix_ == rhs.matrix_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return matrix_ != rhs.matrix_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two expression iterators.
      //
      // \param rhs The right-hand side expression iterator.
      // \return The number of elements between the two expression iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return matrix_ - rhs.matrix_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType matrix_;  //!< Iterator over the elements of the left-hand side sparse matrix expression.
      RightOperand scalar_;  //!< Right-hand side scalar of the division expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatScalarDivExpr class.
   //
   // \param matrix The left-hand side sparse matrix of the division expression.
   // \param scalar The right-hand side scalar of the division expression.
   */
   inline SMatScalarDivExpr( const MT& matrix, ST scalar ) noexcept
      : matrix_( matrix )  // Left-hand side sparse matrix of the division expression
      , scalar_( scalar )  // Right-hand side scalar of the division expression
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
      BLAZE_INTERNAL_ASSERT( i < matrix_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < matrix_.columns(), "Invalid column access index" );
      return matrix_(i,j) / scalar_;
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
      if( i >= matrix_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= matrix_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( matrix_.begin(i), scalar_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( matrix_.end(i), scalar_ );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return matrix_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return matrix_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return matrix_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return matrix_.nonZeros(i);
   }
   //**********************************************************************************************

   //**Find function*******************************************************************************
   /*!\brief Searches for a specific matrix element.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   */
   inline ConstIterator find( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( matrix_.find( i, j ), scalar_ );
   }
   //**********************************************************************************************

   //**LowerBound function*************************************************************************
   /*!\brief Returns an iterator to the first index not less then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index not less then the given index, end() iterator otherwise.
   */
   inline ConstIterator lowerBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( matrix_.lowerBound( i, j ), scalar_ );
   }
   //**********************************************************************************************

   //**UpperBound function*************************************************************************
   /*!\brief Returns an iterator to the first index greater then the given index.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the first index greater then the given index, end() iterator otherwise.
   */
   inline ConstIterator upperBound( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( matrix_.upperBound( i, j ), scalar_ );
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse matrix operand.
   //
   // \return The left-hand side sparse matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return matrix_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side scalar operand.
   //
   // \return The right-hand side scalar operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return scalar_;
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
      return matrix_.canAlias( alias );
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
      return matrix_.isAliased( alias );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  matrix_;  //!< Left-hand side sparse matrix of the division expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the division expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-scalar
   // division expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.matrix_ );
      (*lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-scalar division to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side division expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-scalar
   // division expression to a sparse matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.matrix_ );
      (*lhs) /= rhs.scalar_;
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
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
   /*!\brief Subtraction assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
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
   /*!\brief Schur product assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case
   // the operand requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
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

   //**SMP assignment to dense matrices************************************************************
   // No special implementation for the SMP assignment to dense matrices.
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   // No special implementation for the SMP assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpAddAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSubAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix-scalar division to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side division expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix-scalar division expression to a dense matrix. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatScalarDivExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSchurAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP Schur product assignment to sparse matrices*********************************************
   // No special implementation for the SMP Schur product assignment to sparse matrices.
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_NOT_BE_FLOATING_POINT_TYPE( ElementType );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the sparse matrix/scalar division operator.
// \ingroup sparse_matrix
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , typename ST >  // Type of the right-hand side scalar
using SMatScalarDivExprHelper_t =
   If_t< IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > ||
         IsFloatingPoint_v< UnderlyingBuiltin_t<ST> >
       , If_t< IsBuiltin_v<ST>
             , DivTrait_t< UnderlyingBuiltin_t<MT>, ST >
             , decltype( inv( std::declval<ST>() ) ) >
       , ST >;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the division between a sparse matrix and a scalar value
//        (\f$ A=B/s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side sparse matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This function implements the default treatment of the sparse matrix/scalar division.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , bool SO      // Storage order of the left-hand side sparse matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< !IsZero_v<MT> &&
                      !IsInvertible_v< SMatScalarDivExprHelper_t<MT,ST> > >* = nullptr >
inline decltype(auto) smatscalardiv( const SparseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = SMatScalarDivExprHelper_t<MT,ST>;
   using ReturnType = const SMatScalarDivExpr<MT,ScalarType,SO>;

   return ReturnType( *mat, scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the division between a sparse matrix and a scalar value
//        (\f$ A=B/s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side sparse matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This function implements a performance optimized treatment of the sparse matrix/scalar division.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , bool SO      // Storage order of the left-hand side sparse matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< !IsZero_v<MT> &&
                      IsInvertible_v< SMatScalarDivExprHelper_t<MT,ST> > >* = nullptr >
inline decltype(auto) smatscalardiv( const SparseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = SMatScalarDivExprHelper_t<MT,ST>;
   using ReturnType = const SMatScalarMultExpr<MT,ScalarType,SO>;

   return ReturnType( *mat, inv(scalar) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the division between a zero matrix and a scalar value
//        (\f$ A=B/s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side zero matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the division between a zero
// matrix and a scalar value. It returns a zero matrix.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , bool SO      // Storage order of the left-hand side zero matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsZero_v<MT> >* = nullptr >
inline decltype(auto)
   smatscalardiv( const SparseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( scalar );

   using ReturnType = const DivTrait_t< ResultType_t<MT>, ST >;

   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ReturnType, SO );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*mat).rows(), (*mat).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division operator for the division of a sparse matrix by a scalar value (\f$ A=B/s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side sparse matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This operator represents the division of a sparse matrix by a scalar value:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = A / 0.24;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the involved data types \a MT::ElementType and \a ST. Note that this operator only
// works for scalar values of built-in data type.
//
// \note A division by zero is only checked by a user assert.
*/
template< typename MT  // Type of the left-hand side sparse matrix
        , bool SO      // Storage order of the left-hand side sparse matrix
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator/( const SparseMatrix<MT,SO>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != ST{}, "Division by zero detected" );

   return smatscalardiv( *mat, scalar );
}
//*************************************************************************************************

} // namespace blaze

#endif
