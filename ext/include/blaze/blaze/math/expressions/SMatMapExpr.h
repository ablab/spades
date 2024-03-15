//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatMapExpr.h
//  \brief Header file for the sparse matrix map expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATMAPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATMAPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Functors.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATMAPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the sparse matrix map() function.
// \ingroup sparse_matrix_expression
//
// The SMatMapExpr class represents the compile time expression for the evaluation of a custom
// operation on each element of a sparse matrix via the map() function.
*/
template< typename MT  // Type of the sparse matrix
        , typename OP  // Type of the custom operation
        , bool SO >    // Storage order
class SMatMapExpr
   : public MatMapExpr< SparseMatrix< SMatMapExpr<MT,OP,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;    //!< Result type of the sparse matrix expression.
   using OT = OppositeType_t<MT>;  //!< Opposite type of the sparse matrix expression.
   using RN = ReturnType_t<MT>;    //!< Return type of the sparse matrix expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the map expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the map expression. In case the given sparse matrix
       expression of type \a MT requires an intermediate evaluation, \a useAssign will be
       set to 1 and the map expression will be evaluated via the \a assign function family.
       Otherwise \a useAssign will be set to 0 and the expression will be evaluated via the
       subscript operator. */
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
       In case either the target matrix or the sparse matrix operand is not SMP assignable or
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
   //! Type of this SMatMapExpr instance.
   using This = SMatMapExpr<MT,OP,SO>;

   //! Base type of this SMatMapExpr instance.
   using BaseType = MatMapExpr< SparseMatrix<This,SO> >;

   using ResultType    = MapTrait_t<RT,OP>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<RN>() ) );

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SMatMapExpr& >;

   //! Composite data type of the sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Data type of the custom unary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the sparse matrix map expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse matrix expression.
      using Element = ValueIndexPair<ElementType>;

      //! Iterator type of the sparse matrix expression.
      using IteratorType = ConstIterator_t< RemoveReference_t<Operand> >;

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
      //
      // \param it Iterator to the initial matrix element.
      // \param op The custom unary operation.
      */
      inline ConstIterator( IteratorType it, OP op )
         : it_( it )             // Iterator over the elements of the sparse matrix expression
         , op_( std::move(op) )  // The custom unary operation
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented expression iterator.
      */
      inline ConstIterator& operator++() {
         ++it_;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline const Element operator*() const {
         return Element( op_( it_->value() ), it_->index() );
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
         return op_( it_->value() );
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return it_->index();
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return it_ == rhs.it_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return it_ != rhs.it_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two expression iterators.
      //
      // \param rhs The right-hand side expression iterator.
      // \return The number of elements between the two expression iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return it_ - rhs.it_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType it_;  //!< Iterator over the elements of the sparse matrix expression.
      OP           op_;  //!< The custom unary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatMapExpr class.
   //
   // \param sm The sparse matrix operand of the map expression.
   // \param op The custom unary operation.
   */
   inline SMatMapExpr( const MT& sm, OP op ) noexcept
      : sm_( sm )             // Sparse matrix of the map expression
      , op_( std::move(op) )  // The custom unary operation
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
      return op_( sm_(i,j) );
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
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( sm_.begin(i), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( sm_.end(i), op_ );
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

   //**Find function*******************************************************************************
   /*!\brief Searches for a specific matrix element.
   //
   // \param i The row index of the search element.
   // \param j The column index of the search element.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   */
   inline ConstIterator find( size_t i, size_t j ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );
      return ConstIterator( sm_.find( i, j ), op_ );
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
      return ConstIterator( sm_.lowerBound( i, j ), op_ );
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
      return ConstIterator( sm_.upperBound( i, j ), op_ );
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

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the custom operation.
   //
   // \return A copy of the custom operation.
   */
   inline Operation operation() const {
      return op_;
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
   Operand   sm_;  //!< Sparse matrix of the map expression.
   Operation op_;  //!< The custom unary operation.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix map expression to a row-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix map
   // expression to a row-major sparse matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand requires
   // an intermediate evaluation and the underlying numeric data type of the operand and the
   // target matrix are identical.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,false>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> &&
                     IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.sm_ );

      const size_t m( rhs.rows() );

      for( size_t i=0UL; i<m; ++i ) {
         const auto end( (*lhs).end(i) );
         for( auto element=(*lhs).begin(i); element!=end; ++element ) {
            element->value() = rhs.op_( element->value() );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix map expression to a column-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix map
   // expression to a column-major sparse matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand requires
   // an intermediate evaluation and the underlying numeric data type of the operand and the
   // target matrix are identical.
   */
   template< typename MT2 >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,true>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> &&
                     IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, rhs.sm_ );

      const size_t n( rhs.columns() );

      for( size_t j=0UL; j<n; ++j ) {
         const auto end( (*lhs).end(j) );
         for( auto element=(*lhs).begin(j); element!=end; ++element ) {
            element->value() = rhs.op_( element->value() );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix map expression to a sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix map
   // expression to a sparse matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation and the underlying numeric data type of the operand and the
   // target vector differ.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> &&
                     !IsSame_v< UnderlyingScalar_t<MT>, UnderlyingScalar_t<MT2> > >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      (*lhs).reserve( tmp.nonZeros() );
      assign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix
   // map expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the operand requires an
   // intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      addAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto subAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      subAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the operand
   // requires an intermediate evaluation.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto schurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( serial( rhs.sm_ ) );
      schurAssign( *lhs, map( tmp, rhs.op_ ) );
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
   /*!\brief SMP assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a sparse matrix map
   // expression to a dense matrix. Due to the explicit application of the SFINAE principle,
   // this function can only be selected by the compiler in case the the expression specific
   // parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   // No special implementation for the SMP assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a sparse
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpAddAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a sparse
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline auto smpSubAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpSubAssign( *lhs, map( tmp, rhs.op_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix map expression to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side map expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a sparse
   // matrix map expression to a dense matrix. Due to the explicit application of the SFINAE
   // principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline auto smpSchurAssign( DenseMatrix<MT2,SO2>& lhs, const SMatMapExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RT );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( RT, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( RT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const RT tmp( rhs.sm_ );
      smpSchurAssign( *lhs, map( tmp, rhs.op_ ) );
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
/*!\brief Evaluates the given custom operation on each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a sm.
//
// The \a map() function evaluates the given custom operation on each non-zero element of the
// input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a map() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = map( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto) map( const SparseMatrix<MT,SO>& sm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const SMatMapExpr<MT,OP,SO>;
   return ReturnType( *sm, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluates the given custom operation on each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \param op The custom operation.
// \return The custom operation applied to each single element of \a sm.
//
// The \a forEach() function evaluates the given custom operation on each non-zero element of the
// input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a forEach() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = forEach( A, []( double a ){ return std::sqrt( a ); } );
   \endcode
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , typename OP >  // Type of the custom operation
inline decltype(auto) forEach( const SparseMatrix<MT,SO>& sm, OP op )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, std::move(op) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a abs() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the abs() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a abs() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = abs( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) abs( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Abs() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a sign() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the sign() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sign() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = sign( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sign( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Sign() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a floor() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the floor() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a floor() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = floor( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) floor( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Floor() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a ceil() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the ceil() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a ceil() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = ceil( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) ceil( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Ceil() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a trunc() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the trunc() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a trunc() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = trunc( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) trunc( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Trunc() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Applies the \a round() function to each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// This function applies the round() function to each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a round() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = round( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) round( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Round() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the complex conjugate of each single element of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The complex conjugate of each single element of \a sm.
//
// The \a conj() function calculates the complex conjugate of each element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a conj() function:

   \code
   blaze::CompressedMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) conj( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Conj() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the conjugate transpose matrix of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The conjugate transpose of \a sm.
//
// The \a ctrans() function returns an expression representing the conjugate transpose (also
// called adjoint matrix, Hermitian conjugate matrix or transjugate matrix) of the given input
// matrix \a sm.\n
// The following example demonstrates the use of the \a ctrans() function:

   \code
   blaze::CompressedMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = ctrans( A );
   \endcode

// Note that the \a ctrans() function has the same effect as manually applying the \a conj() and
// \a trans function in any order:

   \code
   B = trans( conj( A ) );  // Computing the conjugate transpose matrix
   B = conj( trans( A ) );  // Computing the conjugate transpose matrix
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) ctrans( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return trans( conj( *sm ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the real parts of each single element of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The real part of each single element of \a sm.
//
// The \a real() function calculates the real part of each element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a real() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = real( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) real( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Real() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the imaginary parts of each single element of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The imaginary part of each single element of \a sm.
//
// The \a imag() function calculates the imaginary part of each element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a imag() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = imag( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) imag( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Imag() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns a matrix containing the phase angle of each single element of \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The phase angle of each single element of \a sm.
//
// The \a arg() function calculates the phase angle of each element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a arg() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = arg( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) arg( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Arg() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the square root of each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[0..\infty)\f$.
// \return The square root of each single element of \a sm.
//
// The \a sqrt() function computes the square root of each non-zero element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sqrt() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = sqrt( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sqrt( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Sqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse square root of each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$(0..\infty)\f$.
// \return The inverse square root of each single element of \a sm.
//
// The \a invsqrt() function computes the inverse square root of each non-zero element of the
// input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invsqrt() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = invsqrt( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$(0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) invsqrt( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, InvSqrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cubic root of each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[0..\infty)\f$.
// \return The cubic root of each single element of \a sm.
//
// The \a cbrt() function computes the cubic root of each non-zero element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cbrt() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = cbrt( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) cbrt( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Cbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cubic root of each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$(0..\infty)\f$.
// \return The inverse cubic root of each single element of \a sm.
//
// The \a invcbrt() function computes the inverse cubic root of each non-zero element of the
// input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a invcbrt() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = invcbrt( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$(0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) invcbrt( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, InvCbrt() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Restricts each single element of the sparse matrix \a sm to the range \f$[min..max]\f$.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \param min The lower delimiter.
// \param max The upper delimiter.
// \return The matrix with restricted elements.
//
// The \a clamp() function restricts each element of the input matrix \a sm to the range
// \f$[min..max]\f$. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a clamp() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = clamp( A, -1.0, 1.0 );
   \endcode
*/
template< typename MT    // Type of the sparse matrix
        , bool SO        // Storage order
        , typename DT >  // Type of the delimiters
inline decltype(auto) clamp( const SparseMatrix<MT,SO>& sm, const DT& min, const DT& max )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, bind2nd( bind3rd( Clamp(), max ), min ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the exponential value for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \param exp The scalar exponent.
// \return The exponential value of each non-zero element of \a sm.
//
// The \a pow() function computes the exponential value for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a pow() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = pow( A, 4.2 );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO      // Storage order
        , typename ST  // Type of the scalar exponent
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) pow( const SparseMatrix<MT,SO>& sm, ST exp )
{
   BLAZE_FUNCTION_TRACE;

   using ScalarType = MultTrait_t< UnderlyingBuiltin_t<MT>, ST >;
   return map( *sm, blaze::bind2nd( Pow(), ScalarType( exp ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ e^x \f$ for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// The \a exp() function computes \f$ e^x \f$ for each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) exp( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Exp() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 2^x \f$ for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// The \a exp2() function computes \f$ 2^x \f$ for each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp2() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp2( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) exp2( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Exp2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes \f$ 10^x \f$ for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The resulting sparse matrix.
//
// The \a exp10() function computes \f$ 10^x \f$ for each non-zero element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a exp10() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = exp10( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) exp10( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Exp10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of each non-zero element of \a sm.
//
// The \a log() function computes the natural logarithm for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = log( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) log( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Log() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the binary logarithm for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[0..\infty)\f$.
// \return The binary logarithm of each non-zero element of \a sm.
//
// The \a log2() function computes the binary logarithm for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log2() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = log2( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) log2( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Log2() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the common logarithm for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[0..\infty)\f$.
// \return The common logarithm of each non-zero element of \a sm.
//
// The \a log10() function computes the common logarithm for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a log10() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = log10( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[0..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) log10( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Log10() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of x+1 for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all elements must be in the range \f$[-1..\infty)\f$.
// \return The natural logarithm of x+1 for each non-zero element of \a sm.
//
// The \a log1p() function computes the natural logarithm of x+1 for each non-zero element of
// the input matrix \a sm. This may be preferred over the natural logarithm for higher precision
// computing the natural logarithm of a quantity very close to 1. The function returns an
// expression representing this operation.\n
// The following example demonstrates the use of the \a log1p() function:

   \code
   blaze::DynamicMatrix<double> A, B;
   // ... Resizing and initialization
   B = log1p( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) log1p( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Log1p() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the natural logarithm of the absolute value of the gamma function for each
//        non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all elements must be in the range \f$[0..\infty)\f$.
// \return The natural logarithm of the absolute value of the gamma function of each non-zero element of \a sm.
//
// The \a lgamma() function computes the natural logarithm of the absolute value of the gamma
// function for each non-zero element of the input matrix \a sm. The function returns an
// expression representing this operation.\n
// The following example demonstrates the use of the \a lgamma() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = lgamma( A );
   \endcode

// \note All elements are expected to be in the range \f$[0..\infty)\f$. No runtime checks are
// performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) lgamma( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, LGamma() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the sine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The sine of each non-zero element of \a sm.
//
// The \a sin() function computes the sine for each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sin() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = sin( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sin( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Sin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse sine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[-1..1]\f$.
// \return The inverse sine of each non-zero element of \a sm.
//
// The \a asin() function computes the inverse sine for each non-zero element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asin() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = asin( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[-1..1]\f$. No runtime checks
// are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) asin( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Asin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic sine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The hyperbolic sine of each non-zero element of \a sm.
//
// The \a sinh() function computes the hyperbolic sine for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a sinh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = sinh( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sinh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Sinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic sine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The inverse hyperbolic sine of each non-zero element of \a sm.
//
// The \a asinh() function computes the inverse hyperbolic sine for each non-zero element of
// the input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a asinh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = asinh( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) asinh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Asinh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the cosine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The cosine of each non-zero element of \a sm.
//
// The \a cos() function computes the cosine for each non-zero element of the input matrix
// \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cos() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = cos( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) cos( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Cos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse cosine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[-1..1]\f$.
// \return The inverse cosine of each non-zero element of \a sm.
//
// The \a acos() function computes the inverse cosine for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acos() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = acos( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[-1..1]\f$. No runtime checks
// are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) acos( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Acos() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic cosine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The hyperbolic cosine of each non-zero element of \a sm.
//
// The \a cosh() function computes the hyperbolic cosine for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a cosh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = cosh( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) cosh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Cosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic cosine for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[1..\infty)\f$.
// \return The inverse hyperbolic cosine of each non-zero element of \a sm.
//
// The \a acosh() function computes the inverse hyperbolic cosine for each non-zero element of
// the input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a acosh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = acosh( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[1..\infty)\f$. No runtime
// checks are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) acosh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Acosh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the tangent for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The tangent of each non-zero element of \a sm.
//
// The \a tan() function computes the tangent for each non-zero element of the input matrix \a sm.
// The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tan() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = tan( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) tan( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Tan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse tangent for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The inverse tangent of each non-zero element of \a sm.
//
// The \a atan() function computes the inverse tangent for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atan() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = atan( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) atan( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Atan() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hyperbolic tangent for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[-1..1]\f$.
// \return The hyperbolic tangent of each non-zero element of \a sm.
//
// The \a tanh() function computes the hyperbolic tangent for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a tanh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = tanh( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[-1..1]\f$. No runtime checks
// are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) tanh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Tanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the inverse hyperbolic tangent for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix; all non-zero elements must be in the range \f$[-1..1]\f$.
// \return The inverse hyperbolic tangent of each non-zero element of \a sm.
//
// The \a atanh() function computes the inverse hyperbolic tangent for each non-zero element of
// the input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a atanh() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = atanh( A );
   \endcode

// \note All non-zero elements are expected to be in the range \f$[-1..1]\f$. No runtime checks
// are performed to assert this precondition!
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) atanh( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Atanh() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the error function for each non-zero element of the sparse matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The error function of each non-zero element of \a sm.
//
// The \a erf() function computes the error function for each non-zero element of the input
// matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erf() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = erf( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) erf( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Erf() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the complementary error function for each non-zero element of the sparse
//        matrix \a sm.
// \ingroup sparse_matrix
//
// \param sm The input matrix.
// \return The complementary error function of each non-zero element of \a sm.
//
// The \a erfc() function computes the complementary error function for each non-zero element of
// the input matrix \a sm. The function returns an expression representing this operation.\n
// The following example demonstrates the use of the \a erfc() function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = erfc( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) erfc( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return map( *sm, Erfc() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Absolute value function for absolute value sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param sm The absolute value sparse matrix expression.
// \return The absolute value of each single element of \a sm.
//
// This function implements a performance optimized treatment of the absolute value operation
// on a sparse matrix absolute value expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) abs( const SMatMapExpr<MT,Abs,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a sign() function to a sparse matrix \a sign() expressions.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix \a sign() expression.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the \a sign() operation on
// a sparse matrix \a sign() expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) sign( const SMatMapExpr<MT,Sign,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a floor() function to a sparse matrix \a floor() expressions.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix \a floor() expression.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the \a floor() operation on
// a sparse matrix \a floor() expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) floor( const SMatMapExpr<MT,Floor,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a ceil() function to a sparse matrix \a ceil() expressions.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix \a ceil() expression.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the \a ceil() operation on
// a sparse matrix \a ceil() expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) ceil( const SMatMapExpr<MT,Ceil,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a trunc() function to a sparse matrix \a trunc() expressions.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix \a trunc() expression.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the \a trunc() operation on
// a sparse matrix \a trunc() expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) trunc( const SMatMapExpr<MT,Trunc,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Applies the \a round() function to a sparse matrix \a round() expressions.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix \a round() expression.
// \return The resulting sparse matrix.
//
// This function implements a performance optimized treatment of the \a round() operation on
// a sparse matrix \a round() expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) round( const SMatMapExpr<MT,Round,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for complex conjugate sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param sm The complex conjugate sparse matrix expression.
// \return The original sparse matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a sparse matrix complex conjugate expression. It returns an expression representing the
// original sparse matrix:

   \code
   blaze::CompressedMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( conj( A ) );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool TF >    // Transpose flag
inline decltype(auto) conj( const SMatMapExpr<MT,Conj,TF>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Complex conjugate function for conjugate transpose sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param dm The conjugate transpose sparse matrix expression.
// \return The transpose sparse matrix.
//
// This function implements a performance optimized treatment of the complex conjugate operation
// on a sparse matrix conjugate transpose expression. It returns an expression representing the
// transpose of the sparse matrix:

   \code
   blaze::CompressedMatrix< complex<double> > A, B;
   // ... Resizing and initialization
   B = conj( ctrans( A ) );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) conj( const SMatTransExpr<SMatMapExpr<MT,Conj,SO>,!SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return trans( sm.operand().operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Real part function for real part sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param sm The real part sparse matrix expression.
// \return The real part of each single element of \a sm.
//
// This function implements a performance optimized treatment of the real part operation on
// a sparse matrix real part expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) real( const SMatMapExpr<MT,Real,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Imaginary part function for imaginary part sparse matrix expressions.
// \ingroup sparse_matrix
//
// \param sm The imaginary part sparse matrix expression.
// \return The imaginary part of each single element of \a sm.
//
// This function implements a performance optimized treatment of the imaginary part operation
// on a sparse matrix imaginary part expression.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) imag( const SMatMapExpr<MT,Imag,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   return sm;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
