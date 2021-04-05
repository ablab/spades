//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatDMatSchurExpr.h
//  \brief Header file for the sparse matrix/dense matrix Schur product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATDMATSCHUREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATDMATSCHUREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SchurExpr.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Max.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS SMATDMATSCHUREXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for sparse matrix-dense matrix Schur product.
// \ingroup sparse_matrix_expression
//
// The SMatDMatSchurExpr class represents the compile time expression for Schur products between
// a row-major sparse matrix and a dense matrix.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
class SMatDMatSchurExpr
   : public SchurExpr< SparseMatrix< SMatDMatSchurExpr<MT1,MT2>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side sparse matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side dense matrix expression.
   using RN1 = ReturnType_t<MT1>;     //!< Return type of the left-hand side sparse matrix expression.
   using RN2 = ReturnType_t<MT2>;     //!< Return type of the right-hand side dense matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side sparse matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either matrix operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = ( !IsTemporary_v<RN1> && !IsTemporary_v<RN2> );

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN1>() * std::declval<RN2>() );
   //**********************************************************************************************

   //**Evaluation strategy*************************************************************************
   //! Compilation switch for the evaluation strategy of the Schur product expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the Schur product expression. In case either the sparse or
       the dense matrix operand requires an intermediate evaluation, \a useAssign will be set
       to \a true and the Schur product expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to \a false and the expression will be
       evaluated via the function call operator. */
   static constexpr bool useAssign = ( RequiresEvaluation_v<MT1> || RequiresEvaluation_v<MT2> );

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatDMatSchurExpr instance.
   using This = SMatDMatSchurExpr<MT1,MT2>;

   //! Base type of this SMatDMatSchurExpr instance.
   using BaseType = SchurExpr< SparseMatrix<This,false> >;

   using ResultType    = SchurTrait_t<RT1,RT2>;        //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SMatDMatSchurExpr& >;

   //! Composite type of the left-hand side sparse matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side dense matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the sparse matrix/dense matrix Schur product expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse matrix expression.
      using Element = ValueIndexPair<ElementType>;

      //! Iterator type of the sparse matrix expression.
      using LeftIterator = ConstIterator_t< RemoveReference_t<LeftOperand> >;

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
      // \param left Iterator to the current position in the left-hand side sparse matrix expression.
      // \param right Handle to the right-hand side dense matrix expression.
      // \param row The row index of the given iterator.
      */
      inline ConstIterator( LeftIterator left, RightOperand right, size_t row )
         : left_ ( left  )  // Iterator over the elements of the left-hand side sparse matrix expression
         , right_( right )  // Right-hand side dense matrix expression
         , row_  ( row   )  // The row index of the iterator
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented expression iterator.
      */
      inline ConstIterator& operator++() {
         ++left_;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The element at the current iterator position.
      */
      inline const Element operator*() const {
         return Element( left_->value() * right_(row_,left_->index()), left_->index() );
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
         return left_->value() * right_(row_,left_->index());
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return left_->index();
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return left_ == rhs.left_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side expression iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return left_ != rhs.left_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two expression iterators.
      //
      // \param rhs The right-hand side expression iterator.
      // \return The number of elements between the two expression iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return left_ - rhs.left_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      LeftIterator left_;   //!< Iterator over the elements of the left-hand side sparse matrix expression.
      RightOperand right_;  //!< Right-hand side dense matrix expression.
      size_t       row_;    //!< The row index of the iterator.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatDMatSchurExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the Schur product expression.
   // \param rhs The right-hand side dense matrix operand of the Schur product expression.
   */
   inline SMatDMatSchurExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse matrix of the Schur product expression
      , rhs_( rhs )  // Right-hand side dense matrix of the Schur product expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.columns(), "Invalid number of columns" );
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
      BLAZE_INTERNAL_ASSERT( i < lhs_.rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < lhs_.columns(), "Invalid column access index" );
      return lhs_(i,j) * rhs_(i,j);
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
      if( i >= lhs_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= lhs_.columns() ) {
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
      return ConstIterator( lhs_.begin(i), rhs_, i );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row \a i.
   //
   // \param i The row index.
   // \return Iterator just past the last non-zero element of row \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( lhs_.end(i), rhs_, i );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return lhs_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return lhs_.columns();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   inline size_t nonZeros() const {
      return lhs_.nonZeros();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const {
      return lhs_.nonZeros(i);
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
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT1 );
      return ConstIterator( lhs_.find( i, j ), rhs_, i );
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
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT1 );
      return ConstIterator( lhs_.lowerBound( i, j ), rhs_, i );
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
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT1 );
      return ConstIterator( lhs_.upperBound( i, j ), rhs_, i );
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side sparse matrix operand.
   //
   // \return The left-hand side sparse matrix operand.
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
   // \return \a true in case the expression can alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( lhs_.canAlias( alias ) || rhs_.isAliased( alias ) );
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
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the Schur product expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the Schur product expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-dense matrix Schur product to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-dense
   // matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO2>& lhs, const SMatDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      for( size_t i=0UL; i<(*lhs).rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            (*lhs)(i,element->index()) = element->value() * B(i,element->index());
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-dense matrix Schur product to a row-major sparse matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-dense
   // matrix Schur product expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,false>& lhs, const SMatDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      // Final memory allocation (based on the evaluated operands)
      (*lhs).reserve( A.nonZeros() );

      // Performing the Schur product
      for( size_t i=0UL; i<(*lhs).rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            (*lhs).append( i, element->index(), element->value() * B(i,element->index()) );
         (*lhs).finalize( i );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a sparse matrix-dense matrix Schur product to a column-major sparse
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side Schur product expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a sparse matrix-dense
   // matrix Schur product expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,true>& lhs, const SMatDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      const size_t m( rhs.rows()    );
      const size_t n( rhs.columns() );

      // Counting the number of elements per column
      std::vector<size_t> nonzeros( n, 0UL );
      for( size_t i=0UL; i<m; ++i ) {
         const auto lend( A.end(i) );
         for( auto l=A.begin(i); l!=lend; ++l ) {
            ++nonzeros[l->index()];
         }
      }

      // Resizing the left-hand side sparse matrix
      for( size_t j=0UL; j<n; ++j ) {
         (*lhs).reserve( j, nonzeros[j] );
      }

      // Performing the Schur product
      for( size_t i=0UL; i<m; ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            (*lhs).append( i, element->index(), element->value() * B(i,element->index()) );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a sparse matrix-dense matrix Schur product to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a sparse matrix-
   // dense matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO2>& lhs, const SMatDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      for( size_t i=0UL; i<(*lhs).rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            (*lhs)(i,element->index()) += element->value() * B(i,element->index());
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a sparse matrix-dense matrix Schur product to a dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a sparse matrix-
   // dense matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO2>& lhs, const SMatDMatSchurExpr& rhs )
      -> EnableIf_t< UseAssign_v<MT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      for( size_t i=0UL; i<(*lhs).rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            (*lhs)(i,element->index()) -= element->value() * B(i,element->index());
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a sparse matrix-dense matrix Schur product to a dense
   //        matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a sparse
   // matrix-dense matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO2 >   // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO2>& lhs, const SMatDMatSchurExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).columns()  , "Invalid number of columns" );

      for( size_t i=0UL; i<(*lhs).rows(); ++i )
      {
         const auto end( A.end(i) );
         auto begin( A.begin(i) );
         size_t j( 0UL );

         for( ; begin!=end; ++begin ) {
            const size_t index( begin->index() );
            for( ; j<index; ++j )
               reset( (*lhs)(i,j) );
            (*lhs)(i,index) *= begin->value() * B(i,index);
            ++j;
         }

         for( ; j<(*lhs).columns(); ++j )
            reset( (*lhs)(i,j) );
      }
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
   // No special implementation for the SMP addition assignment to dense matrices.
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   // No special implementation for the SMP subtraction assignment to dense matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a sparse matrix-dense matrix Schur product to a
   //        dense matrix.
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side Schur product expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // sparse matrix-dense matrix Schur product expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const SMatDMatSchurExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSchurAssign( *lhs, rhs.lhs_ );
      smpSchurAssign( *lhs, rhs.rhs_ );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_SCHUREXPR( MT1, MT2 );
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
/*!\brief Backend implementation of the Schur product between a row-major sparse matrix and
//        a row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product between a
// row-major sparse matrix and a row-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , DisableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                       ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) ||
                       ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                       ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                       ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                       ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                       IsZero_v<MT1> >* = nullptr >
inline const SMatDMatSchurExpr<MT1,MT2>
   smatdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return SMatDMatSchurExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a unitriangular row-major sparse
//        matrix and a unitriangular row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the Schur product between
// a unitriangular row-major sparse matrix and a unitriangular row-major dense matrix. It
// returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                      ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   smatdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE( ReturnType );

   return ReturnType( (*lhs).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a (strictly) triangular row-major
//        sparse matrix and a (strictly) triangular row-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// (strictly) triangular row-major sparse matrix and a (strictly) triangular row-major dense
// matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                      ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                      ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                      ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                      IsZero_v<MT1> >* = nullptr >
inline decltype(auto)
   smatdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).rows(), (*lhs).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Operator for the Schur product of a row-major sparse matrix and a row-major dense
//        matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the Schur product of a row-major sparse matrix and a row-major dense
// matrix:

   \code
   using blaze::rowMajor;

   blaze::CompressedMatrix<double,rowMajor> A, C;
   blaze::DynamicMatrix<double,rowMajor> B;
   // ... Resizing and initialization
   C = A % B;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator%( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return smatdmatschur( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a row-major sparse matrix and
//        a column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
//
// This function implements a performance optimized treatment of the Schur product between a
// row-major sparse matrix and a column-major dense matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , DisableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                       ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) ||
                       ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                       ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                       ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                       ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                       IsZero_v<MT1> >* = nullptr >
inline const SMatDMatSchurExpr<MT1,MT2>
   smattdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   return SMatDMatSchurExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a unitriangular row-major sparse
//        matrix and a unitriangular column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// unitriangular row-major sparse matrix and a unitriangular column-major dense matrix. It
// returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsUniLower_v<MT1> && IsUniUpper_v<MT2> ) ||
                      ( IsUniUpper_v<MT1> && IsUniLower_v<MT2> ) >* = nullptr >
inline decltype(auto)
   smattdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE( ReturnType );

   return ReturnType( (*lhs).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the Schur product between a (strictly) triangular row-major
//        sparse matrix and a (strictly) triangular column-major dense matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the Schur product between a
// (strictly) triangular row-major sparse matrix and a (strictly) triangular column-major dense
// matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , EnableIf_t< ( IsStrictlyLower_v<MT1> && IsUpper_v<MT2> ) ||
                      ( IsStrictlyUpper_v<MT1> && IsLower_v<MT2> ) ||
                      ( IsLower_v<MT1> && IsStrictlyUpper_v<MT2> ) ||
                      ( IsUpper_v<MT1> && IsStrictlyLower_v<MT2> ) ||
                      IsZero_v<MT1> >* = nullptr >
inline decltype(auto)
   smattdmatschur( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   using ReturnType = const SchurTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).rows(), (*lhs).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Operator for the Schur product of a row-major sparse matrix and a column-major dense
//        matrix (\f$ A=B \circ C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the Schur product.
// \param rhs The right-hand side dense matrix for the Schur product.
// \return The Schur product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the Schur product of a row-major sparse matrix and a column-major
// dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,rowMajor> A, C;
   blaze::DynamicMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = A % B;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator%( const SparseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return smattdmatschur( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2 >
struct Size< SMatDMatSchurExpr<MT1,MT2>, 0UL >
   : public Max_t< Size<MT1,0UL>, Size<MT2,0UL> >
{};

template< typename MT1, typename MT2 >
struct Size< SMatDMatSchurExpr<MT1,MT2>, 1UL >
   : public Max_t< Size<MT1,1UL>, Size<MT2,1UL> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
