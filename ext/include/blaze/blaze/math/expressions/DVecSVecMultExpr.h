//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecSVecMultExpr.h
//  \brief Header file for the dense vector/sparse vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/constraints/VecVecMultExpr.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/VecVecMultExpr.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsTemporary.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Max.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DVECSVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense vector-sparse vector multiplications.
// \ingroup sparse_vector_expression
//
// The DVecSVecMultExpr class represents the compile time expression for componentwise
// multiplications between a dense vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
class DVecSVecMultExpr
   : public VecVecMultExpr< SparseVector< DVecSVecMultExpr<VT1,VT2,TF>, TF > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<VT1>;     //!< Result type of the left-hand side dense vector expression.
   using RT2 = ResultType_t<VT2>;     //!< Result type of the right-hand side sparse vector expression.
   using RN1 = ReturnType_t<VT1>;     //!< Return type of the left-hand side dense vector expression.
   using RN2 = ReturnType_t<VT2>;     //!< Return type of the right-hand side sparse vector expression.
   using CT1 = CompositeType_t<VT1>;  //!< Composite type of the left-hand side dense vector expression.
   using CT2 = CompositeType_t<VT2>;  //!< Composite type of the right-hand side sparse vector expression.
   using TT1 = TransposeType_t<VT1>;  //!< Transpose type of the left-hand side dense vector expression.
   using TT2 = TransposeType_t<VT2>;  //!< Transpose type of the right-hand side sparse vector expression.
   //**********************************************************************************************

   //**Return type evaluation**********************************************************************
   //! Compilation switch for the selection of the subscript operator return type.
   /*! The \a returnExpr compile time constant expression is a compilation switch for the
       selection of the \a ReturnType. If either vector operand returns a temporary vector
       or matrix, \a returnExpr will be set to \a false and the subscript operator will
       return it's result by value. Otherwise \a returnExpr will be set to \a true and
       the subscript operator may return it's result as an expression. */
   static constexpr bool returnExpr = !IsTemporary_v<RN1> && !IsTemporary_v<RN2>;

   //! Expression return type for the subscript operator.
   using ExprReturnType = decltype( std::declval<RN1>() * std::declval<RN2>() );
   //**********************************************************************************************

   //**Evaluation strategy*************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case either the dense or
       the sparse vector operand requires an intermediate evaluation, \a useAssign will be set
       to \a true and the multiplication expression will be evaluated via the \a assign function
       family. Otherwise \a useAssign will be set to \a false and the expression will be
       evaluated via the subscript operator. */
   static constexpr bool useAssign = ( RequiresEvaluation_v<VT1> || RequiresEvaluation_v<VT2> );

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecSVecMultExpr instance.
   using This = DVecSVecMultExpr<VT1,VT2,TF>;

   //! Base type of this DVecSVecMultExpr instance.
   using BaseType = VecVecMultExpr< SparseVector<This,TF> >;

   using ResultType    = MultTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = const If_t< returnExpr, ExprReturnType, ElementType >;

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DVecSVecMultExpr& >;

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = If_t< IsExpression_v<VT1>, const VT1, const VT1& >;

   //! Composite type of the right-hand side sparse vector expression.
   using RightOperand = If_t< IsExpression_v<VT2>, const VT2, const VT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector-sparse vector multiplication expression.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the sparse vector expression.
      using Element = ValueIndexPair<ElementType>;

      //! Iterator type of the sparse vector expression.
      using IteratorType = ConstIterator_t< RemoveReference_t<RightOperand> >;

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
      // \param vec Handle to the left-hand side dense vector expression.
      // \param it Iterator to the current position in the right-hand side sparse vector expression.
      */
      inline ConstIterator( LeftOperand vec, IteratorType it )
         : vec_( vec )  // Left-hand side dense vector expression
         , it_ ( it  )  // Iterator over the elements of the right-hand side sparse vector expression
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
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline const Element operator*() const {
         return Element( vec_[it_->index()] * it_->value(), it_->index() );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
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
         return vec_[it_->index()] * it_->value();
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
      LeftOperand vec_;  //!< Left-hand side dense vector expression.
      IteratorType it_;  //!< Iterator over the elements of the right-hand side sparse vector expression.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecSVecMultExpr class.
   //
   // \param lhs The left-hand side dense vector operand of the multiplication expression.
   // \param rhs The right-hand side sparse vector operand of the multiplication expression.
   */
   inline DVecSVecMultExpr( const VT1& lhs, const VT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense vector of the multiplication expression
      , rhs_( rhs )  // Right-hand side sparse vector of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.size() == rhs.size(), "Invalid vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < lhs_.size(), "Invalid vector access index" );
      return lhs_[index] * rhs_[index];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= lhs_.size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the sparse vector.
   //
   // \return Iterator to the first non-zero element of the sparse vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( lhs_, rhs_.begin() );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the sparse vector.
   //
   // \return Iterator just past the last non-zero element of the sparse vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( lhs_, rhs_.end() );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return lhs_.size();
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse vector.
   //
   // \return The number of non-zero elements in the sparse vector.
   */
   inline size_t nonZeros() const {
      return rhs_.nonZeros();
   }
   //**********************************************************************************************

   //**Find function*******************************************************************************
   /*!\brief Searches for a specific vector element.
   //
   // \param index The index of the search element.
   // \return Iterator to the element in case the index is found, end() iterator otherwise.
   */
   inline ConstIterator find( size_t index ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT2 );
      return ConstIterator( lhs_, rhs_.find( index ) );
   }
   //**********************************************************************************************

   //**LowerBound function*************************************************************************
   /*!\brief Returns an iterator to the first index not less then the given index.
   //
   // \param index The index of the search element.
   // \return Iterator to the first index not less then the given index, end() iterator otherwise.
   */
   inline ConstIterator lowerBound( size_t index ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT2 );
      return ConstIterator( lhs_, rhs_.lowerBound( index ) );
   }
   //**********************************************************************************************

   //**UpperBound function*************************************************************************
   /*!\brief Returns an iterator to the first index greater then the given index.
   //
   // \param index The index of the search element.
   // \return Iterator to the first index greater then the given index, end() iterator otherwise.
   */
   inline ConstIterator upperBound( size_t index ) const {
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT2 );
      return ConstIterator( lhs_, rhs_.upperBound( index ) );
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side sparse vector operand.
   //
   // \return The right-hand side sparse vector operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return rhs_;
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
      return ( lhs_.canAlias( alias ) || rhs_.canAlias( alias ) );
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
   LeftOperand  lhs_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector multiplication expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT,TF>& lhs, const DVecSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      for( auto element=y.begin(); element!=y.end(); ++element )
         (*lhs)[element->index()] = x[element->index()] * element->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense vector-sparse vector multiplication to a sparse vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense vector-sparse
   // vector multiplication expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT,TF>& lhs, const DVecSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      // Final memory allocation (based on the evaluated operands)
      (*lhs).reserve( y.nonZeros() );

      // Performing the vector multiplication
      for( auto element=y.begin(); element!=y.end(); ++element )
         (*lhs).append( element->index(), x[element->index()] * element->value() );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense vector-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense vector-
   // sparse vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT,TF>& lhs, const DVecSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      for( auto element=y.begin(); element!=y.end(); ++element )
         (*lhs)[element->index()] += x[element->index()] * element->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense vector-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense vector-
   // sparse vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT,TF>& lhs, const DVecSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      for( auto element=y.begin(); element!=y.end(); ++element )
         (*lhs)[element->index()] -= x[element->index()] * element->value();
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense vector-sparse vector multiplication to a dense vector.
   // \ingroup sparse_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // vector-sparse vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either of the two operands requires an intermediate evaluation.
   */
   template< typename VT >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT,TF>& lhs, const DVecSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      CT1 x( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense vector operand
      CT2 y( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse vector operand

      BLAZE_INTERNAL_ASSERT( x.size() == rhs.lhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( y.size() == rhs.rhs_.size(), "Invalid vector size" );
      BLAZE_INTERNAL_ASSERT( x.size() == (*lhs).size()  , "Invalid vector size" );

      const auto end( y.end() );
      auto begin( y.begin() );
      size_t i( 0UL );

      for( ; begin!=end; ++begin ) {
         const size_t index( begin->index() );
         for( ; i<index; ++i )
            reset( (*lhs)[i] );
         (*lhs)[index] *= x[index] * begin->value();
         ++i;
      }

      for( ; i<(*lhs).size(); ++i )
         reset( (*lhs)[i] );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT1, TF );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT2, TF );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_VECVECMULTEXPR( VT1, VT2 );
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
/*!\brief Backend implementation of the componentwise product of a dense vector and a sparse vector
//        (\f$ \vec{a}=\vec{b}*\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side dense vector for the component product.
// \param rhs The right-hand side sparse vector for the component product.
// \return The product of the two vectors.
//
// This function implements a performance optimized treatment of the componentwise multiplication
// of a dense vector and a sparse vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , DisableIf_t< IsZero_v<VT2> >* = nullptr >
inline const DVecSVecMultExpr<VT1,VT2,TF>
   dvecsvecmult( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   return DVecSVecMultExpr<VT1,VT2,TF>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the componentwise product of a dense vector and a zero vector
//        (\f$ \vec{a}=\vec{b}*\vec{c} \f$).
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the component product.
// \param rhs The right-hand side sparse vector for the component product.
// \return The resulting zero vector.
//
// This function implements a performance optimized treatment of the componentwise multiplication
// of a dense vector and a zero vector. It returns a zero vector.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF       // Transpose flag
        , EnableIf_t< IsZero_v<VT2> >* = nullptr >
inline decltype(auto)
   dvecsvecmult( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   using ReturnType = const MultTrait_t< ResultType_t<VT1>, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ReturnType, TF );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the componentwise product of a dense vector and a sparse
//        vector (\f$ \vec{a}=\vec{b}*\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side dense vector for the component product.
// \param rhs The right-hand side sparse vector for the component product.
// \return The product of the two vectors.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the componentwise multiplication of a dense vector and a sparse
// vector:

   \code
   blaze::DynamicVector<double> a;
   blaze::CompressedVector<double> b, c;
   // ... Resizing and initialization
   c = a * b;
   \endcode

// The operator returns an expression representing a sparse vector of the higher-order element
// type of the two involved vector element types \a VT1::ElementType and \a VT2::ElementType.
// Both vector types \a VT1 and \a VT2 as well as the two element types \a VT1::ElementType
// and \a VT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag
inline decltype(auto)
   operator*( const DenseVector<VT1,TF>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   return dvecsvecmult( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT1, typename VT2, bool TF >
struct Size< DVecSVecMultExpr<VT1,VT2,TF>, 0UL >
   : public Max_t< Size<VT1,0UL>, Size<VT2,0UL> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
