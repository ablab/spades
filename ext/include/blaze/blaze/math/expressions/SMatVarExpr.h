//=================================================================================================
/*!
//  \file blaze/math/expressions/SMatVarExpr.h
//  \brief Header file for the sparse matrix variance expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SMATVAREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SMATVAREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatReduceExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/functors/InvAdd.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base template for the row-/column-wise variance computations on row-major sparse matrices.
// \ingroup dense_vector_expression
//
// The SMatVarExpr class represents the compile time expression for the computation of the
// row-/column-wise variance function on row-major sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , ReductionFlag RF >  // Reduction flag
class SMatVarExpr
{};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR THE COLUMN-WISE VARIANCE FUNCTION ON ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the column-wise variance function on row-major sparse matrices.
// \ingroup dense_vector_expression
//
// This specialization of the SMatVarExpr class template represents the compile time expression
// for the column-wise variance function on row-major sparse matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SMatVarExpr<MT,columnwise>
   : public MatReduceExpr< DenseVector< SMatVarExpr<MT,columnwise>, true >, columnwise >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;     //!< Result type of the sparse matrix expression.
   using OT = OppositeType_t<MT>;   //!< Opposite type of the sparse matrix expression.
   using CT = CompositeType_t<MT>;  //!< Composite type of the sparse matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatVarExpr instance.
   using This = SMatVarExpr<MT,columnwise>;

   //! Base type of this SMatVarExpr instance.
   using BaseType = MatReduceExpr< DenseVector<This,true>, columnwise >;

   using ResultType    = ReduceTrait_t<RT,InvAdd,columnwise>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;          //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;            //!< Resulting element type.
   using ReturnType    = const ElementType;                    //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;                     //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatVarExpr class.
   //
   // \param sm The matrix operand of the variance expression.
   */
   explicit inline SMatVarExpr( const MT& sm ) noexcept
      : sm_( sm )  // Sparse matrix of the variance expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < sm_.columns(), "Invalid vector access index" );
      return var( column( sm_, index, unchecked ) );
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
      if( index >= sm_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return sm_.columns();
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
   // \return \a true in case an aliasing effect is possible, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
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
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sm_;  //!< Sparse matrix of the variance expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-wise row-major sparse matrix variance operation to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a column-wise row-major
   // sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      assign( *lhs, var<columnwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a column-wise row-major sparse matrix variance operation to a dense
   //        vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a column-wise row-major
   // sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( SparseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a column-wise row-major sparse matrix variance operation to a
   //        dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      addAssign( *lhs, var<columnwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( SparseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      addAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      subAssign( *lhs, var<columnwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( SparseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      subAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      multAssign( *lhs, var<columnwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( SparseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( OT );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( OT );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const OT tmp( serial( rhs.sm_ ) );
      divAssign( *lhs, var<columnwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a column-wise row-major sparse matrix variance operation
   //        to a dense vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side variance expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a column-wise
   // row-major sparse matrix variance expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( SparseVector<VT1,true>& lhs, const SMatVarExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      divAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR THE ROW-WISE VARIANCE FUNCTION ON ROW-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the row-wise variance function on row-major sparse matrices.
// \ingroup dense_vector_expression
//
// This specialization of the SMatVarExpr class template represents the compile time expression
// for the row-wise variance function on row-major sparse matrices.
*/
template< typename MT >  // Type of the sparse matrix
class SMatVarExpr<MT,rowwise>
   : public MatReduceExpr< DenseVector< SMatVarExpr<MT,rowwise>, false >, rowwise >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;  //!< Result type of the sparse matrix expression.
   //**********************************************************************************************

   //**Serial evaluation strategy******************************************************************
   //! Compilation switch for the serial evaluation strategy of the variance expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the serial evaluation strategy of the variance expression. In case the sparse matrix
       operand requires an intermediate evaluation, \a useAssign will be set to 1 and the
       variance expression will be evaluated via the \a assign function family. Otherwise
       \a useAssign will be set to 0 and the expression will be evaluated via the subscript
       operator. */
   static constexpr bool useAssign = RequiresEvaluation_v<MT>;

   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this SMatVarExpr instance.
   using This = SMatVarExpr<MT,rowwise>;

   //! Base type of this SMatVarExpr instance.
   using BaseType = MatReduceExpr< DenseVector<This,false>, rowwise >;

   using ResultType    = ReduceTrait_t<RT,InvAdd,rowwise>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;       //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;         //!< Resulting element type.
   using ReturnType    = const ElementType;                 //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const SMatVarExpr& >;

   //! Composite type of the left-hand side sparse matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense vector.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
      using ValueType        = ElementType;                      //!< Type of the underlying elements.
      using PointerType      = ElementType*;                     //!< Pointer return type.
      using ReferenceType    = ElementType&;                     //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param sm The sparse matrix operand of the variance expression.
      // \param index Index to the initial matrix row.
      */
      inline ConstIterator( Operand sm, size_t index )
         : sm_   ( sm    )  // Sparse matrix of the variance expression
         , index_( index )  // Index to the current matrix row
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ConstIterator& operator+=( size_t inc ) {
         index_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline ConstIterator& operator-=( size_t dec ) {
         index_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ConstIterator& operator++() {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator++( int ) {
         return ConstIterator( index_++ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ConstIterator& operator--() {
         --index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ConstIterator operator--( int ) {
         return ConstIterator( index_-- );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReturnType operator*() const {
         return var( row( sm_, index_, unchecked ) );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ConstIterator& rhs ) const {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ConstIterator& rhs ) const {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ConstIterator& rhs ) const {
         return index_ < rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ConstIterator& rhs ) const {
         return index_ > rhs.index_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ConstIterator& rhs ) const {
         return index_ <= rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ConstIterator& rhs ) const {
         return index_ >= rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ConstIterator& rhs ) const {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         return ConstIterator( it.index_ + inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         return ConstIterator( it.index_ + inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         return ConstIterator( it.index_ - dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      Operand sm_;     //!< Sparse matrix of the variance expression.
      size_t  index_;  //!< Index to the current matrix row.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the SMatVarExpr class.
   //
   // \param sm The matrix operand of the variance expression.
   */
   explicit inline SMatVarExpr( const MT& sm ) noexcept
      : sm_( sm )  // Sparse matrix of the variance expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < sm_.rows(), "Invalid vector access index" );
      return var( row( sm_, index, unchecked ) );
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
      if( index >= sm_.rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of the dense vector.
   //
   // \return Iterator to the first non-zero element of the dense vector.
   */
   inline ConstIterator begin() const {
      return ConstIterator( sm_, 0UL );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of the dense vector.
   //
   // \return Iterator just past the last non-zero element of the dense vector.
   */
   inline ConstIterator end() const {
      return ConstIterator( sm_, size() );
   }
   //**********************************************************************************************

   //**Size function*******************************************************************************
   /*!\brief Returns the current size/dimension of the vector.
   //
   // \return The size of the vector.
   */
   inline size_t size() const noexcept {
      return sm_.rows();
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
   // \return \a true in case an aliasing effect is possible, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return ( sm_.isAliased( alias ) );
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
      return ( sm_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand sm_;  //!< Sparse matrix of the variance expression.
   //**********************************************************************************************

   //**Assignment to vectors***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a row-wise row-major sparse matrix variance operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side variance expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a row-wise row-major
   // sparse matrix variance expression to a vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto assign( Vector<VT1,false>& lhs, const SMatVarExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      assign( *lhs, var<rowwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to vectors**************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a row-wise row-major sparse matrix variance operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side variance expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a row-wise
   // row-major sparse matrix variance expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto addAssign( Vector<VT1,false>& lhs, const SMatVarExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      addAssign( *lhs, var<rowwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to vectors***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a row-wise row-major sparse matrix variance operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side variance expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a row-wise
   // row-major sparse matrix variance expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto subAssign( Vector<VT1,false>& lhs, const SMatVarExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      subAssign( *lhs, var<rowwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a row-wise row-major sparse matrix variance operation
   //        to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side variance expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a row-wise
   // row-major sparse matrix variance expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto multAssign( Vector<VT1,false>& lhs, const SMatVarExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      multAssign( *lhs, var<rowwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Division assignment to vectors**************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Division assignment of a row-wise row-major sparse matrix variance operation to a vector.
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side vector.
   // \param rhs The right-hand side variance expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a row-wise
   // row-major sparse matrix variance expression to a vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target vector
   friend inline auto divAssign( Vector<VT1,false>& lhs, const SMatVarExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const RT tmp( serial( rhs.sm_ ) );  // Evaluation of the sparse matrix operand
      divAssign( *lhs, var<rowwise>( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief Backend implementation of the \c var() function for general sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given general sparse matrix for the variance computation.
// \return The variance of the given matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) var_backend( const SparseMatrix<MT,SO>& sm, FalseType )
{
   using BT = UnderlyingBuiltin_t<MT>;

   const size_t n ( size( *sm ) );
   const size_t nz( nonZeros( *sm ) );

   BLAZE_INTERNAL_ASSERT( n > 1UL, "Invalid matrix size detected" );
   BLAZE_INTERNAL_ASSERT( n >= nz, "Invalid number of non-zero elements detected" );

   const auto meanValue( mean( *sm ) );
   auto variance( ( n - nz ) * pow2( meanValue ) );

   const size_t iend( SO ? (*sm).columns() : (*sm).rows() );
   for( size_t i=0UL; i<iend; ++i ) {
      const auto end( (*sm).end(i) );
      for( auto element=(*sm).begin(i); element!=end; ++element ) {
         variance += pow2( element->value() - meanValue );
      }
   }

   return variance * inv( BT( n-1UL ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the \c var() function for uniform sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given uniform sparse matrix for the variance computation.
// \return The var of the given matrix.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) var_backend( const SparseMatrix<MT,SO>& sm, TrueType )
{
   MAYBE_UNUSED( sm );

   BLAZE_INTERNAL_ASSERT( size( *sm ) > 1UL, "Invalid matrix size detected" );

   return ElementType_t<MT>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the variance for the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the variance computation.
// \return The variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for
// the given sparse matrix \a sm. Both the non-zero and zero elements of the sparse matrix are
// taken into account. Example:

   \code
   using blaze::CompressedMatrix;

   CompressedMatrix<int> A{ { 1, 3, 2 }
                          , { 2, 6, 4 }
                          , { 9, 6, 3 } };

   const double v = var( A );  // Results in 6.5
   \endcode

// In case the size of the given matrix is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline decltype(auto) var( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   if( size( *sm ) < 2UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid input matrix" );
   }

   return var_backend( *sm, IsZero<MT>() );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for variance operations on row-major sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given row-major sparse matrix for the variance computation.
// \return The result of the variance operation.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT >     // Type of the sparse matrix
inline const SMatVarExpr<MT,RF> var_backend( const SparseMatrix<MT,false>& sm )
{
   using ReturnType = const SMatVarExpr<MT,RF>;
   return ReturnType( *sm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation for variance operations on column-major sparse matrices.
// \ingroup sparse_matrix
//
// \param sm The given column-major sparse matrix for the variance computation.
// \return The result of the variance operation.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT >     // Type of the sparse matrix
inline decltype(auto) var_backend( const SparseMatrix<MT,true>& sm )
{
   constexpr ReductionFlag RF2( RF == rowwise ? columnwise : rowwise );
   return trans( var<RF2>( trans( *sm ) ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the row-/column-wise variance function for the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix for the variance computation.
// \return The row-/column-wise variance of the given matrix.
// \exception std::invalid_argument Invalid input matrix.
//
// This function computes the row-/column-wise
// <a href="https://en.wikipedia.org/wiki/Variance">variance</a> for the given sparse matrix
// \a sm. In case \a RF is set to \a rowwise, the function returns a column vector containing
// the variance of each row of \a sm. In case \a RF is set to \a columnwise, the function
// returns a row vector containing the variance of each column of \a dm. Both the non-zero
// and zero elements of the sparse matrix are taken into account. Example:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicVector;
   using blaze::columnVector;
   using blaze::rowVector;

   CompressedMatrix<int> A{ { 1, 3, 2 }
                          , { 2, 6, 4 }
                          , { 9, 6, 3 } };

   DynamicVector<double,columnVector> rv;
   DynamicVector<double,rowVector> cv;

   rv = var<rowwise>( A );     // Results in ( 1  4  9 )
   cv = var<columnwise>( A );  // Results in ( 19  3  1 )
   \endcode

// In case \a RF is set to \a rowwise and the number of columns of the given matrix is smaller
// than 2 or in case \a RF is set to \a columnwise and the number of rows of the given matrix is
// smaller than 2, a \a std::invalid_argument is thrown.
*/
template< ReductionFlag RF  // Reduction flag
        , typename MT       // Type of the sparse matrix
        , bool SO >         // Storage order
inline decltype(auto) var( const SparseMatrix<MT,SO>& sm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_STATIC_ASSERT_MSG( RF < 2UL, "Invalid reduction flag" );

   return var_backend<RF>( *sm );
}
//*************************************************************************************************

} // namespace blaze

#endif
