//=================================================================================================
/*!
//  \file blaze/math/expressions/TDVecTSMatMultExpr.h
//  \brief Header file for the transpose dense vector/transpose sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDVECTSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDVECTSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/TVecMatMultExpr.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDVECSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense vector-transpose sparse matrix multiplications.
// \ingroup dense_vector_expression
//
// The TDVecTSMatMultExpr class represents the compile time expression for multiplications
// between transpose dense vectors and column-major sparse matrices.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Type of the right-hand side sparse matrix
class TDVecTSMatMultExpr
   : public TVecMatMultExpr< DenseVector< TDVecTSMatMultExpr<VT,MT>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using VRT = ResultType_t<VT>;     //!< Result type of the left-hand side dense vector expression.
   using MRT = ResultType_t<MT>;     //!< Result type of the right-hand side sparse matrix expression.
   using VCT = CompositeType_t<VT>;  //!< Composite type of the left-hand side dense vector expression.
   using MCT = CompositeType_t<MT>;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense vector expression.
   static constexpr bool evaluateVector = ( IsComputation_v<VT> || RequiresEvaluation_v<VT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
   static constexpr bool evaluateMatrix = RequiresEvaluation_v<MT>;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case either the vector or
       the matrix operand requires an intermediate evaluation or the dense vector expression
       is a compound expression, \a useAssign will be set to \a true and the multiplication
       expression will be evaluated via the \a assign function family. Otherwise \a useAssign
       will be set to \a false and the expression will be evaluated via the subscript operator. */
   static constexpr bool useAssign = ( evaluateVector || evaluateMatrix );
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool UseAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either the vector or the matrix operand requires an intermediate evaluation, the
       variable will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   static constexpr bool UseSMPAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TDVecTSMatMultExpr instance.
   using This = TDVecTSMatMultExpr<VT,MT>;

   //! Base type of this TDVecTSMatMultExpr instance.
   using BaseType = TVecMatMultExpr< DenseVector<This,true> >;

   using ResultType    = MultTrait_t<VRT,MRT>;         //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const TDVecTSMatMultExpr& >;

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = If_t< IsExpression_v<VT>, const VT, const VT& >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Composite type of the left-hand side dense vector expression.
   using LT = If_t< evaluateVector, const VRT, VCT >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RT = If_t< evaluateMatrix, const MRT, MCT >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateVector && VT::smpAssignable && !evaluateMatrix && MT::smpAssignable );
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDVecTSMatMultExpr class.
   */
   inline TDVecTSMatMultExpr( const VT& vec, const MT& mat ) noexcept
      : vec_( vec )  // Left-hand side dense vector of the multiplication expression
      , mat_( mat )  // Right-hand side sparse matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( vec_.size() == mat.rows(), "Invalid vector and matrix sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < mat_.columns(), "Invalid vector access index" );
      return vec_ * column( mat_, index, unchecked );
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
      if( index >= mat_.columns() ) {
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
      return mat_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return vec_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side transpose sparse matrix operand.
   //
   // \return The right-hand side transpose sparse matrix operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return mat_;
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
      return ( vec_.isAliased( alias ) || mat_.isAliased( alias ) );
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
      return ( vec_.isAliased( alias ) || mat_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return vec_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( size() > SMP_TDVECTSMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vec_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand mat_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a transpose dense vector-transpose sparse matrix multiplication to a
   //        dense vector (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense vector-
   // transpose sparse matrix multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case either the left-hand side vector operand is a compound expression or the right-hand
   // side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT x( serial( rhs.vec_ ) );  // Evaluation of the left-hand side dense vector operator
      RT A( serial( rhs.mat_ ) );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      assign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a transpose dense vector-transpose sparse matrix multiplication to a
   //        sparse vector (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense vector-
   // transpose sparse matrix multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case either the left-hand side vector operand is a compound expression or the right-hand
   // side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a transpose dense vector-transpose sparse matrix multiplication
   //        to a dense vector (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due
   // to the explicit application of the SFINAE principle, this function can only be selected
   // by the compiler in case either the left-hand side vector operand is a compound expression
   // or the right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( serial( rhs.vec_ ) );  // Evaluation of the left-hand side dense vector operator
      RT A( serial( rhs.mat_ ) );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      addAssign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either the left-hand side vector operand is a compound expression
   // or the right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( serial( rhs.vec_ ) );  // Evaluation of the left-hand side dense vector operator
      RT A( serial( rhs.mat_ ) );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      subAssign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T*=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case either the left-hand side vector operand is a compound expression or the
   // right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      multAssign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**Multiplication assignment to sparse vectors*************************************************
   // No special implementation for the multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**Division assignment to dense vectors********************************************************
   /*!\brief Division assignment of a transpose dense vector-transpose sparse matrix multiplication
   //        to a dense vector (\f$ \vec{y}^T/=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due
   // to the explicit application of the SFINAE principle, this function can only be selected
   // by the compiler in case either the left-hand side vector operand is a compound expression
   // or the right-hand side matrix operand requires an intermediate evaluation.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto divAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      divAssign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**Division assignment to sparse vectors*******************************************************
   // No special implementation for the division assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP assignment to dense vectors*************************************************************
   /*!\brief SMP assignment of a transpose dense vector-transpose sparse matrix multiplication
   //        to a dense vector (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // vector-transpose sparse matrix multiplication expression to a dense vector. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      smpAssign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*!\brief SMP assignment of a transpose dense vector-transpose sparse matrix multiplication
   //        to a sparse vector (\f$ \vec{y}^T=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // vector-transpose sparse matrix multiplication expression to a sparse vector. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*!\brief SMP addition assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T+=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      smpAddAssign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*!\brief SMP subtraction assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T-=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }

      LT x( rhs.vec_ );  // Evaluation of the left-hand side dense vector operator
      RT A( rhs.mat_ );  // Evaluation of the right-hand side sparse matrix operator

      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.columns() == (*lhs).size()     , "Invalid vector size"       );

      smpSubAssign( *lhs, x * A );
   }
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*!\brief SMP multiplication assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T*=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // transpose dense vector-transpose sparse matrix multiplication expression to a dense vector.
   // Due to the explicit application of the SFINAE principle, this function can only be selected
   // by the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*!\brief SMP division assignment of a transpose dense vector-transpose sparse matrix
   //        multiplication to a dense vector (\f$ \vec{y}^T/=\vec{x}^T*A \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a transpose
   // dense vector-transpose sparse matrix multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT2 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT2,true>& lhs, const TDVecTSMatMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpDivAssign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   // No special implementation for the SMP division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_TVECMATMULTEXPR( VT, MT );
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
/*!\brief Backend implementation of the multiplication of a transpose dense vector
//        and a column-major sparse matrix (\f$ \vec{a}=B*\vec{c} \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side column-major sparse matrix for the multiplication.
// \return The resulting transpose vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// transpose dense vector and a column-major sparse matrix.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename MT  // Type of the right-hand side sparse matrix
        , DisableIf_t< ( IsIdentity_v<MT> &&
                         IsSame_v< ElementType_t<VT>, ElementType_t<MT> > ) ||
                       IsZero_v<MT> >* = nullptr >
inline const TDVecTSMatMultExpr<VT,MT>
   tdvectsmatmult( const DenseVector<VT,true>& vec, const SparseMatrix<MT,true>& mat )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*vec).size() == (*mat).rows(), "Invalid vector and matrix sizes" );

   return TDVecTSMatMultExpr<VT,MT>( *vec, *mat );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication of a transpose dense vector
//        and a column-major identity matrix (\f$ \vec{a}=B*\vec{c} \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side column-major identity matrix for the multiplication.
// \return Reference to the given dense vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// transpose dense vector and a column-major identity matrix. It returns a reference to the
// given dense vector.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename MT  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsIdentity_v<MT> &&
                      IsSame_v< ElementType_t<VT>, ElementType_t<MT> > >* = nullptr >
inline const VT&
   tdvectsmatmult( const DenseVector<VT,true>& vec, const SparseMatrix<MT,true>& mat )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( mat );

   BLAZE_INTERNAL_ASSERT( (*vec).size() == (*mat).rows(), "Invalid vector and matrix sizes" );

   return (*vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication of a transpose dense vector
//        and a column-major zero matrix (\f$ \vec{a}=B*\vec{c} \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side column-major zero matrix for the multiplication.
// \return The resulting zero vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// transpose dense vector and a column-major zero matrix. It returns a zero vector.
*/
template< typename VT  // Type of the left-hand side dense vector
        , typename MT  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsZero_v<MT> >* = nullptr >
inline decltype(auto)
   tdvectsmatmult( const DenseVector<VT,true>& vec, const SparseMatrix<MT,true>& mat )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( vec );

   BLAZE_INTERNAL_ASSERT( (*vec).size() == (*mat).rows(), "Invalid vector and matrix sizes" );

   using ReturnType = const MultTrait_t< ResultType_t<VT>, ResultType_t<MT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*mat).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a transpose dense vector and a
//        column-major sparse matrix (\f$ \vec{y}^T=\vec{x}^T*A \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side column-major sparse matrix for the multiplication.
// \return The resulting transpose vector.
// \exception std::invalid_argument Vector and matrix sizes do not match.
//
// This operator represents the multiplication between a transpose dense vector and a column-major
// sparse matrix:

   \code
   using blaze::rowVector;
   using blaze::columnMajor;

   blaze::DynamicVector<double,rowVector> x, y;
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   y = x * A;
   \endcode

// The operator returns an expression representing a transpose dense vector of the higher-order
// element type of the two involved element types \a VT::ElementType and \a MT::ElementType.
// Both the dense matrix type \a VT and the dense vector type \a MT as well as the two element
// types \a VT::ElementType and \a MT::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of rows of
// the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const DenseVector<VT,true>& vec, const SparseMatrix<MT,true>& mat )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATMATMULTEXPR_TYPE( MT );

   if( (*vec).size() != (*mat).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector and matrix sizes do not match" );
   }

   return tdvectsmatmult( *vec, *mat );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, typename MT >
struct IsAligned< TDVecTSMatMultExpr<VT,MT> >
   : public IsAligned<VT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
