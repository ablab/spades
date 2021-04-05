//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatSVecMultExpr.h
//  \brief Header file for the dense matrix/sparse vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATSVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATSVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/MatVecMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Vector.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUpper.h>
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
//  CLASS DMATSVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix-sparse vector multiplications.
// \ingroup dense_vector_expression
//
// The DMatSVecMultExpr class represents the compile time expression for multiplications
// between row-major dense matrices and sparse vectors.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
class DMatSVecMultExpr
   : public MatVecMultExpr< DenseVector< DMatSVecMultExpr<MT,VT>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using MRT = ResultType_t<MT>;     //!< Result type of the left-hand side dense matrix expression.
   using VRT = ResultType_t<VT>;     //!< Result type of the right-hand side sparse vector expression.
   using MCT = CompositeType_t<MT>;  //!< Composite type of the left-hand side dense matrix expression.
   using VCT = CompositeType_t<VT>;  //!< Composite type of the right-hand side sparse vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateMatrix = RequiresEvaluation_v<MT>;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse vector expression.
   static constexpr bool evaluateVector = ( IsComputation_v<VT> || RequiresEvaluation_v<VT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the evaluation strategy of the multiplication expression.
   /*! The \a useAssign compile time constant expression represents a compilation switch for
       the evaluation strategy of the multiplication expression. In case either the matrix or
       the vector operand requires an intermediate evaluation or the sparse vector expression
       is a compound expression, \a useAssign will be set to \a true and the multiplication
       expression will be evaluated via the \a assign function family. Otherwise \a useAssign
       will be set to \a false and the expression will be evaluated via the subscript operator. */
   static constexpr bool useAssign = ( evaluateMatrix || evaluateVector );
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
       In case either the matrix or the vector operand requires an intermediate evaluation, the
       variable will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   static constexpr bool UseSMPAssign_v = useAssign;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatSVecMultExpr instance.
   using This = DMatSVecMultExpr<MT,VT>;

   //! Base type of this DMatSVecMultExpr instance.
   using BaseType = MatVecMultExpr< DenseVector<This,false> >;

   using ResultType    = MultTrait_t<MRT,VRT>;         //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = If_t< useAssign, const ResultType, const DMatSVecMultExpr& >;

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Composite type of the right-hand side dense vector expression.
   using RightOperand = If_t< IsExpression_v<VT>, const VT, const VT& >;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateMatrix, const MRT, MCT >;

   //! Type for the assignment of the right-hand side dense matrix operand.
   using RT = If_t< evaluateVector, const VRT, VCT >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateMatrix && MT::smpAssignable && !evaluateVector && VT::smpAssignable );
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatSVecMultExpr class.
   //
   // \param mat The left-hand side dense matrix operand of the multiplication expression.
   // \param vec The right-hand side sparse vector operand of the multiplication expression.
   */
   inline DMatSVecMultExpr( const MT& mat, const VT& vec ) noexcept
      : mat_( mat )  // Left-hand side dense matrix of the multiplication expression
      , vec_( vec )  // Right-hand side sparse vector of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( mat_.columns() == vec_.size(), "Invalid matrix and vector sizes" );
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < mat_.rows(), "Invalid vector access index" );

      if( IsDiagonal_v<MT> )
      {
         return mat_(index,index) * vec_[index];
      }
      else if( IsLower_v<MT> )
      {
         const size_t n( IsStrictlyLower_v<MT> ? index : index+1UL );
         return subvector( row( mat_, index, unchecked ), 0UL, n, unchecked ) *
                subvector( vec_, 0UL, n, unchecked );
      }
      else if( IsUpper_v<MT> )
      {
         const size_t begin( IsStrictlyUpper_v<MT> ? index+1UL : index );
         const size_t n    ( mat_.columns() - begin );
         return subvector( row( mat_, index, unchecked ), begin, n, unchecked ) *
                subvector( vec_, begin, n, unchecked );
      }
      else
      {
         return row( mat_, index, unchecked ) * vec_;
      }
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
      if( index >= mat_.rows() ) {
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
   /*!\brief Returns the right-hand side sparse vector operand.
   //
   // \return The right-hand side sparse vector operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return vec_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the given alias is contained in this expression, \a false if not.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return mat_.isAliased( alias ) || vec_.isAliased( alias );
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
      return mat_.isAliased( alias ) || vec_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return mat_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( size() > SMP_DMATSVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side sparse vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-sparse vector multiplication to a dense vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-sparse
   // vector multiplication expression to a dense vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto assign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( serial( rhs.vec_ ) );
      if( x.nonZeros() == 0UL ) {
         reset( *lhs );
         return;
      }

      // Evaluation of the left-hand side dense matrix operand
      LT A( serial( rhs.mat_ ) );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      assign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix-sparse vector multiplication to a sparse vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix-sparse
   // vector multiplication expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline auto assign( SparseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix-sparse vector multiplication to a dense vector
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix-
   // sparse vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto addAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( serial( rhs.vec_ ) );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side dense matrix operand
      LT A( serial( rhs.mat_ ) );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      addAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix-sparse vector multiplication to a dense
   //        vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix-
   // sparse vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto subAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( serial( rhs.vec_ ) );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side dense matrix operand
      LT A( serial( rhs.mat_ ) );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      subAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a dense matrix-sparse vector multiplication to a dense
   //        vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a dense
   // matrix-sparse vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case either the left-hand side matrix operand requires an intermediate evaluation or
   // the right-hand side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto multAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   /*!\brief Division assignment of a dense matrix-sparse vector multiplication to a dense vector
   //        (\f$ \vec{y}/=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisior.
   // \return void
   //
   // This function implements the performance optimized division assignment of a dense matrix-
   // sparse vector multiplication expression to a dense vector. Due to the explicit application
   // of the SFINAE principle, this function can only be selected by the compiler in case either
   // the left-hand side matrix operand requires an intermediate evaluation or the right-hand
   // side vector operand is a compound expression.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto divAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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

   //**SMP assignment to dense vectors*************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-sparse vector multiplication to a dense vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-sparse
   // vector multiplication expression to a dense vector. Due to the explicit application of the
   // SFINAE principle, this function can only be selected by the compiler in case the expression
   // specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) {
         reset( *lhs );
         return;
      }

      // Evaluation of the left-hand side dense matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      smpAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a dense matrix-sparse vector multiplication to a sparse vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a dense matrix-sparse
   // vector multiplication expression to a sparse vector. Due to the explicit application of
   // the SFINAE principle, this function can only be selected by the compiler in case the
   // expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a dense matrix-sparse vector multiplication to a dense
   //        vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a dense
   // matrix-sparse vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side dense matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      smpAddAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a dense matrix-sparse vector multiplication to a dense
   //        vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a dense
   // matrix-sparse vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      // Evaluation of the right-hand side sparse vector operand
      RT x( rhs.vec_ );
      if( x.nonZeros() == 0UL ) return;

      // Evaluation of the left-hand side dense matrix operand
      LT A( rhs.mat_ );

      // Checking the evaluated operands
      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      // Performing the dense matrix-sparse vector multiplication
      smpSubAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a dense matrix-sparse vector multiplication to a
   //        dense vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of
   // a dense matrix-sparse vector multiplication expression to a dense vector. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpMultAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP dvision assignment of a dense matrix-sparse vector multiplication to a dense
   //        vector (\f$ \vec{y}/=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a dense
   // matrix-sparse vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT1,false>& lhs, const DMatSVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( rhs );
      smpDivAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   // No special implementation for the SMP division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATVECMULTEXPR( MT, VT );
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
/*!\brief Backend implementation of the multiplication of a row-major dense matrix and a sparse
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side row-major dense matrix for the multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// row-major dense matrix and a sparse vector.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename VT  // Type of the right-hand side sparse vector
        , DisableIf_t< IsSymmetric_v<MT> || IsZero_v<VT> >* = nullptr >
inline const DMatSVecMultExpr<MT,VT>
   dmatsvecmult( const DenseMatrix<MT,false>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*mat).columns() == (*vec).size(), "Invalid matrix and vector sizes" );

   return DMatSVecMultExpr<MT,VT>( *mat, *vec );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication of a symmetric row-major dense matrix
//        and a sparse vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side row-major dense matrix for the multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// symmetric row-major dense matrix and a sparse vector. It restructures the expression
// \f$ \vec{y}=A*\vec{x} \f$ to the expression \f$ \vec{y}=A^T*\vec{x} \f$.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename VT  // Type of the right-hand side sparse vector
        , EnableIf_t< IsSymmetric_v<MT> && !IsZero_v<VT> >* = nullptr >
inline decltype(auto)
   dmatsvecmult( const DenseMatrix<MT,false>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*mat).columns() == (*vec).size(), "Invalid matrix and vector sizes" );

   return trans( *mat ) * (*vec);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication of a row-major dense matrix and a zero
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side row-major dense matrix for the multiplication.
// \param vec The right-hand side zero vector for the multiplication.
// \return The resulting zero vector.
//
// This function implements the performance optimized treatment of the multiplication of a
// row-major dense matrix and a zero vector. It returns a zero vector.
*/
template< typename MT  // Type of the left-hand side dense matrix
        , typename VT  // Type of the right-hand side sparse vector
        , EnableIf_t< IsZero_v<VT> >* = nullptr >
inline decltype(auto)
   dmatsvecmult( const DenseMatrix<MT,false>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( vec );

   BLAZE_INTERNAL_ASSERT( (*mat).columns() == (*vec).size(), "Invalid matrix and vector sizes" );

   using ReturnType = const MultTrait_t< ResultType_t<MT>, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*mat).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a row-major dense matrix and a sparse
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side row-major dense matrix for the multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a row-major dense matrix and a sparse
// vector:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,rowMajor> A;
   blaze::CompressedVector<double,columnVector> x;
   blaze::DynamicVector<double,columnVector> y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved element types \a MT::ElementType and \a VT::ElementType. Both the
// dense matrix type \a MT and the sparse vector type \a VT as well as the two element types
// \a MT::ElementType and \a VT::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side sparse vector
inline decltype(auto)
   operator*( const DenseMatrix<MT,false>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATMATMULTEXPR_TYPE( MT );

   if( (*mat).columns() != (*vec).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix and vector sizes do not match" );
   }

   return dmatsvecmult( *mat, *vec );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, typename VT >
struct IsAligned< DMatSVecMultExpr<MT,VT> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
