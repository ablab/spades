//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatDVecMultExpr.h
//  \brief Header file for the transpose dense matrix/dense vector multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATDVECMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/gemv.h>
#include <blaze/math/blas/trmv.h>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/MatVecMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/VecScalarMultExpr.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsBLASCompatible.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyTriangular.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/BLAS.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/Complex.h>
#include <blaze/util/constraints/SameType.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsComplexDouble.h>
#include <blaze/util/typetraits/IsComplexFloat.h>
#include <blaze/util/typetraits/IsDouble.h>
#include <blaze/util/typetraits/IsFloat.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATDVECMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// The TDMatDVecMultExpr class represents the compile time expression for multiplications
// between column-major dense matrices and dense vectors.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side dense vector
class TDMatDVecMultExpr
   : public MatVecMultExpr< DenseVector< TDMatDVecMultExpr<MT,VT>, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using MRT = ResultType_t<MT>;     //!< Result type of the left-hand side dense matrix expression.
   using VRT = ResultType_t<VT>;     //!< Result type of the right-hand side dense vector expression.
   using MET = ElementType_t<MRT>;   //!< Element type of the left-hand side dense matrix expression.
   using VET = ElementType_t<VRT>;   //!< Element type of the right-hand side dense vector expression.
   using MCT = CompositeType_t<MT>;  //!< Composite type of the left-hand side dense matrix expression.
   using VCT = CompositeType_t<VT>;  //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateMatrix =
      ( ( IsComputation_v<MT> && IsSame_v<MET,VET> &&
          IsBLASCompatible_v<MET> ) || RequiresEvaluation_v<MT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense vector expression.
   static constexpr bool evaluateVector = ( IsComputation_v<VT> || RequiresEvaluation_v<VT> );
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either the matrix or the vector operand requires an intermediate evaluation, the
       variable will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   static constexpr bool UseSMPAssign_v = ( evaluateMatrix || evaluateVector );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a BLAS kernel,
       the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseBlasKernel_v =
      ( BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION &&
        IsContiguous_v<T1> && HasMutableDataAccess_v<T1> &&
        IsContiguous_v<T2> && HasConstDataAccess_v<T2> &&
        IsContiguous_v<T3> && HasConstDataAccess_v<T3> &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsBLASCompatible_v< ElementType_t<T1> > &&
        IsBLASCompatible_v< ElementType_t<T2> > &&
        IsBLASCompatible_v< ElementType_t<T3> > &&
        IsSame_v< ElementType_t<T1>, ElementType_t<T2> > &&
        IsSame_v< ElementType_t<T1>, ElementType_t<T3> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the matrix type and the two involved vector types are suited for a vectorized
       computation of the matrix/vector multiplication, the variable will be set to 1, otherwise
       it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseVectorizedDefaultKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3> > &&
        HasSIMDAdd_v< ElementType_t<T2>, ElementType_t<T3> > &&
        HasSIMDMult_v< ElementType_t<T2>, ElementType_t<T3> > );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TDMatDVecMultExpr instance.
   using This = TDMatDVecMultExpr<MT,VT>;

   //! Base type of this TDMatDVecMultExpr instance.
   using BaseType = MatVecMultExpr< DenseVector<This,false> >;

   using ResultType    = MultTrait_t<MRT,VRT>;         //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT>, const MT, const MT& >;

   //! Composite type of the right-hand side dense vector expression.
   using RightOperand = If_t< IsExpression_v<VT>, const VT, const VT& >;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateMatrix, const MRT, MCT >;

   //! Type for the assignment of the right-hand side dense vector operand.
   using RT = If_t< evaluateVector, const VRT, VCT >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( !IsDiagonal_v<MT> &&
        MT::simdEnabled && VT::simdEnabled &&
        HasSIMDAdd_v<MET,VET> &&
        HasSIMDMult_v<MET,VET> );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateMatrix && MT::smpAssignable && !evaluateVector && VT::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatDVecMultExpr class.
   //
   // \param mat The left-hand side matrix operand of the multiplication expression.
   // \param vec The right-hand side vector operand of the multiplication expression.
   */
   inline TDMatDVecMultExpr( const MT& mat, const VT& vec ) noexcept
      : mat_( mat )  // Left-hand side dense matrix of the multiplication expression
      , vec_( vec )  // Right-hand side dense vector of the multiplication expression
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
      else if( IsLower_v<MT> && ( index + 8UL < mat_.rows() ) )
      {
         const size_t n( IsStrictlyLower_v<MT> ? index : index+1UL );
         return subvector( row( mat_, index, unchecked ), 0UL, n, unchecked ) *
                subvector( vec_, 0UL, n, unchecked );
      }
      else if( IsUpper_v<MT> && ( index > 8UL ) )
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
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return mat_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side dense vector operand.
   //
   // \return The right-hand side dense vector operand.
   */
   inline RightOperand rightOperand() const noexcept {
      return vec_;
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
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
      return ( mat_.isAliased( alias ) || vec_.isAliased( alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return mat_.isAligned() && vec_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( !BLAZE_BLAS_MODE ||
               !BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION ||
               !BLAZE_BLAS_IS_PARALLEL ||
               ( IsComputation_v<MT> && !evaluateMatrix ) ||
               ( mat_.rows() * mat_.columns() < TDMATDVECMULT_THRESHOLD ) ) &&
             ( size() > SMP_TDMATDVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  mat_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand vec_;  //!< Right-hand side dense vector of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a dense vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }
      else if( rhs.mat_.columns() == 0UL ||
               ( IsStrictlyTriangular_v<MT> && rhs.mat_.columns() == 1UL ) ) {
         reset( *lhs );
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectAssignKernel( *lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense vectors (kernel selection)**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an assignment of a transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAssignKernel( y, A, x );
      else
         selectBlasAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose dense matrix-dense
   // vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      if( IsStrictlyLower_v<MT1> ) {
         reset( y[0] );
      }

      if( !IsUpper_v<MT1> )
      {
         for( size_t i=( IsStrictlyLower_v<MT1> ? 1UL : 0UL ); i<M; ++i ) {
            y[i] = A(i,0UL) * x[0UL];
         }
      }

      for( size_t j=( IsUpper_v<MT1> && !IsStrictlyUpper_v<MT1> ? 0UL : 1UL ); j<N; ++j )
      {
         if( IsDiagonal_v<MT1> )
         {
            y[j] = A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower_v<MT1> )
                                 ?( IsStrictlyLower_v<MT1> ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( IsStrictlyUpper_v<MT1> ? j-1UL : j )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
            if( IsUpper_v<MT1> ) {
               y[iend] = A(iend,j) * x[j];
            }
         }
      }

      if( IsStrictlyUpper_v<MT1> ) {
         reset( y[M-1UL] );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors (small matrices)****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a small transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (small matrices)*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the transpose dense
   // matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );
         SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jbegin) * x1 );
         SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jbegin) * x1 );
         SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jbegin) * x1 );
         SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
         y.store( i+SIMDSIZE*4UL, xmm5 );
         y.store( i+SIMDSIZE*5UL, xmm6 );
         y.store( i+SIMDSIZE*6UL, xmm7 );
         y.store( i+SIMDSIZE*7UL, xmm8 );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i         ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i         ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+SIMDSIZE, xmm2 );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( A.load(i,jbegin) * set( x[jbegin] ) );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            xmm1 += A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] = value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense vectors (large matrices)****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a large transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (large matrices)*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the transpose dense
   // matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      reset( y );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4 );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) + xmm5 );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) + xmm6 );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) + xmm7 );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) + xmm8 );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4 );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1 );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) + xmm2 );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      selectLargeAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors******************************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         assign( y, x );
         trmv( y, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
      }
      else {
         gemv( y, A, x, ET(1), ET(0) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-dense vector multiplication to a sparse vector
   //        (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief Addition assignment of a transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && rhs.mat_.rows() == 1UL ) ) {
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectAddAssignKernel( *lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense vectors (kernel selection)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an addition assignment of a transpose dense matrix-dense
   //        vector multiplication to a dense vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( y, A, x );
      else
         selectBlasAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         if( IsDiagonal_v<MT1> )
         {
            y[j] += A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower_v<MT1> )
                                 ?( IsStrictlyLower_v<MT1> ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( IsStrictlyUpper_v<MT1> ? j : j+1UL )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (small matrices)*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (small matrices)********************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );
         SIMDType xmm4( y.load(i+SIMDSIZE*3UL) );
         SIMDType xmm5( y.load(i+SIMDSIZE*4UL) );
         SIMDType xmm6( y.load(i+SIMDSIZE*5UL) );
         SIMDType xmm7( y.load(i+SIMDSIZE*6UL) );
         SIMDType xmm8( y.load(i+SIMDSIZE*7UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
         y.store( i+SIMDSIZE*4UL, xmm5 );
         y.store( i+SIMDSIZE*5UL, xmm6 );
         y.store( i+SIMDSIZE*6UL, xmm7 );
         y.store( i+SIMDSIZE*7UL, xmm8 );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );
         SIMDType xmm4( y.load(i+SIMDSIZE*3UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i         ) );
         SIMDType xmm2( y.load(i+SIMDSIZE) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 += A.load(i         ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+SIMDSIZE, xmm2 );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i) );

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 += A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] += value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (large matrices)*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (large matrices)********************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4 );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) + xmm5 );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) + xmm6 );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) + xmm7 );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) + xmm8 );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4 );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3 );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1 );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) + xmm2 );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      selectLargeAddAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors*********************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose matrix-vector multiplication
   //        (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         ResultType_t<VT1> tmp( serial( x ) );
         trmv( tmp, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
         addAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && rhs.mat_.rows() == 1UL ) ) {
         return;
      }

      LT A( serial( rhs.mat_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( rhs.vec_ ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      TDMatDVecMultExpr::selectSubAssignKernel( *lhs, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors (kernel selection)**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for a subtraction assignment of a transpose dense matrix-
   //        dense vector multiplication to a dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( y, A, x );
      else
         selectBlasSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline void selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         if( IsDiagonal_v<MT1> )
         {
            y[j] -= A(j,j) * x[j];
         }
         else
         {
            const size_t ibegin( ( IsLower_v<MT1> )
                                 ?( IsStrictlyLower_v<MT1> ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( IsStrictlyUpper_v<MT1> ? j : j+1UL )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] -= A(i    ,j) * x[j];
               y[i+1UL] -= A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] -= A(ipos,j) * x[j];
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (small matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a small transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (small matrices)*****************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a small transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for small
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );
         SIMDType xmm4( y.load(i+SIMDSIZE*3UL) );
         SIMDType xmm5( y.load(i+SIMDSIZE*4UL) );
         SIMDType xmm6( y.load(i+SIMDSIZE*5UL) );
         SIMDType xmm7( y.load(i+SIMDSIZE*6UL) );
         SIMDType xmm8( y.load(i+SIMDSIZE*7UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 -= A.load(i             ,j) * x1;
            xmm2 -= A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 -= A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 -= A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 -= A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 -= A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 -= A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 -= A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
         y.store( i+SIMDSIZE*4UL, xmm5 );
         y.store( i+SIMDSIZE*5UL, xmm6 );
         y.store( i+SIMDSIZE*6UL, xmm7 );
         y.store( i+SIMDSIZE*7UL, xmm8 );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );
         SIMDType xmm4( y.load(i+SIMDSIZE*3UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 -= A.load(i             ,j) * x1;
            xmm2 -= A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 -= A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 -= A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
         y.store( i+SIMDSIZE*3UL, xmm4 );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i             ) );
         SIMDType xmm2( y.load(i+SIMDSIZE    ) );
         SIMDType xmm3( y.load(i+SIMDSIZE*2UL) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 -= A.load(i             ,j) * x1;
            xmm2 -= A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 -= A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , xmm1 );
         y.store( i+SIMDSIZE    , xmm2 );
         y.store( i+SIMDSIZE*2UL, xmm3 );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i         ) );
         SIMDType xmm2( y.load(i+SIMDSIZE) );

         for( size_t j=jbegin; j<jend; ++j ) {
            const SIMDType x1( set( x[j] ) );
            xmm1 -= A.load(i         ,j) * x1;
            xmm2 -= A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , xmm1 );
         y.store( i+SIMDSIZE, xmm2 );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( y.load(i) );

         for( size_t j=jbegin; j<jend; ++j ) {
            xmm1 -= A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1 );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] -= value;
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (large matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a large transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      selectDefaultSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (large matrices)*****************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a large transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for large
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4 );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) - xmm5 );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) - xmm6 );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) - xmm7 );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) - xmm8 );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3 );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4 );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1 );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2 );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3 );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) - xmm1 );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) - xmm2 );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) - xmm1 );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] -= value;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      selectLargeSubAssignKernel( y, A, x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors******************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subtraction assignment of a transpose matrix-vector multiplication
   //        (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \return void
   //
   // This function performs the transpose dense matrix-dense vector multiplication based on the
   // according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2 >  // Type of the right-hand side vector operand
   static inline auto selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         ResultType_t<VT1> tmp( serial( x ) );
         trmv( tmp, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
         subAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(-1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Multiplication assignment of a transpose dense matrix-dense vector multiplication to
   //        a dense vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief Division assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}/=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief SMP assignment of a transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL ) {
         return;
      }
      else if( rhs.mat_.columns() == 0UL ||
               ( IsStrictlyTriangular_v<MT> && rhs.mat_.columns() == 1UL ) ) {
         reset( *lhs );
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      smpAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-dense vector multiplication to a sparse
   //        vector (\f$ \vec{y}=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-dense vector multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler
   // in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief SMP addition assignment of a transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}+=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && rhs.mat_.rows() == 1UL ) ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      smpAddAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}-=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      if( rhs.mat_.rows() == 0UL || rhs.mat_.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && rhs.mat_.rows() == 1UL ) ) {
         return;
      }

      LT A( rhs.mat_ );  // Evaluation of the left-hand side dense matrix operand
      RT x( rhs.vec_ );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.mat_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.mat_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == rhs.vec_.size()   , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size()     , "Invalid vector size"       );

      smpSubAssign( *lhs, A * x );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP multiplication assignment of a transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}*=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // transpose dense matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   /*!\brief SMP division assignment of a transpose dense matrix-dense vector multiplication to
   //        a dense vector (\f$ \vec{y}/=A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT1,false>& lhs, const TDMatDVecMultExpr& rhs )
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
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATVECMULTEXPR( MT, VT );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DVECSCALARMULTEXPR SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expression object for scaled transpose dense matrix-dense vector multiplications.
// \ingroup dense_vector_expression
//
// This specialization of the DVecScalarMultExpr class represents the compile time expression
// for scaled multiplications between a column-major dense matrix and a non-transpose dense
// vector.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT    // Type of the right-hand side dense vector
        , typename ST >  // Type of the side scalar value
class DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >
   : public VecScalarMultExpr< DenseVector< DVecScalarMultExpr< TDMatDVecMultExpr<MT,VT>, ST, false >, false > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using MVM = TDMatDVecMultExpr<MT,VT>;  //!< Type of the transpose dense matrix-dense vector multiplication expression.
   using RES = ResultType_t<MVM>;         //!< Result type of the dense matrix-dense vector multiplication expression.
   using MRT = ResultType_t<MT>;          //!< Result type of the left-hand side dense matrix expression.
   using VRT = ResultType_t<VT>;          //!< Result type of the right-hand side dense vector expression.
   using MET = ElementType_t<MRT>;        //!< Element type of the left-hand side dense matrix expression.
   using VET = ElementType_t<VRT>;        //!< Element type of the right-hand side dense vector expression.
   using MCT = CompositeType_t<MT>;       //!< Composite type of the left-hand side dense matrix expression.
   using VCT = CompositeType_t<VT>;       //!< Composite type of the right-hand side dense vector expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   static constexpr bool evaluateMatrix =
      ( ( IsComputation_v<MT> && IsSame_v<MET,VET> &&
          IsBLASCompatible_v<MET> ) || RequiresEvaluation_v<MT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense vector expression.
   static constexpr bool evaluateVector = ( IsComputation_v<VT> || RequiresEvaluation_v<VT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either the matrix or the vector operand requires an intermediate evaluation, the
       variable will be set to 1, otherwise it will be 0. */
   template< typename T1 >
   static constexpr bool UseSMPAssign_v = ( evaluateMatrix || evaluateVector );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the matrix type, the two involved vector types, and the scalar type are suited
       for a BLAS kernel, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   static constexpr bool UseBlasKernel_v =
      ( BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION &&
        IsContiguous_v<T1> && HasMutableDataAccess_v<T1> &&
        IsContiguous_v<T2> && HasConstDataAccess_v<T2> &&
        IsContiguous_v<T3> && HasConstDataAccess_v<T3> &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsBLASCompatible_v< ElementType_t<T1> > &&
        IsBLASCompatible_v< ElementType_t<T2> > &&
        IsBLASCompatible_v< ElementType_t<T3> > &&
        IsSame_v< ElementType_t<T1>, ElementType_t<T2> > &&
        IsSame_v< ElementType_t<T1>, ElementType_t<T3> > &&
        !( IsBuiltin_v< ElementType_t<T1> > && IsComplex_v<T4> ) );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the two involved vector types, the matrix type, and the scalar type are suited
       for a vectorized computation of the scaled vector/matrix multiplication, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   static constexpr bool UseVectorizedDefaultKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3>
                          , T4 > &&
        HasSIMDAdd_v< ElementType_t<T2>, ElementType_t<T3> > &&
        HasSIMDMult_v< ElementType_t<T2>, ElementType_t<T3> > );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DVecScalarMultExpr instance.
   using This = DVecScalarMultExpr<MVM,ST,false>;

   //! Base type of this DVecScalarMultExpr instance.
   using BaseType = VecScalarMultExpr< DenseVector<This,false> >;

   using ResultType    = MultTrait_t<RES,ST>;          //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense vector expression.
   using LeftOperand = const TDMatDVecMultExpr<MT,VT>;

   //! Composite type of the right-hand side scalar value.
   using RightOperand = ST;

   //! Type for the assignment of the dense matrix operand of the left-hand side expression.
   using LT = If_t< evaluateMatrix, const MRT, MCT >;

   //! Type for the assignment of the dense vector operand of the left-hand side expression.
   using RT = If_t< evaluateVector, const VRT, VCT >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( !IsDiagonal_v<MT> &&
        MT::simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v<MET,VET,ST> &&
        HasSIMDAdd_v<MET,VET> &&
        HasSIMDMult_v<MET,VET> );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateMatrix && MT::smpAssignable && !evaluateVector && VT::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DVecScalarMultExpr class.
   //
   // \param vector The left-hand side dense vector of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   inline DVecScalarMultExpr( const MVM& vector, ST scalar )
      : vector_( vector )  // Left-hand side dense vector of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
   {}
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < vector_.size(), "Invalid vector access index" );
      return vector_[index] * scalar_;
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
      if( index >= vector_.size() ) {
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
   inline size_t size() const {
      return vector_.size();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense vector operand.
   //
   // \return The left-hand side dense vector operand.
   */
   inline LeftOperand leftOperand() const {
      return vector_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side scalar operand.
   //
   // \return The right-hand side scalar operand.
   */
   inline RightOperand rightOperand() const {
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
   inline bool canAlias( const T* alias ) const {
      return vector_.canAlias( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const {
      return vector_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return vector_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      LeftOperand_t<MVM> A( vector_.leftOperand() );
      return ( !BLAZE_BLAS_MODE ||
               !BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION ||
               !BLAZE_BLAS_IS_PARALLEL ||
               ( IsComputation_v<MT> && !evaluateMatrix ) ||
               ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) ) &&
             ( size() > SMP_TDMATDVECMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  vector_;  //!< Left-hand side dense vector of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense vectors*****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a dense
   //        vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void assign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ||
               ( IsStrictlyTriangular_v<MT> && left.columns() == 1UL ) ) {
         reset( *lhs );
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectAssignKernel( *lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Assignment to dense vectors (kernel selection)**********************************************
   /*!\brief Selection of the kernel for an assignment of a scaled transpose dense matrix-dense
   //        vector multiplication to a dense vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAssignKernel( y, A, x, scalar );
      else
         selectBlasAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors*********************************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the default assignment kernel for the scaled transpose dense
   // matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      if( IsStrictlyLower_v<MT1> ) {
         reset( y[0] );
      }

      if( !IsUpper_v<MT1> )
      {
         for( size_t i=( IsStrictlyLower_v<MT1> ? 1UL : 0UL ); i<M; ++i ) {
            y[i] = A(i,0UL) * x[0UL];
         }
      }

      for( size_t j=( IsUpper_v<MT1> && !IsStrictlyUpper_v<MT1> ? 0UL : 1UL ); j<N; ++j )
      {
         if( IsDiagonal_v<MT1> )
         {
            y[j] = A(j,j) * x[j] * scalar;
         }
         else
         {
            const size_t ibegin( ( IsLower_v<MT1> )
                                 ?( IsStrictlyLower_v<MT1> ? j+1UL : j )
                                 :( 0UL ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( IsStrictlyUpper_v<MT1> ? j-1UL : j )
                               :( M ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               y[i    ] += A(i    ,j) * x[j];
               y[i+1UL] += A(i+1UL,j) * x[j];
            }
            if( ipos < iend ) {
               y[ipos] += A(ipos,j) * x[j];
            }
            if( IsUpper_v<MT1> ) {
               y[iend] = A(iend,j) * x[j];
            }
         }
      }

      if( IsStrictlyUpper_v<MT1> ) {
         reset( y[M-1UL] );
      }

      if( !IsDiagonal_v<MT1> )
      {
         const size_t iend( IsStrictlyUpper_v<MT1> ? M-1UL : M );
         for( size_t i=( IsStrictlyLower_v<MT1> ? 1UL : 0UL ); i<iend; ++i ) {
            y[i] *= scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors (small matrices)****************************************
   /*!\brief Default assignment of a small scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (small matrices)*****************************
   /*!\brief Vectorized default assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );
         SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jbegin) * x1 );
         SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jbegin) * x1 );
         SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jbegin) * x1 );
         SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+SIMDSIZE    , xmm2*factor );
         y.store( i+SIMDSIZE*2UL, xmm3*factor );
         y.store( i+SIMDSIZE*3UL, xmm4*factor );
         y.store( i+SIMDSIZE*4UL, xmm5*factor );
         y.store( i+SIMDSIZE*5UL, xmm6*factor );
         y.store( i+SIMDSIZE*6UL, xmm7*factor );
         y.store( i+SIMDSIZE*7UL, xmm8*factor );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+SIMDSIZE    , xmm2*factor );
         y.store( i+SIMDSIZE*2UL, xmm3*factor );
         y.store( i+SIMDSIZE*3UL, xmm4*factor );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , xmm1*factor );
         y.store( i+SIMDSIZE    , xmm2*factor );
         y.store( i+SIMDSIZE*2UL, xmm3*factor );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i         ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i         ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , xmm1*factor );
         y.store( i+SIMDSIZE, xmm2*factor );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( A.load(i,jbegin) * set( x[jbegin] ) );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            xmm1 += A.load(i,j) * set( x[j] );
         }

         y.store( i, xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] = value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense vectors (large matrices)****************************************
   /*!\brief Default assignment of a large scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense vectors (large matrices)*****************************
   /*!\brief Vectorized default assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor
   // \return void
   //
   // This function implements the vectorized default assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication. This kernel is optimized for large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      const SIMDType factor( set( scalar ) );

      reset( y );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) + xmm5*factor );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) + xmm6*factor );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) + xmm7*factor );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) + xmm8*factor );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1*factor );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) + xmm2*factor );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors (default)********************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectLargeAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense vectors******************************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-dense vector multiplication
   //        (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         assign( y, scalar * x );
         trmv( y, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
      }
      else {
         gemv( y, A, x, ET(scalar), ET(0) );
      }
   }
#endif
   //**********************************************************************************************

   //**Assignment to sparse vectors****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-dense vector multiplication to a sparse
   //        vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-dense vector multiplication expression to a sparse vector.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline void assign( SparseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors********************************************************
   /*!\brief Addition assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void addAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && left.rows() == 1UL ) ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectAddAssignKernel( *lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Addition assignment to dense vectors (kernel selection)*************************************
   /*!\brief Selection of the kernel for an addition assignment of a scaled transpose dense
   //        matrix-dense vector multiplication to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( y, A, x, scalar );
      else
         selectBlasAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors************************************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.addAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (small matrices)*******************************
   /*!\brief Default addition assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (small matrices)********************
   /*!\brief Vectorized default addition assignment of a small scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for small
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );
         SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jbegin) * x1 );
         SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jbegin) * x1 );
         SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jbegin) * x1 );
         SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
         y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
         y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) + xmm5*factor );
         y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) + xmm6*factor );
         y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) + xmm7*factor );
         y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) + xmm8*factor );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
         y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , y.load(i             ) + xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i         ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i         ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , y.load(i         ) + xmm1*factor );
         y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) + xmm2*factor );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( A.load(i,jbegin) * set( x[jbegin] ) );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            xmm1 += A.load(i,j) * set( x[j] );
         }

         y.store( i, y.load(i) + xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] += value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense vectors (large matrices)*******************************
   /*!\brief Default addition assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense vectors (large matrices)********************
   /*!\brief Vectorized default addition assignment of a large scaled transpose dense matrix-dense
   //        vector multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment kernel for the scaled
   // transpose dense matrix-dense vector multiplication. This kernel is optimized for large
   // matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      const SIMDType factor( set( scalar ) );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) + xmm5*factor );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) + xmm6*factor );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) + xmm7*factor );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) + xmm8*factor );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) + xmm4*factor );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) + xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) + xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) + xmm3*factor );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) + xmm1*factor );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) + xmm2*factor );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) + xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] += value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors (default)***********************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectLargeAddAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense vectors*********************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAddAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         ResultType_t<VT1> tmp( serial( scalar * x ) );
         trmv( tmp, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
         addAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Addition assignment to sparse vectors*******************************************************
   // No special implementation for the addition assignment to sparse vectors.
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors*****************************************************
   /*!\brief Subtraction assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a scaled
   // dense transpose matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void subAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && left.rows() == 1UL ) ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT x( serial( right ) );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      DVecScalarMultExpr::selectSubAssignKernel( *lhs, A, x, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Subtraction assignment to dense vectors (kernel selection)**********************************
   /*!\brief Selection of the kernel for a subtraction assignment of a scaled transpose dense
   //        matrix-dense vector multiplication to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT1> ) ||
          ( IsComputation_v<MT> && !evaluateMatrix ) ||
          ( A.rows() * A.columns() < TDMATDVECMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( y, A, x, scalar );
      else
         selectBlasSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors*********************************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the scaled transpose
   // dense matrix-dense vector multiplication.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectDefaultSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
   {
      y.subAssign( A * x * scalar );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (small matrices)****************************
   /*!\brief Default subtraction assignment of a small scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (small matrices)*****************
   /*!\brief Vectorized default subtraction assignment of a small scaled transpose dense matrix-
   //        dense vector multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // scaled transpose dense matrix-dense vector multiplication. This kernel is optimized for
   // small matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*8UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );
         SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jbegin) * x1 );
         SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jbegin) * x1 );
         SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jbegin) * x1 );
         SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
            xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
            xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
            xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
            xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
         y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4*factor );
         y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) - xmm5*factor );
         y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) - xmm6*factor );
         y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) - xmm7*factor );
         y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) - xmm8*factor );
      }

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*4UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );
         SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
            xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
         y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4*factor );
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*3UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i             ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE    ,jbegin) * x1 );
         SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i             ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
            xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
         }

         y.store( i             , y.load(i             ) - xmm1*factor );
         y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
         y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE*2UL, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType x1( set( x[jbegin] ) );
         SIMDType xmm1( A.load(i         ,jbegin) * x1 );
         SIMDType xmm2( A.load(i+SIMDSIZE,jbegin) * x1 );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            x1 = set( x[j] );
            xmm1 += A.load(i         ,j) * x1;
            xmm2 += A.load(i+SIMDSIZE,j) * x1;
         }

         y.store( i         , y.load(i         ) - xmm1*factor );
         y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) - xmm2*factor );
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )
                              ?( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )
                            ?( min( i+SIMDSIZE, N ) - ( IsStrictlyLower_v<MT1> ? 1UL : 0UL ) )
                            :( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         SIMDType xmm1( A.load(i,jbegin) * set( x[jbegin] ) );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            xmm1 += A.load(i,j) * set( x[j] );
         }

         y.store( i, y.load(i) - xmm1*factor );
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jbegin( ( IsUpper_v<MT1> )?( i ):( 0UL ) );
         const size_t jend( ( IsLower_v<MT1> )?( min( i+1UL, N ) ):( N ) );
         BLAZE_INTERNAL_ASSERT( jbegin < jend, "Invalid loop indices detected" );

         ElementType value( A(i,jbegin) * x[jbegin] );

         for( size_t j=jbegin+1UL; j<jend; ++j ) {
            value += A(i,j) * x[j];
         }

         y[i] -= value * scalar;
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense vectors (large matrices)****************************
   /*!\brief Default subtraction assignment of a large scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectDefaultSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense vectors (large matrices)*****************
   /*!\brief Vectorized default subtraction assignment of a large scaled transpose dense matrix-
   //        dense vector multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment kernel for the
   // scaled transpose dense matrix-dense vector multiplication. This kernel is optimized for
   // large matrices.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<VT1,MT1,VT2,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT1> || !IsPadded_v<VT1> );

      const size_t M( A.rows()    );
      const size_t N( A.columns() );

      const size_t iblock( 32768UL / sizeof( ElementType ) );
      const size_t jblock( ( N < iblock )?( 8UL ):( 4UL ) );

      BLAZE_INTERNAL_ASSERT( ( iblock % SIMDSIZE ) == 0UL, "Invalid block size detected" );

      const SIMDType factor( set( scalar ) );

      for( size_t ii=0U; ii<M; ii+=iblock ) {
         for( size_t jj=0UL; jj<N; jj+=jblock )
         {
            const size_t jend( min( jj+jblock, N ) );
            const size_t itmp( min( ii+iblock, M ) );
            const size_t iend( ( IsUpper_v<MT1> )
                               ?( min( itmp, ( IsStrictlyUpper_v<MT1> ? jend-1UL : jend ) ) )
                               :( itmp ) );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ( IsLower_v<MT1> )
                      ?( max( ii, prevMultiple( ( IsStrictlyLower_v<MT1> ? jj+1UL : jj ), SIMDSIZE ) ) )
                      :( ii ) );

            for( ; (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,jj) * x1 );
               SIMDType xmm6( A.load(i+SIMDSIZE*5UL,jj) * x1 );
               SIMDType xmm7( A.load(i+SIMDSIZE*6UL,jj) * x1 );
               SIMDType xmm8( A.load(i+SIMDSIZE*7UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,j) * x1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,j) * x1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,j) * x1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4*factor );
               y.store( i+SIMDSIZE*4UL, y.load(i+SIMDSIZE*4UL) - xmm5*factor );
               y.store( i+SIMDSIZE*5UL, y.load(i+SIMDSIZE*5UL) - xmm6*factor );
               y.store( i+SIMDSIZE*6UL, y.load(i+SIMDSIZE*6UL) - xmm7*factor );
               y.store( i+SIMDSIZE*7UL, y.load(i+SIMDSIZE*7UL) - xmm8*factor );
            }

            for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
               y.store( i+SIMDSIZE*3UL, y.load(i+SIMDSIZE*3UL) - xmm4*factor );
            }

            for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i             ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,jj) * x1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i             ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE    ,j) * x1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,j) * x1;
               }

               y.store( i             , y.load(i             ) - xmm1*factor );
               y.store( i+SIMDSIZE    , y.load(i+SIMDSIZE    ) - xmm2*factor );
               y.store( i+SIMDSIZE*2UL, y.load(i+SIMDSIZE*2UL) - xmm3*factor );
            }

            for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
            {
               SIMDType x1( set( x[jj] ) );
               SIMDType xmm1( A.load(i         ,jj) * x1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,jj) * x1 );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  x1 = set( x[j] );
                  xmm1 += A.load(i         ,j) * x1;
                  xmm2 += A.load(i+SIMDSIZE,j) * x1;
               }

               y.store( i         , y.load(i         ) - xmm1*factor );
               y.store( i+SIMDSIZE, y.load(i+SIMDSIZE) - xmm2*factor );
            }

            for( ; i<ipos; i+=SIMDSIZE )
            {
               SIMDType xmm1( A.load(i,jj) * set( x[jj] ) );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  xmm1 += A.load(i,j) * set( x[j] );
               }

               y.store( i, y.load(i) - xmm1*factor );
            }

            for( ; remainder && i<iend; ++i )
            {
               ElementType value( A(i,jj) * x[jj] );

               for( size_t j=jj+1UL; j<jend; ++j ) {
                  value += A(i,j) * x[j];
               }

               y[i] -= value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors (default)********************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      selectLargeSubAssignKernel( y, A, x, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense vectors******************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_VECTOR_MULTIPLICATION
   /*!\brief BLAS-based subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param y The target left-hand side dense vector.
   // \param A The left-hand side dense matrix operand.
   // \param x The right-hand side dense vector operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-dense vector multiplication based
   // on the according BLAS functionality.
   */
   template< typename VT1    // Type of the left-hand side target vector
           , typename MT1    // Type of the left-hand side matrix operand
           , typename VT2    // Type of the right-hand side vector operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasSubAssignKernel( VT1& y, const MT1& A, const VT2& x, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<VT1,MT1,VT2,ST2> >
   {
      using ET = ElementType_t<VT1>;

      if( IsTriangular_v<MT1> ) {
         ResultType_t<VT1> tmp( serial( scalar * x ) );
         trmv( tmp, A, ( IsLower_v<MT1> )?( CblasLower ):( CblasUpper ) );
         subAssign( y, tmp );
      }
      else {
         gemv( y, A, x, ET(-scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Subtraction assignment to sparse vectors****************************************************
   // No special implementation for the subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**Multiplication assignment to dense vectors**************************************************
   /*!\brief Multiplication assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}*=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized multiplication assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void multAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   /*!\brief Division assignment of a scaled transpose dense matrix-dense vector multiplication to
   //        a dense vector (\f$ \vec{y}/=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized division assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline void divAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( ResultType );
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
   /*!\brief SMP assignment of a scaled transpose dense matrix-dense vector multiplication to a
   //        dense vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a dense vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ||
               ( IsStrictlyTriangular_v<MT> && left.columns() == 1UL ) ) {
         reset( *lhs );
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      smpAssign( *lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP assignment to sparse vectors************************************************************
   /*!\brief SMP assignment of a scaled transpose dense matrix-dense vector multiplication to a
   //        sparse vector (\f$ \vec{y}=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side sparse vector.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-dense vector multiplication expression to a sparse vector. Due to the explicit
   // application of the SFINAE principle, this function can only be selected by the compiler in
   // case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target sparse vector
   friend inline auto smpAssign( SparseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
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
   //**********************************************************************************************

   //**SMP addition assignment to dense vectors****************************************************
   /*!\brief SMP addition assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}+=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpAddAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && left.rows() == 1UL ) ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      smpAddAssign( *lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP addition assignment to sparse vectors***************************************************
   // No special implementation for the SMP addition assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense vectors*************************************************
   /*!\brief SMP subtraction assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}-=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a scaled
   // dense transpose matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpSubAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
      -> EnableIf_t< UseSMPAssign_v<VT1> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).size() == rhs.size(), "Invalid vector sizes" );

      LeftOperand_t<MVM>  left ( rhs.vector_.leftOperand()  );
      RightOperand_t<MVM> right( rhs.vector_.rightOperand() );

      if( left.rows() == 0UL || left.columns() == 0UL ||
          ( IsStrictlyTriangular_v<MT> && left.rows() == 1UL ) ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT x( right );  // Evaluation of the right-hand side dense vector operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( x.size()    == right.size()  , "Invalid vector size"       );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).size() , "Invalid vector size"       );

      smpSubAssign( *lhs, A * x * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse vectors************************************************
   // No special implementation for the SMP subtraction assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP multiplication assignment to dense vectors**********************************************
   /*!\brief SMP multiplication assignment of a scaled transpose dense matrix-dense vector
   //        multiplication to a dense vector (\f$ \vec{y}*=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression to be multiplied.
   // \return void
   //
   // This function implements the performance optimized SMP multiplication assignment of a
   // scaled transpose dense matrix-dense vector multiplication expression to a dense vector.
   // Due to the explicit application of the SFINAE principle, this function can only be
   // selected by the compiler in case the expression specific parallel evaluation strategy
   // is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpMultAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
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
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse vectors*********************************************
   // No special implementation for the SMP multiplication assignment to sparse vectors.
   //**********************************************************************************************

   //**SMP division assignment to dense vectors****************************************************
   /*!\brief SMP division assignment of a scaled transpose dense matrix-dense vector multiplication
   //        to a dense vector (\f$ \vec{y}/=s*A*\vec{x} \f$).
   // \ingroup dense_vector
   //
   // \param lhs The target left-hand side dense vector.
   // \param rhs The right-hand side scaled multiplication expression divisor.
   // \return void
   //
   // This function implements the performance optimized SMP division assignment of a scaled
   // transpose dense matrix-dense vector multiplication expression to a dense vector. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case the expression specific parallel evaluation strategy is selected.
   */
   template< typename VT1 >  // Type of the target dense vector
   friend inline auto smpDivAssign( DenseVector<VT1,false>& lhs, const DVecScalarMultExpr& rhs )
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
   //**********************************************************************************************

   //**SMP division assignment to sparse vectors***************************************************
   // No special implementation for the SMP division assignment to sparse vectors.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( MVM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( MVM );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( ST );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TYPE( ST, RightOperand );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a dense
//        vector (\f$ \vec{y}=A*\vec{x} \f$).
// \ingroup dense_vector
//
// \param mat The left-hand side column-major dense matrix for the multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
// \exception std::invalid_argument Matrix and vector sizes do not match.
//
// This operator represents the multiplication between a column-major dense matrix and a dense
// vector:

   \code
   using blaze::columnMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<double,columnMajor> A;
   blaze::DynamicVector<double,columnVector> x, y;
   // ... Resizing and initialization
   y = A * x;
   \endcode

// The operator returns an expression representing a dense vector of the higher-order element
// type of the two involved element types \a MT::ElementType and \a VT::ElementType. Both the
// dense matrix type \a MT and the dense vector type \a VT as well as the two element types
// \a MT::ElementType and \a VT::ElementType have to be supported by the MultTrait class
// template.\n
// In case the current size of the vector \a vec doesn't match the current number of columns
// of the matrix \a mat, a \a std::invalid_argument is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , typename VT >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator*( const DenseMatrix<MT,true>& mat, const DenseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_MATMATMULTEXPR_TYPE( MT );

   if( (*mat).columns() != (*vec).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix and vector sizes do not match" );
   }

   using ReturnType = const TDMatDVecMultExpr<MT,VT>;
   return ReturnType( *mat, *vec );
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
struct IsAligned< TDMatDVecMultExpr<MT,VT> >
   : public BoolConstant< IsAligned_v<MT> && IsAligned_v<VT> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
