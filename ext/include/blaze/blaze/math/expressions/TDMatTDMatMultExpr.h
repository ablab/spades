//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatTDMatMultExpr.h
//  \brief Header file for the transpose dense matrix/transpose dense matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATTDMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATTDMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/blas/gemm.h>
#include <blaze/math/blas/trmm.h>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/dense/MMM.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/functors/DeclDiag.h>
#include <blaze/math/functors/DeclHerm.h>
#include <blaze/math/functors/DeclLow.h>
#include <blaze/math/functors/DeclSym.h>
#include <blaze/math/functors/DeclUpp.h>
#include <blaze/math/functors/Noop.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsBLASCompatible.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyTriangular.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/BLAS.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/Debugging.h>
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
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATTDMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-transpose dense matrix multiplications.
// \ingroup dense_matrix_expression
//
// The TDMatTDMatMultExpr class represents the compile time expression for multiplications between
// two column-major dense matrices.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
class TDMatTDMatMultExpr
   : public MatMatMultExpr< DenseMatrix< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side dense matrix expression.
   using ET1 = ElementType_t<RT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<RT2>;    //!< Element type of the right-hand side dense matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateLeft = ( IsComputation_v<MT1> || RequiresEvaluation_v<MT1> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   static constexpr bool evaluateRight = ( IsComputation_v<MT2> || RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool SYM  = ( SF && !( HF || LF || UF )    );  //!< Flag for symmetric matrices.
   static constexpr bool HERM = ( HF && !( LF || UF )          );  //!< Flag for Hermitian matrices.
   static constexpr bool LOW  = ( LF || ( ( SF || HF ) && UF ) );  //!< Flag for lower matrices.
   static constexpr bool UPP  = ( UF || ( ( SF || HF ) && LF ) );  //!< Flag for upper matrices.
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the optimal evaluation strategy.
       In case the target matrix is row-major and either of the two matrix operands is symmetric,
       the variable is set to 1 and an optimized evaluation strategy is selected. Otherwise the
       variable is set to 0 and the default strategy is chosen. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool CanExploitSymmetry_v =
      ( IsRowMajorMatrix_v<T1> && ( IsSymmetric_v<T2> || IsSymmetric_v<T3> ) );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either of the two matrix operands requires an intermediate evaluation, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool IsEvaluationRequired_v =
      ( ( evaluateLeft || evaluateRight ) && CanExploitSymmetry_v<T1,T2,T3> );
   /*! \endcond */
   //**********************************************************************************************

      //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the types of all three involved matrices are suited for a BLAS kernel, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseBlasKernel_v =
      ( BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION &&
        !SYM && !HERM && !LOW && !UPP &&
        IsContiguous_v<T1> && HasMutableDataAccess_v<T1> &&
        IsContiguous_v<T2> && HasConstDataAccess_v<T2> &&
        IsContiguous_v<T3> && HasConstDataAccess_v<T3> &&
        !IsDiagonal_v<T2> && !IsDiagonal_v<T3> &&
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
   /*! In case all three involved data types are suited for a vectorized computation of the
       matrix multiplication, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseVectorizedDefaultKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> && !IsDiagonal_v<T3> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3> > &&
        HasSIMDAdd_v< ElementType_t<T2>, ElementType_t<T3> > &&
        HasSIMDMult_v< ElementType_t<T2>, ElementType_t<T3> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Type of the functor for forwarding an expression to another assign kernel.
   /*! In case a temporary matrix needs to be created, this functor is used to forward the
       resulting expression to another assign kernel. */
   using ForwardFunctor = If_t< HERM
                              , DeclHerm
                              , If_t< SYM
                                    , DeclSym
                                    , If_t< LOW
                                          , If_t< UPP
                                                , DeclDiag
                                                , DeclLow >
                                          , If_t< UPP
                                                , DeclUpp
                                                , Noop > > > >;
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TDMatTDMatMultExpr instance.
   using This = TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>;

   //! Base type of this TDMatTDMatMultExpr instance.
   using BaseType = MatMatMultExpr< DenseMatrix<This,true> >;

   //! Result type for expression template evaluations.
   using ResultType = typename If_t< HERM
                                   , DeclHermTrait< MultTrait_t<RT1,RT2> >
                                   , If_t< SYM
                                         , DeclSymTrait< MultTrait_t<RT1,RT2> >
                                         , If_t< LOW
                                               , If_t< UPP
                                                     , DeclDiagTrait< MultTrait_t<RT1,RT2> >
                                                     , DeclLowTrait< MultTrait_t<RT1,RT2> > >
                                               , If_t< UPP
                                                     , DeclUppTrait< MultTrait_t<RT1,RT2> >
                                                     , MultTrait<RT1,RT2> > > > >::Type;

   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side dense matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateLeft, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side dense matrix operand.
   using RT = If_t< evaluateRight, const RT2, CT2 >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( !IsDiagonal_v<MT1> &&
        MT1::simdEnabled && MT2::simdEnabled &&
        HasSIMDAdd_v<ET1,ET2> &&
        HasSIMDMult_v<ET1,ET2> );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateLeft  && MT1::smpAssignable && !evaluateRight && MT2::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TDMatTDMatMultExpr class.
   //
   // \param lhs The left-hand side operand of the multiplication expression.
   // \param rhs The right-hand side operand of the multiplication expression.
   */
   inline TDMatTDMatMultExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the multiplication expression
      , rhs_( rhs )  // Right-hand side dense matrix of the multiplication expression
   {
      BLAZE_INTERNAL_ASSERT( lhs.columns() == rhs.rows(), "Invalid matrix sizes" );
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
      BLAZE_INTERNAL_ASSERT( j < rhs_.columns(), "Invalid column access index" );

      if( IsDiagonal_v<MT1> ) {
         return lhs_(i,i) * rhs_(i,j);
      }
      else if( IsDiagonal_v<MT2> ) {
         return lhs_(i,j) * rhs_(j,j);
      }
      else if( IsTriangular_v<MT1> || IsTriangular_v<MT2> ) {
         const size_t begin( ( IsUpper_v<MT1> )
                             ?( ( IsLower_v<MT2> )
                                ?( max( ( IsStrictlyUpper_v<MT1> ? i+1UL : i )
                                      , ( IsStrictlyLower_v<MT2> ? j+1UL : j ) ) )
                                :( IsStrictlyUpper_v<MT1> ? i+1UL : i ) )
                             :( ( IsLower_v<MT2> )
                                ?( IsStrictlyLower_v<MT2> ? j+1UL : j )
                                :( 0UL ) ) );
         const size_t end( ( IsLower_v<MT1> )
                           ?( ( IsUpper_v<MT2> )
                              ?( min( ( IsStrictlyLower_v<MT1> ? i : i+1UL )
                                    , ( IsStrictlyUpper_v<MT2> ? j : j+1UL ) ) )
                              :( IsStrictlyLower_v<MT1> ? i : i+1UL ) )
                           :( ( IsUpper_v<MT2> )
                              ?( IsStrictlyUpper_v<MT2> ? j : j+1UL )
                              :( lhs_.columns() ) ) );

         if( begin >= end ) return ElementType();

         const size_t n( end - begin );

         return subvector( row( lhs_, i, unchecked ), begin, n, unchecked ) *
                subvector( column( rhs_, j, unchecked ), begin, n, unchecked );
      }
      else {
         return row( lhs_, i, unchecked ) * column( rhs_, j, unchecked );
      }
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
      if( j >= rhs_.columns() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
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
      return rhs_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side transpose dense matrix operand.
   //
   // \return The left-hand side transpose dense matrix operand.
   */
   inline LeftOperand leftOperand() const noexcept {
      return lhs_;
   }
   //**********************************************************************************************

   //**Right operand access************************************************************************
   /*!\brief Returns the right-hand side transpose dense matrix operand.
   //
   // \return The right-hand side transpose dense matrix operand.
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
      return ( lhs_.isAliased( alias ) || rhs_.isAliased( alias ) );
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

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return lhs_.isAligned() && rhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( !BLAZE_BLAS_MODE ||
               !BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION ||
               !BLAZE_BLAS_IS_PARALLEL ||
               ( rows() * columns() < TDMATTDMATMULT_THRESHOLD ) ) &&
             ( rows() * columns() >= SMP_TDMATTDMATMULT_THRESHOLD ) &&
             !IsDiagonal_v<MT1> && !IsDiagonal_v<MT2>;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side dense matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-transpose dense matrix multiplication to a
   //        dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL ) {
         return;
      }
      else if( rhs.lhs_.columns() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTDMatMultExpr::selectAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to dense matrices (kernel selection)*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an assignment of a transpose dense matrix-transpose
   //        dense matrix multiplication to a dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallAssignKernel( C, A, B );
      else
         selectBlasAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (general/general)**************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a general transpose dense matrix-general transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a general transpose dense matrix-general
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( SYM || HERM || LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t kbegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t kend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( K ) );
         BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

         if( IsStrictlyTriangular_v<MT5> && kbegin == kend ) {
            for( size_t i=0UL; i<M; ++i ) {
               reset( C(i,j) );
            }
            continue;
         }

         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j,kbegin+1UL) : kbegin+1UL )
                                    :( LOW ? max(j,kbegin) : kbegin ) )
                                 :( LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,kbegin) : kbegin )
                                  :( UPP ? min(j,kbegin)+1UL : kbegin+1UL ) )
                               :( UPP ? j+1UL : M ) );

            if( ( IsLower_v<MT4> && IsLower_v<MT5> ) || LOW ) {
               for( size_t i=0UL; i<ibegin; ++i ) {
                  reset( C(i,j) );
               }
            }
            else if( IsStrictlyLower_v<MT4> ) {
               reset( C(0UL,j) );
            }
            for( size_t i=ibegin; i<iend; ++i ) {
               C(i,j) = A(i,kbegin) * B(kbegin,j);
            }
            if( ( IsUpper_v<MT4> && IsUpper_v<MT5> ) || UPP ) {
               for( size_t i=iend; i<M; ++i ) {
                  reset( C(i,j) );
               }
            }
            else if( IsStrictlyUpper_v<MT4> ) {
               reset( C(M-1UL,j) );
            }
         }

         for( size_t k=kbegin+1UL; k<kend; ++k )
         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( SYM || HERM || LOW ? max( j, k+1UL ) : k+1UL )
                                    :( SYM || HERM || LOW ? max( j, k ) : k ) )
                                 :( SYM || HERM || LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,k-1UL) : k-1UL )
                                  :( UPP ? min(j+1UL,k) : k ) )
                               :( UPP ? j+1UL : M ) );

            if( ( SYM || HERM || LOW || UPP ) && ( ibegin > iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            for( size_t i=ibegin; i<iend; ++i ) {
               C(i,j) += A(i,k) * B(k,j);
            }
            if( IsUpper_v<MT4> ) {
               C(iend,j) = A(iend,k) * B(k,j);
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=1UL; j<N; ++j ) {
            for( size_t i=0UL; i<j; ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (general/diagonal)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a general transpose dense matrix-diagonal transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a general transpose dense matrix-diagonal
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         if( IsLower_v<MT4> ) {
            for( size_t i=0UL; i<ibegin; ++i ) {
               reset( C(i,j) );
            }
         }
         for( size_t i=ibegin; i<iend; ++i ) {
            C(i,j) = A(i,j) * B(j,j);
         }
         if( IsUpper_v<MT4> ) {
            for( size_t i=iend; i<M; ++i ) {
               reset( C(i,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (diagonal/general)*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a diagonal transpose dense matrix-general transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a diagonal transpose dense matrix-general
   // transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         if( IsLower_v<MT4> ) {
            for( size_t i=0UL; i<ibegin; ++i ) {
               reset( C(i,j) );
            }
         }
         for( size_t i=ibegin; i<iend; ++i ) {
            C(i,j) = A(i,i) * B(i,j);
         }
         if( IsUpper_v<MT4> ) {
            for( size_t i=iend; i<M; ++i ) {
               reset( C(i,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (diagonal/diagonal)************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a diagonal transpose dense matrix-diagonal transpose dense
   //        matrix multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default assignment of a diagonal transpose dense matrix-
   // diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      reset( C );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) = A(i,i) * B(i,i);
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (small matrices)***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a small transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to row-major dense matrices (small matrices)******************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a small transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default assignment of a transpose dense matrix-
   // transpose dense matrix multiplication expression to a row-major dense matrix. This kernel
   // is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         assign( C, fwd( A * tmp ) );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         assign( C, fwd( tmp * B ) );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         assign( C, fwd( A * tmp ) );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         assign( C, fwd( tmp * B ) );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to column-major dense matrices (small matrices)***************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a small transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default assignment of a transpose dense matrix-
   // transpose dense matrix multiplication expression to a column-major dense matrix. This
   // kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( SYM || HERM || LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !SYM && !HERM && !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               size_t k( kbegin );

               if( k < kend )
               {
                  SIMDType b1( set( B(k,j) ) );
                  SIMDType xmm1( A.load(i             ,k) * b1 );
                  SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
                  SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
                  SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
                  SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );
                  SIMDType xmm6( A.load(i+SIMDSIZE*5UL,k) * b1 );
                  SIMDType xmm7( A.load(i+SIMDSIZE*6UL,k) * b1 );
                  SIMDType xmm8( A.load(i+SIMDSIZE*7UL,k) * b1 );

                  for( ++k; k<kend; ++k ) {
                     b1 = set( B(k,j) );
                     xmm1 += A.load(i             ,k) * b1;
                     xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                     xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                     xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                     xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
                     xmm6 += A.load(i+SIMDSIZE*5UL,k) * b1;
                     xmm7 += A.load(i+SIMDSIZE*6UL,k) * b1;
                     xmm8 += A.load(i+SIMDSIZE*7UL,k) * b1;
                  }

                  C.store( i             , j, xmm1 );
                  C.store( i+SIMDSIZE    , j, xmm2 );
                  C.store( i+SIMDSIZE*2UL, j, xmm3 );
                  C.store( i+SIMDSIZE*3UL, j, xmm4 );
                  C.store( i+SIMDSIZE*4UL, j, xmm5 );
                  C.store( i+SIMDSIZE*5UL, j, xmm6 );
                  C.store( i+SIMDSIZE*6UL, j, xmm7 );
                  C.store( i+SIMDSIZE*7UL, j, xmm8 );
               }
               else
               {
                  const SIMDType zero;
                  C.store( i             , j, zero );
                  C.store( i+SIMDSIZE    , j, zero );
                  C.store( i+SIMDSIZE*2UL, j, zero );
                  C.store( i+SIMDSIZE*3UL, j, zero );
                  C.store( i+SIMDSIZE*4UL, j, zero );
                  C.store( i+SIMDSIZE*5UL, j, zero );
                  C.store( i+SIMDSIZE*6UL, j, zero );
                  C.store( i+SIMDSIZE*7UL, j, zero );
               }
            }
         }
      }

      for( ; !SYM && !HERM && !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1 ( a1 * b1 );
               SIMDType xmm2 ( a2 * b1 );
               SIMDType xmm3 ( a3 * b1 );
               SIMDType xmm4 ( a4 * b1 );
               SIMDType xmm5 ( a5 * b1 );
               SIMDType xmm6 ( a1 * b2 );
               SIMDType xmm7 ( a2 * b2 );
               SIMDType xmm8 ( a3 * b2 );
               SIMDType xmm9 ( a4 * b2 );
               SIMDType xmm10( a5 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  a5 = A.load(i+SIMDSIZE*4UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1  += a1 * b1;
                  xmm2  += a2 * b1;
                  xmm3  += a3 * b1;
                  xmm4  += a4 * b1;
                  xmm5  += a5 * b1;
                  xmm6  += a1 * b2;
                  xmm7  += a2 * b2;
                  xmm8  += a3 * b2;
                  xmm9  += a4 * b2;
                  xmm10 += a5 * b2;
               }

               C.store( i             , j    , xmm1  );
               C.store( i+SIMDSIZE    , j    , xmm2  );
               C.store( i+SIMDSIZE*2UL, j    , xmm3  );
               C.store( i+SIMDSIZE*3UL, j    , xmm4  );
               C.store( i+SIMDSIZE*4UL, j    , xmm5  );
               C.store( i             , j+1UL, xmm6  );
               C.store( i+SIMDSIZE    , j+1UL, xmm7  );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm8  );
               C.store( i+SIMDSIZE*3UL, j+1UL, xmm9  );
               C.store( i+SIMDSIZE*4UL, j+1UL, xmm10 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i+SIMDSIZE*3UL, j    , zero );
               C.store( i+SIMDSIZE*4UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
               C.store( i+SIMDSIZE*3UL, j+1UL, zero );
               C.store( i+SIMDSIZE*4UL, j+1UL, zero );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
               }

               C.store( i             , j, xmm1 );
               C.store( i+SIMDSIZE    , j, xmm2 );
               C.store( i+SIMDSIZE*2UL, j, xmm3 );
               C.store( i+SIMDSIZE*3UL, j, xmm4 );
               C.store( i+SIMDSIZE*4UL, j, xmm5 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
               C.store( i+SIMDSIZE*3UL, j, zero );
               C.store( i+SIMDSIZE*4UL, j, zero );
            }
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*4UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a4 * b1 );
               SIMDType xmm5( a1 * b2 );
               SIMDType xmm6( a2 * b2 );
               SIMDType xmm7( a3 * b2 );
               SIMDType xmm8( a4 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a4 * b1;
                  xmm5 += a1 * b2;
                  xmm6 += a2 * b2;
                  xmm7 += a3 * b2;
                  xmm8 += a4 * b2;
               }

               C.store( i             , j    , xmm1 );
               C.store( i+SIMDSIZE    , j    , xmm2 );
               C.store( i+SIMDSIZE*2UL, j    , xmm3 );
               C.store( i+SIMDSIZE*3UL, j    , xmm4 );
               C.store( i             , j+1UL, xmm5 );
               C.store( i+SIMDSIZE    , j+1UL, xmm6 );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm7 );
               C.store( i+SIMDSIZE*3UL, j+1UL, xmm8 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i+SIMDSIZE*3UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
               C.store( i+SIMDSIZE*3UL, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
               }

               C.store( i             , j, xmm1 );
               C.store( i+SIMDSIZE    , j, xmm2 );
               C.store( i+SIMDSIZE*2UL, j, xmm3 );
               C.store( i+SIMDSIZE*3UL, j, xmm4 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
               C.store( i+SIMDSIZE*3UL, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*3UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a1 * b2 );
               SIMDType xmm5( a2 * b2 );
               SIMDType xmm6( a3 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a1 * b2;
                  xmm5 += a2 * b2;
                  xmm6 += a3 * b2;
               }

               C.store( i             , j    , xmm1 );
               C.store( i+SIMDSIZE    , j    , xmm2 );
               C.store( i+SIMDSIZE*2UL, j    , xmm3 );
               C.store( i             , j+1UL, xmm4 );
               C.store( i+SIMDSIZE    , j+1UL, xmm5 );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm6 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               }

               C.store( i             , j, xmm1 );
               C.store( i+SIMDSIZE    , j, xmm2 );
               C.store( i+SIMDSIZE*2UL, j, xmm3 );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType b4( set( B(k,j+3UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );
               SIMDType xmm7( a1 * b4 );
               SIMDType xmm8( a2 * b4 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  b4 = set( B(k,j+3UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
                  xmm7 += a1 * b4;
                  xmm8 += a2 * b4;
               }

               C.store( i         , j    , xmm1 );
               C.store( i+SIMDSIZE, j    , xmm2 );
               C.store( i         , j+1UL, xmm3 );
               C.store( i+SIMDSIZE, j+1UL, xmm4 );
               C.store( i         , j+2UL, xmm5 );
               C.store( i+SIMDSIZE, j+2UL, xmm6 );
               C.store( i         , j+3UL, xmm7 );
               C.store( i+SIMDSIZE, j+3UL, xmm8 );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
               C.store( i         , j+2UL, zero );
               C.store( i+SIMDSIZE, j+2UL, zero );
               C.store( i         , j+3UL, zero );
               C.store( i+SIMDSIZE, j+3UL, zero );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
               }

               C.store( i         , j    , xmm1 );
               C.store( i+SIMDSIZE, j    , xmm2 );
               C.store( i         , j+1UL, xmm3 );
               C.store( i+SIMDSIZE, j+1UL, xmm4 );
               C.store( i         , j+2UL, xmm5 );
               C.store( i+SIMDSIZE, j+2UL, xmm6 );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
               C.store( i         , j+2UL, zero );
               C.store( i+SIMDSIZE, j+2UL, zero );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
               }

               C.store( i         , j    , xmm1 );
               C.store( i+SIMDSIZE, j    , xmm2 );
               C.store( i         , j+1UL, xmm3 );
               C.store( i+SIMDSIZE, j+1UL, xmm4 );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i         ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i         ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE,k) * b1;
               }

               C.store( i         , j, xmm1 );
               C.store( i+SIMDSIZE, j, xmm2 );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j, zero );
               C.store( i+SIMDSIZE, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );
               SIMDType xmm4( a1 * set( B(k,j+3UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
                  xmm4 += a1 * set( B(k,j+3UL) );
               }

               C.store( i, j    , xmm1 );
               C.store( i, j+1UL, xmm2 );
               C.store( i, j+2UL, xmm3 );
               C.store( i, j+3UL, xmm4 );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
               C.store( i, j+2UL, zero );
               C.store( i, j+3UL, zero );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
               }

               C.store( i, j    , xmm1 );
               C.store( i, j+1UL, xmm2 );
               C.store( i, j+2UL, xmm3 );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
               C.store( i, j+2UL, zero );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
               }

               C.store( i, j    , xmm1 );
               C.store( i, j+1UL, xmm2 );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               SIMDType xmm1( A.load(i,k) * set( B(k,j) ) );

               for( ++k; k<K; ++k ) {
                  xmm1 += A.load(i,k) * set( B(k,j) );
               }

               C.store( i, j, xmm1 );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; remainder && i<M; ++i )
      {
         size_t j( 0UL );

         if( SYM || HERM ) {
            for( ; j<i; ++j ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
         else if( UPP ) {
            for( ; j<i; ++j ) {
               reset( C(i,j) );
            }
         }

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               ElementType value1( A(i,k) * B(k,j    ) );
               ElementType value2( A(i,k) * B(k,j+1UL) );

               for( ++k; k<kend; ++k ) {
                  value1 += A(i,k) * B(k,j    );
                  value2 += A(i,k) * B(k,j+1UL);
               }

               C(i,j    ) = value1;
               C(i,j+1UL) = value2;
            }
            else
            {
               reset( C(i,j    ) );
               reset( C(i,j+1UL) );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               ElementType value( A(i,k) * B(k,j) );

               for( ++k; k<K; ++k ) {
                  value += A(i,k) * B(k,j);
               }

               C(i,j) = value;
            }
            else
            {
               reset( C(i,j) );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices (large matrices)***************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a large transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default assignment to dense matrices (large matrices)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default assignment of a large transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default assignment of a transpose dense matrix-
   // transpose dense matrix multiplication expression to a dense matrix. This kernel is
   // optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      if( SYM )
         smmm( C, A, B, ElementType(1) );
      else if( HERM )
         hmmm( C, A, B, ElementType(1) );
      else if( LOW )
         lmmm( C, A, B, ElementType(1), ElementType(0) );
      else if( UPP )
         ummm( C, A, B, ElementType(1), ElementType(0) );
      else
         mmm( C, A, B, ElementType(1), ElementType(0) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense matrices (default)*******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-transpose dense matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      selectLargeAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based assignment to dense matrices*****************************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function performs the transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         assign( C, B );
         trmm( C, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(1) );
      }
      else if( IsTriangular_v<MT5> ) {
         assign( C, A );
         trmm( C, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(1) );
      }
      else {
         gemm( C, A, B, ET(1), ET(0) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-transpose dense matrix multiplication to a
   //        sparse matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // transpose dense matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( serial( rhs ) );
      assign( *lhs, fwd( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment to row-major matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a row-major matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a row-major matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto assign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      assign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-transpose dense matrix multiplication
   //        to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT   // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || rhs.lhs_.columns() == 0UL ) {
         return;
      }

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTDMatMultExpr::selectAddAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices (kernel selection)************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for an addition assignment of a transpose dense matrix-
   //        transpose dense matrix multiplication to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( C, A, B );
      else
         selectBlasAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (general/general)*****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a general transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default addition assignment of a general transpose dense matrix-
   // general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t kbegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t kend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( K ) );
         BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

         for( size_t k=kbegin; k<kend; ++k )
         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j,k+1UL) : k+1UL )
                                    :( LOW ? max(j,k) : k ) )
                                 :( LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,k) : k )
                                  :( UPP ? min(j,k)+1UL : k+1UL ) )
                               :( UPP ? j+1UL : M ) );

            if( ( LOW || UPP ) && ibegin >= iend ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               C(i    ,j) += A(i    ,k) * B(k,j);
               C(i+1UL,j) += A(i+1UL,k) * B(k,j);
            }
            if( ipos < iend ) {
               C(ipos,j) += A(ipos,k) * B(k,j);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (general/diagonal)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a general transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default addition assignment of a general transpose dense matrix-
   // diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) += A(i    ,j) * B(j,j);
            C(i+1UL,j) += A(i+1UL,j) * B(j,j);
         }
         if( ipos < iend ) {
            C(ipos,j) += A(ipos,j) * B(j,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (diagonal/general)****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a diagonal transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default addition assignment of a diagonal transpose dense
   // matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) += A(i    ,i    ) * B(i    ,j);
            C(i+1UL,j) += A(i+1UL,i+1UL) * B(i+1UL,j);
         }
         if( ipos < iend ) {
            C(ipos,j) += A(ipos,ipos) * B(ipos,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (diagonal/diagonal)***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a diagonal transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default addition assignment of a diagonal transpose dense
   // matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) += A(i,i) * B(i,i);
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (small matrices)******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a small transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to row-major dense matrices (small matrices)*********
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a small transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a row-major dense matrix. This
   // kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         addAssign( C, fwd( A * tmp ) );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         addAssign( C, fwd( tmp * B ) );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         addAssign( C, fwd( A * tmp ) );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         addAssign( C, fwd( tmp * B ) );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to column-major dense matrices (small matrices)******
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a small transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a column-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               SIMDType xmm1( C.load(i             ,j) );
               SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
               SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
               SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );
               SIMDType xmm5( C.load(i+SIMDSIZE*4UL,j) );
               SIMDType xmm6( C.load(i+SIMDSIZE*5UL,j) );
               SIMDType xmm7( C.load(i+SIMDSIZE*6UL,j) );
               SIMDType xmm8( C.load(i+SIMDSIZE*7UL,j) );

               for( size_t k=kbegin; k<kend; ++k ) {
                  const SIMDType b1( set( B(k,j) ) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
                  xmm6 += A.load(i+SIMDSIZE*5UL,k) * b1;
                  xmm7 += A.load(i+SIMDSIZE*6UL,k) * b1;
                  xmm8 += A.load(i+SIMDSIZE*7UL,k) * b1;
               }

               C.store( i             , j, xmm1 );
               C.store( i+SIMDSIZE    , j, xmm2 );
               C.store( i+SIMDSIZE*2UL, j, xmm3 );
               C.store( i+SIMDSIZE*3UL, j, xmm4 );
               C.store( i+SIMDSIZE*4UL, j, xmm5 );
               C.store( i+SIMDSIZE*5UL, j, xmm6 );
               C.store( i+SIMDSIZE*6UL, j, xmm7 );
               C.store( i+SIMDSIZE*7UL, j, xmm8 );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            SIMDType xmm1 ( C.load(i             ,j    ) );
            SIMDType xmm2 ( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3 ( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4 ( C.load(i+SIMDSIZE*3UL,j    ) );
            SIMDType xmm5 ( C.load(i+SIMDSIZE*4UL,j    ) );
            SIMDType xmm6 ( C.load(i             ,j+1UL) );
            SIMDType xmm7 ( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm8 ( C.load(i+SIMDSIZE*2UL,j+1UL) );
            SIMDType xmm9 ( C.load(i+SIMDSIZE*3UL,j+1UL) );
            SIMDType xmm10( C.load(i+SIMDSIZE*4UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               const SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1  += a1 * b1;
               xmm2  += a2 * b1;
               xmm3  += a3 * b1;
               xmm4  += a4 * b1;
               xmm5  += a5 * b1;
               xmm6  += a1 * b2;
               xmm7  += a2 * b2;
               xmm8  += a3 * b2;
               xmm9  += a4 * b2;
               xmm10 += a5 * b2;
            }

            C.store( i             , j    , xmm1  );
            C.store( i+SIMDSIZE    , j    , xmm2  );
            C.store( i+SIMDSIZE*2UL, j    , xmm3  );
            C.store( i+SIMDSIZE*3UL, j    , xmm4  );
            C.store( i+SIMDSIZE*4UL, j    , xmm5  );
            C.store( i             , j+1UL, xmm6  );
            C.store( i+SIMDSIZE    , j+1UL, xmm7  );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm8  );
            C.store( i+SIMDSIZE*3UL, j+1UL, xmm9  );
            C.store( i+SIMDSIZE*4UL, j+1UL, xmm10 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );
            SIMDType xmm5( C.load(i+SIMDSIZE*4UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 += A.load(i             ,k) * b1;
               xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
               xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
            C.store( i+SIMDSIZE*3UL, j, xmm4 );
            C.store( i+SIMDSIZE*4UL, j, xmm5 );
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            SIMDType xmm1( C.load(i             ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j    ) );
            SIMDType xmm5( C.load(i             ,j+1UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm7( C.load(i+SIMDSIZE*2UL,j+1UL) );
            SIMDType xmm8( C.load(i+SIMDSIZE*3UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 += a1 * b1;
               xmm2 += a2 * b1;
               xmm3 += a3 * b1;
               xmm4 += a4 * b1;
               xmm5 += a1 * b2;
               xmm6 += a2 * b2;
               xmm7 += a3 * b2;
               xmm8 += a4 * b2;
            }

            C.store( i             , j    , xmm1 );
            C.store( i+SIMDSIZE    , j    , xmm2 );
            C.store( i+SIMDSIZE*2UL, j    , xmm3 );
            C.store( i+SIMDSIZE*3UL, j    , xmm4 );
            C.store( i             , j+1UL, xmm5 );
            C.store( i+SIMDSIZE    , j+1UL, xmm6 );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm7 );
            C.store( i+SIMDSIZE*3UL, j+1UL, xmm8 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 += A.load(i             ,k) * b1;
               xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
            C.store( i+SIMDSIZE*3UL, j, xmm4 );
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            SIMDType xmm1( C.load(i             ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4( C.load(i             ,j+1UL) );
            SIMDType xmm5( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE*2UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 += a1 * b1;
               xmm2 += a2 * b1;
               xmm3 += a3 * b1;
               xmm4 += a1 * b2;
               xmm5 += a2 * b2;
               xmm6 += a3 * b2;
            }

            C.store( i             , j    , xmm1 );
            C.store( i+SIMDSIZE    , j    , xmm2 );
            C.store( i+SIMDSIZE*2UL, j    , xmm3 );
            C.store( i             , j+1UL, xmm4 );
            C.store( i+SIMDSIZE    , j+1UL, xmm5 );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm6 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 += A.load(i             ,k) * b1;
               xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );
            SIMDType xmm5( C.load(i         ,j+2UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE,j+2UL) );
            SIMDType xmm7( C.load(i         ,j+3UL) );
            SIMDType xmm8( C.load(i+SIMDSIZE,j+3UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               const SIMDType b3( set( B(k,j+2UL) ) );
               const SIMDType b4( set( B(k,j+3UL) ) );
               xmm1 += a1 * b1;
               xmm2 += a2 * b1;
               xmm3 += a1 * b2;
               xmm4 += a2 * b2;
               xmm5 += a1 * b3;
               xmm6 += a2 * b3;
               xmm7 += a1 * b4;
               xmm8 += a2 * b4;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
            C.store( i         , j+2UL, xmm5 );
            C.store( i+SIMDSIZE, j+2UL, xmm6 );
            C.store( i         , j+3UL, xmm7 );
            C.store( i+SIMDSIZE, j+3UL, xmm8 );
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );
            SIMDType xmm5( C.load(i         ,j+2UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE,j+2UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               const SIMDType b3( set( B(k,j+2UL) ) );
               xmm1 += a1 * b1;
               xmm2 += a2 * b1;
               xmm3 += a1 * b2;
               xmm4 += a2 * b2;
               xmm5 += a1 * b3;
               xmm6 += a2 * b3;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
            C.store( i         , j+2UL, xmm5 );
            C.store( i+SIMDSIZE, j+2UL, xmm6 );
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 += a1 * b1;
               xmm2 += a2 * b1;
               xmm3 += a1 * b2;
               xmm4 += a2 * b2;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i         ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 += A.load(i         ,k) * b1;
               xmm2 += A.load(i+SIMDSIZE,k) * b1;
            }

            C.store( i         , j, xmm1 );
            C.store( i+SIMDSIZE, j, xmm2 );
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW && UPP ? min(i+SIMDSIZE,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );
            SIMDType xmm3( C.load(i,j+2UL) );
            SIMDType xmm4( C.load(i,j+3UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 += a1 * set( B(k,j    ) );
               xmm2 += a1 * set( B(k,j+1UL) );
               xmm3 += a1 * set( B(k,j+2UL) );
               xmm4 += a1 * set( B(k,j+3UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
            C.store( i, j+2UL, xmm3 );
            C.store( i, j+3UL, xmm4 );
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );
            SIMDType xmm3( C.load(i,j+2UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 += a1 * set( B(k,j    ) );
               xmm2 += a1 * set( B(k,j+1UL) );
               xmm3 += a1 * set( B(k,j+2UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
            C.store( i, j+2UL, xmm3 );
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 += a1 * set( B(k,j    ) );
               xmm2 += a1 * set( B(k,j+1UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            SIMDType xmm1( C.load(i,j) );

            for( size_t k=kbegin; k<K; ++k ) {
               xmm1 += A.load(i,k) * set( B(k,j) );
            }

            C.store( i, j, xmm1 );
         }
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jend( LOW ? i+1UL : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            ElementType value1( C(i,j    ) );
            ElementType value2( C(i,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               value1 += A(i,k) * B(k,j    );
               value2 += A(i,k) * B(k,j+1UL);
            }

            C(i,j    ) = value1;
            C(i,j+1UL) = value2;
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            ElementType value( C(i,j) );

            for( size_t k=kbegin; k<K; ++k ) {
               value += A(i,k) * B(k,j);
            }

            C(i,j) = value;
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (large matrices)******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a large transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense matrices (large matrices)*******************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default addition assignment of a large transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix. This kernel
   // is optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      if( LOW )
         lmmm( C, A, B, ElementType(1), ElementType(1) );
      else if( UPP )
         ummm( C, A, B, ElementType(1), ElementType(1) );
      else
         mmm( C, A, B, ElementType(1), ElementType(1) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense matrices (default)**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      selectLargeAddAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense matrices********************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based addition assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function performs the transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         ResultType_t<MT3> tmp( serial( B ) );
         trmm( tmp, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(1) );
         addAssign( C, tmp );
      }
      else if( IsTriangular_v<MT5> ) {
         ResultType_t<MT3> tmp( serial( A ) );
         trmm( tmp, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(1) );
         addAssign( C, tmp );
      }
      else {
         gemm( C, A, B, ET(1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Restructuring addition assignment to row-major matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a row-major matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring addition assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a row-major matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto addAssign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      addAssign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || rhs.lhs_.columns() == 0UL ) {
         return;
      }

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTDMatMultExpr::selectSubAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices (kernel selection)*********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Selection of the kernel for a subtraction assignment of a transpose dense matrix-
   //        transpose dense matrix multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( C, A, B );
      else
         selectBlasSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (general/general)**************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a general transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default subtraction assignment of a general transpose dense
   // matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t kbegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t kend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( K ) );
         BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

         for( size_t k=kbegin; k<kend; ++k )
         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j,k+1UL) : k+1UL )
                                    :( LOW ? max(j,k) : k ) )
                                 :( LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,k) : k )
                                  :( UPP ? min(j,k)+1UL : k+1UL ) )
                               :( UPP ? j+1UL : M ) );

            if( ( LOW || UPP ) && ( ibegin >= iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t inum( iend - ibegin );
            const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
            BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

            for( size_t i=ibegin; i<ipos; i+=2UL ) {
               C(i    ,j) -= A(i    ,k) * B(k,j);
               C(i+1UL,j) -= A(i+1UL,k) * B(k,j);
            }
            if( ipos < iend ) {
               C(ipos,j) -= A(ipos,k) * B(k,j);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (general/diagonal)*************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a general transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default subtraction assignment of a general transpose dense
   // matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) -= A(i    ,j) * B(j,j);
            C(i+1UL,j) -= A(i+1UL,j) * B(j,j);
         }
         if( ipos < iend ) {
            C(ipos,j) -= A(ipos,j) * B(j,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (diagonal/general)*************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a diagonal transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default subtraction assignment of a diagonal transpose dense
   // matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) -= A(i    ,i    ) * B(i    ,j);
            C(i+1UL,j) -= A(i+1UL,i+1UL) * B(i+1UL,j);
         }
         if( ipos < iend ) {
            C(ipos,j) -= A(ipos,ipos) * B(ipos,j);
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (diagonal/diagonal)************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a diagonal transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the default subtraction assignment of a diagonal transpose dense
   // matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) -= A(i,i) * B(i,i);
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (small matrices)***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a small transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to row-major dense matrices (small matrices)******
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a small transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a row-major dense matrix. This
   // kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         subAssign( C, fwd( A * tmp ) );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         subAssign( C, fwd( tmp * B ) );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         subAssign( C, fwd( A * tmp ) );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         subAssign( C, fwd( tmp * B ) );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to column-major dense matrices (small matrices)***
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a small transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a column-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               SIMDType xmm1( C.load(i             ,j) );
               SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
               SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
               SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );
               SIMDType xmm5( C.load(i+SIMDSIZE*4UL,j) );
               SIMDType xmm6( C.load(i+SIMDSIZE*5UL,j) );
               SIMDType xmm7( C.load(i+SIMDSIZE*6UL,j) );
               SIMDType xmm8( C.load(i+SIMDSIZE*7UL,j) );

               for( size_t k=kbegin; k<kend; ++k ) {
                  const SIMDType b1( set( B(k,j) ) );
                  xmm1 -= A.load(i             ,k) * b1;
                  xmm2 -= A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 -= A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 -= A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 -= A.load(i+SIMDSIZE*4UL,k) * b1;
                  xmm6 -= A.load(i+SIMDSIZE*5UL,k) * b1;
                  xmm7 -= A.load(i+SIMDSIZE*6UL,k) * b1;
                  xmm8 -= A.load(i+SIMDSIZE*7UL,k) * b1;
               }

               C.store( i             , j, xmm1 );
               C.store( i+SIMDSIZE    , j, xmm2 );
               C.store( i+SIMDSIZE*2UL, j, xmm3 );
               C.store( i+SIMDSIZE*3UL, j, xmm4 );
               C.store( i+SIMDSIZE*4UL, j, xmm5 );
               C.store( i+SIMDSIZE*5UL, j, xmm6 );
               C.store( i+SIMDSIZE*6UL, j, xmm7 );
               C.store( i+SIMDSIZE*7UL, j, xmm8 );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            SIMDType xmm1 ( C.load(i             ,j    ) );
            SIMDType xmm2 ( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3 ( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4 ( C.load(i+SIMDSIZE*3UL,j    ) );
            SIMDType xmm5 ( C.load(i+SIMDSIZE*4UL,j    ) );
            SIMDType xmm6 ( C.load(i             ,j+1UL) );
            SIMDType xmm7 ( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm8 ( C.load(i+SIMDSIZE*2UL,j+1UL) );
            SIMDType xmm9 ( C.load(i+SIMDSIZE*3UL,j+1UL) );
            SIMDType xmm10( C.load(i+SIMDSIZE*4UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               const SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1  -= a1 * b1;
               xmm2  -= a2 * b1;
               xmm3  -= a3 * b1;
               xmm4  -= a4 * b1;
               xmm5  -= a5 * b1;
               xmm6  -= a1 * b2;
               xmm7  -= a2 * b2;
               xmm8  -= a3 * b2;
               xmm9  -= a4 * b2;
               xmm10 -= a5 * b2;
            }

            C.store( i             , j    , xmm1  );
            C.store( i+SIMDSIZE    , j    , xmm2  );
            C.store( i+SIMDSIZE*2UL, j    , xmm3  );
            C.store( i+SIMDSIZE*3UL, j    , xmm4  );
            C.store( i+SIMDSIZE*4UL, j    , xmm5  );
            C.store( i             , j+1UL, xmm6  );
            C.store( i+SIMDSIZE    , j+1UL, xmm7  );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm8  );
            C.store( i+SIMDSIZE*3UL, j+1UL, xmm9  );
            C.store( i+SIMDSIZE*4UL, j+1UL, xmm10 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );
            SIMDType xmm5( C.load(i+SIMDSIZE*4UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 -= A.load(i             ,k) * b1;
               xmm2 -= A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 -= A.load(i+SIMDSIZE*2UL,k) * b1;
               xmm4 -= A.load(i+SIMDSIZE*3UL,k) * b1;
               xmm5 -= A.load(i+SIMDSIZE*4UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
            C.store( i+SIMDSIZE*3UL, j, xmm4 );
            C.store( i+SIMDSIZE*4UL, j, xmm5 );
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            SIMDType xmm1( C.load(i             ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j    ) );
            SIMDType xmm5( C.load(i             ,j+1UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm7( C.load(i+SIMDSIZE*2UL,j+1UL) );
            SIMDType xmm8( C.load(i+SIMDSIZE*3UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 -= a1 * b1;
               xmm2 -= a2 * b1;
               xmm3 -= a3 * b1;
               xmm4 -= a4 * b1;
               xmm5 -= a1 * b2;
               xmm6 -= a2 * b2;
               xmm7 -= a3 * b2;
               xmm8 -= a4 * b2;
            }

            C.store( i             , j    , xmm1 );
            C.store( i+SIMDSIZE    , j    , xmm2 );
            C.store( i+SIMDSIZE*2UL, j    , xmm3 );
            C.store( i+SIMDSIZE*3UL, j    , xmm4 );
            C.store( i             , j+1UL, xmm5 );
            C.store( i+SIMDSIZE    , j+1UL, xmm6 );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm7 );
            C.store( i+SIMDSIZE*3UL, j+1UL, xmm8 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );
            SIMDType xmm4( C.load(i+SIMDSIZE*3UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 -= A.load(i             ,k) * b1;
               xmm2 -= A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 -= A.load(i+SIMDSIZE*2UL,k) * b1;
               xmm4 -= A.load(i+SIMDSIZE*3UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
            C.store( i+SIMDSIZE*3UL, j, xmm4 );
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            SIMDType xmm1( C.load(i             ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j    ) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j    ) );
            SIMDType xmm4( C.load(i             ,j+1UL) );
            SIMDType xmm5( C.load(i+SIMDSIZE    ,j+1UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE*2UL,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i             ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               const SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 -= a1 * b1;
               xmm2 -= a2 * b1;
               xmm3 -= a3 * b1;
               xmm4 -= a1 * b2;
               xmm5 -= a2 * b2;
               xmm6 -= a3 * b2;
            }

            C.store( i             , j    , xmm1 );
            C.store( i+SIMDSIZE    , j    , xmm2 );
            C.store( i+SIMDSIZE*2UL, j    , xmm3 );
            C.store( i             , j+1UL, xmm4 );
            C.store( i+SIMDSIZE    , j+1UL, xmm5 );
            C.store( i+SIMDSIZE*2UL, j+1UL, xmm6 );
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i             ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE    ,j) );
            SIMDType xmm3( C.load(i+SIMDSIZE*2UL,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 -= A.load(i             ,k) * b1;
               xmm2 -= A.load(i+SIMDSIZE    ,k) * b1;
               xmm3 -= A.load(i+SIMDSIZE*2UL,k) * b1;
            }

            C.store( i             , j, xmm1 );
            C.store( i+SIMDSIZE    , j, xmm2 );
            C.store( i+SIMDSIZE*2UL, j, xmm3 );
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );
            SIMDType xmm5( C.load(i         ,j+2UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE,j+2UL) );
            SIMDType xmm7( C.load(i         ,j+3UL) );
            SIMDType xmm8( C.load(i+SIMDSIZE,j+3UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               const SIMDType b3( set( B(k,j+2UL) ) );
               const SIMDType b4( set( B(k,j+3UL) ) );
               xmm1 -= a1 * b1;
               xmm2 -= a2 * b1;
               xmm3 -= a1 * b2;
               xmm4 -= a2 * b2;
               xmm5 -= a1 * b3;
               xmm6 -= a2 * b3;
               xmm7 -= a1 * b4;
               xmm8 -= a2 * b4;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
            C.store( i         , j+2UL, xmm5 );
            C.store( i+SIMDSIZE, j+2UL, xmm6 );
            C.store( i         , j+3UL, xmm7 );
            C.store( i+SIMDSIZE, j+3UL, xmm8 );
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );
            SIMDType xmm5( C.load(i         ,j+2UL) );
            SIMDType xmm6( C.load(i+SIMDSIZE,j+2UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               const SIMDType b3( set( B(k,j+2UL) ) );
               xmm1 -= a1 * b1;
               xmm2 -= a2 * b1;
               xmm3 -= a1 * b2;
               xmm4 -= a2 * b2;
               xmm5 -= a1 * b3;
               xmm6 -= a2 * b3;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
            C.store( i         , j+2UL, xmm5 );
            C.store( i+SIMDSIZE, j+2UL, xmm6 );
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            SIMDType xmm1( C.load(i         ,j    ) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j    ) );
            SIMDType xmm3( C.load(i         ,j+1UL) );
            SIMDType xmm4( C.load(i+SIMDSIZE,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i         ,k) );
               const SIMDType a2( A.load(i+SIMDSIZE,k) );
               const SIMDType b1( set( B(k,j    ) ) );
               const SIMDType b2( set( B(k,j+1UL) ) );
               xmm1 -= a1 * b1;
               xmm2 -= a2 * b1;
               xmm3 -= a1 * b2;
               xmm4 -= a2 * b2;
            }

            C.store( i         , j    , xmm1 );
            C.store( i+SIMDSIZE, j    , xmm2 );
            C.store( i         , j+1UL, xmm3 );
            C.store( i+SIMDSIZE, j+1UL, xmm4 );
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            SIMDType xmm1( C.load(i         ,j) );
            SIMDType xmm2( C.load(i+SIMDSIZE,j) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType b1( set( B(k,j) ) );
               xmm1 -= A.load(i         ,k) * b1;
               xmm2 -= A.load(i+SIMDSIZE,k) * b1;
            }

            C.store( i         , j, xmm1 );
            C.store( i+SIMDSIZE, j, xmm2 );
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW && UPP ? min(i+SIMDSIZE,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );
            SIMDType xmm3( C.load(i,j+2UL) );
            SIMDType xmm4( C.load(i,j+3UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 -= a1 * set( B(k,j    ) );
               xmm2 -= a1 * set( B(k,j+1UL) );
               xmm3 -= a1 * set( B(k,j+2UL) );
               xmm4 -= a1 * set( B(k,j+3UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
            C.store( i, j+2UL, xmm3 );
            C.store( i, j+3UL, xmm4 );
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );
            SIMDType xmm3( C.load(i,j+2UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 -= a1 * set( B(k,j    ) );
               xmm2 -= a1 * set( B(k,j+1UL) );
               xmm3 -= a1 * set( B(k,j+2UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
            C.store( i, j+2UL, xmm3 );
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            SIMDType xmm1( C.load(i,j    ) );
            SIMDType xmm2( C.load(i,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               const SIMDType a1( A.load(i,k) );
               xmm1 -= a1 * set( B(k,j    ) );
               xmm2 -= a1 * set( B(k,j+1UL) );
            }

            C.store( i, j    , xmm1 );
            C.store( i, j+1UL, xmm2 );
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            SIMDType xmm1( C.load(i,j) );

            for( size_t k=kbegin; k<K; ++k ) {
               xmm1 -= A.load(i,k) * set( B(k,j) );
            }

            C.store( i, j, xmm1 );
         }
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jend( LOW ? i+1UL : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            ElementType value1( C(i,j    ) );
            ElementType value2( C(i,j+1UL) );

            for( size_t k=kbegin; k<kend; ++k ) {
               value1 -= A(i,k) * B(k,j    );
               value2 -= A(i,k) * B(k,j+1UL);
            }

            C(i,j    ) = value1;
            C(i,j+1UL) = value2;
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            ElementType value( C(i,j) );

            for( size_t k=kbegin; k<K; ++k ) {
               value -= A(i,k) * B(k,j);
            }

            C(i,j) = value;
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (large matrices)***************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a large transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      selectDefaultSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense matrices (large matrices)****************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized default subtraction assignment of a large transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix. This kernel is
   // optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5> >
   {
      if( LOW )
         lmmm( C, A, B, ElementType(-1), ElementType(1) );
      else if( UPP )
         ummm( C, A, B, ElementType(-1), ElementType(1) );
      else
         mmm( C, A, B, ElementType(-1), ElementType(1) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense matrices (default)*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      selectLargeSubAssignKernel( C, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**BLAS-based subraction assignment to dense matrices******************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*! \cond BLAZE_INTERNAL */
   /*!\brief BLAS-based subraction assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \return void
   //
   // This function performs the transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectBlasSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         ResultType_t<MT3> tmp( serial( B ) );
         trmm( tmp, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(1) );
         subAssign( C, tmp );
      }
      else if( IsTriangular_v<MT5> ) {
         ResultType_t<MT3> tmp( serial( A ) );
         trmm( tmp, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(1) );
         subAssign( C, tmp );
      }
      else {
         gemm( C, A, B, ET(-1), ET(1) );
      }
   }
   /*! \endcond */
#endif
   //**********************************************************************************************

   //**Restructuring subtraction assignment to row-major matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a transpose dense matrix-transpose dense
   //        matrix multiplication to a row-major matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring subtraction assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only
   // be selected by the compiler in case the symmetry of either of the two matrix operands
   // can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto subAssign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      subAssign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Multiplication assignment to dense matrices*************************************************
   // No special implementation for the multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**Multiplication assignment to sparse matrices************************************************
   // No special implementation for the multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP assignment to dense matrices************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-transpose dense matrix multiplication to a
   //        dense matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL ) {
         return;
      }
      else if( rhs.lhs_.columns() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpAssign( *lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose dense matrix-transpose dense matrix multiplication
   //        to a sparse matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-transpose dense matrix multiplication expression to a sparse matrix. Due to the
   // explicit application of the SFINAE principle, this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( rhs );
      smpAssign( *lhs, fwd( tmp ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP assignment to row-major matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a row-major matrix (\f$ C=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a row-major matrix. Due
   // to the explicit application of the SFINAE principle this function can only be selected by
   // the compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAssign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      smpAssign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT   // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || rhs.lhs_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpAddAssign( *lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP addition assignment to row-major matrices*********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP addition assignment of a transpose dense matrix-transpose dense
   //        matrix multiplication to a row-major matrix (\f$ C+=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP addition assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only be
   // selected by the compiler in case the symmetry of either of the two matrix operands can be
   // exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAddAssign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      smpAddAssign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense
   // matrix. Due to the explicit application of the SFINAE principle, this function can
   // only be selected by the compiler in case either of the two matrix operands requires
   // an intermediate evaluation and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || rhs.lhs_.columns() == 0UL ) {
         return;
      }

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      smpSubAssign( *lhs, A * B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring SMP subtraction assignment to row-major matrices******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP subtraction assignment of a transpose dense matrix-transpose dense
   //        matrix multiplication to a row-major matrix (\f$ C-=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP subtraction assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only be
   // selected by the compiler in case the symmetry of either of the two matrix operands can be
   // exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpSubAssign( Matrix<MT,false>& lhs, const TDMatTDMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.lhs_ ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.rhs_ ) );

      smpSubAssign( *lhs, fwd( A * B ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense
   // matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const TDMatTDMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
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

   //**SMP multiplication assignment to dense matrices*********************************************
   // No special implementation for the SMP multiplication assignment to dense matrices.
   //**********************************************************************************************

   //**SMP multiplication assignment to sparse matrices********************************************
   // No special implementation for the SMP multiplication assignment to sparse matrices.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATMULTEXPR( MT1, MT2 );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DMATSCALARMULTEXPR SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Expression object for scaled transpose dense matrix-transpose dense matrix multiplications.
// \ingroup dense_matrix_expression
//
// This specialization of the DMatScalarMultExpr class represents the compile time expression
// for scaled multiplications between two column-major dense matrices.
*/
template< typename MT1   // Type of the left-hand side dense matrix
        , typename MT2   // Type of the right-hand side dense matrix
        , bool SF        // Symmetry flag
        , bool HF        // Hermitian flag
        , bool LF        // Lower flag
        , bool UF        // Upper flag
        , typename ST >  // Type of the right-hand side scalar value
class DMatScalarMultExpr< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>, ST, true >
   : public MatScalarMultExpr< DenseMatrix< DMatScalarMultExpr< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>, ST, true >, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   //! Type of the dense matrix multiplication expression.
   using MMM = TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>;

   using RES = ResultType_t<MMM>;     //!< Result type of the dense matrix multiplication expression.
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side dense matrix expression.
   using ET1 = ElementType_t<RT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<RT2>;    //!< Element type of the right-hand side dense matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateLeft = ( IsComputation_v<MT1> || RequiresEvaluation_v<MT1> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side dense matrix expression.
   static constexpr bool evaluateRight = ( IsComputation_v<MT2> || RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool SYM  = ( SF && !( HF || LF || UF )    );  //!< Flag for symmetric matrices.
   static constexpr bool HERM = ( HF && !( LF || UF )          );  //!< Flag for Hermitian matrices.
   static constexpr bool LOW  = ( LF || ( ( SF || HF ) && UF ) );  //!< Flag for lower matrices.
   static constexpr bool UPP  = ( UF || ( ( SF || HF ) && LF ) );  //!< Flag for upper matrices.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the optimal evaluation strategy.
       In case the target matrix is row-major and either of the two matrix operands is symmetric,
       the variable is set to 1 and an optimized evaluation strategy is selected. Otherwise the
       variable is set to 0 and the default strategy is chosen. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool CanExploitSymmetry_v =
      ( IsRowMajorMatrix_v<T1> && ( IsSymmetric_v<T2> || IsSymmetric_v<T3> ) );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either of the two matrix operands requires an intermediate evaluation, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool IsEvaluationRequired_v =
      ( ( evaluateLeft || evaluateRight ) && !CanExploitSymmetry_v<T1,T2,T3> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case the types of all three involved matrices and the scalar type are suited for a BLAS
       kernel, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   static constexpr bool UseBlasKernel_v =
      ( BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION &&
        !SYM && !HERM && !LOW && !UPP &&
        IsContiguous_v<T1> && HasMutableDataAccess_v<T1> &&
        IsContiguous_v<T2> && HasConstDataAccess_v<T2> &&
        IsContiguous_v<T3> && HasConstDataAccess_v<T3> &&
        !IsDiagonal_v<T2> && !IsDiagonal_v<T3> &&
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
   /*! In case all four involved data types are suited for a vectorized computation of the
       matrix multiplication, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3, typename T4 >
   static constexpr bool UseVectorizedDefaultKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled && T3::simdEnabled &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3>
                          , T4 > &&
        HasSIMDAdd_v< ElementType_t<T2>, ElementType_t<T2> > &&
        HasSIMDMult_v< ElementType_t<T3>, ElementType_t<T3> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Type of the functor for forwarding an expression to another assign kernel.
   /*! In case a temporary matrix needs to be created, this functor is used to forward the
       resulting expression to another assign kernel. */
   using ForwardFunctor = If_t< HERM
                              , DeclHerm
                              , If_t< SYM
                                    , DeclSym
                                    , If_t< LOW
                                          , If_t< UPP
                                                , DeclDiag
                                                , DeclLow >
                                          , If_t< UPP
                                                , DeclUpp
                                                , Noop > > > >;
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatScalarMultExpr instance.
   using This = DMatScalarMultExpr<MMM,ST,true>;

   //! Base type of this DMatScalarMultExpr instance.
   using BaseType = MatScalarMultExpr< DenseMatrix<This,true> >;

   //! Result type for expression template evaluations.
   using ResultType = typename If_t< HERM
                                   , DeclHermTrait< MultTrait_t<RES,ST> >
                                   , If_t< SYM
                                         , DeclSymTrait< MultTrait_t<RES,ST> >
                                         , If_t< LOW
                                               , If_t< UPP
                                                     , DeclDiagTrait< MultTrait_t<RES,ST> >
                                                     , DeclLowTrait< MultTrait_t<RES,ST> > >
                                               , If_t< UPP
                                                     , DeclUppTrait< MultTrait_t<RES,ST> >
                                                     , MultTrait<RES,ST> > > > >::Type;

   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using SIMDType      = SIMDTrait_t<ElementType>;     //!< Resulting SIMD element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side dense matrix expression.
   using LeftOperand = const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>;

   //! Composite type of the right-hand side scalar value.
   using RightOperand = ST;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateLeft, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side dense matrix operand.
   using RT = If_t< evaluateRight, const RT2, CT2 >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( !IsDiagonal_v<MT1> &&
        MT1::simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v<ET1,ET2,ST> &&
        HasSIMDAdd_v<ET1,ET2> &&
        HasSIMDMult_v<ET1,ET2> );

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateLeft  && MT1::smpAssignable && !evaluateRight && MT2::smpAssignable );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatScalarMultExpr class.
   //
   // \param matrix The left-hand side dense matrix of the multiplication expression.
   // \param scalar The right-hand side scalar of the multiplication expression.
   */
   inline DMatScalarMultExpr( const MMM& matrix, ST scalar )
      : matrix_( matrix )  // Left-hand side dense matrix of the multiplication expression
      , scalar_( scalar )  // Right-hand side scalar of the multiplication expression
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
      return matrix_(i,j) * scalar_;
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

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const {
      return matrix_.rows();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const {
      return matrix_.columns();
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side dense matrix operand.
   //
   // \return The left-hand side dense matrix operand.
   */
   inline LeftOperand leftOperand() const {
      return matrix_;
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
   inline bool isAliased( const T* alias ) const {
      return matrix_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const {
      return matrix_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( !BLAZE_BLAS_MODE ||
               !BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION ||
               !BLAZE_BLAS_IS_PARALLEL ||
               ( rows() * columns() < TDMATTDMATMULT_THRESHOLD ) ) &&
             ( rows() * columns() >= SMP_TDMATTDMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  matrix_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand scalar_;  //!< Right-hand side scalar of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-transpose dense matrix multiplication
   //        to a dense matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( right ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      DMatScalarMultExpr::selectAssignKernel( *lhs, A, B, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Assignment to dense matrices (kernel selection)*********************************************
   /*!\brief Selection of the kernel for an assignment of a scaled transpose dense matrix-
   //        transpose dense matrix multiplication to a dense matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallAssignKernel( C, A, B, scalar );
      else
         selectBlasAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (general/general)**************************************
   /*!\brief Default assignment of a scaled general transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default assignment of a scaled general transpose dense matrix-
   // general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( SYM || HERM || LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t kbegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t kend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( K ) );
         BLAZE_INTERNAL_ASSERT( kbegin <= kend, "Invalid loop indices detected" );

         if( IsStrictlyTriangular_v<MT5> && kbegin == kend ) {
            for( size_t i=0UL; i<M; ++i ) {
               reset( C(i,j) );
            }
            continue;
         }

         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( LOW ? max(j,kbegin+1UL) : kbegin+1UL )
                                    :( LOW ? max(j,kbegin) : kbegin ) )
                                 :( LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,kbegin) : kbegin )
                                  :( UPP ? min(j,kbegin)+1UL : kbegin+1UL ) )
                               :( UPP ? j+1UL : M ) );

            if( ( IsLower_v<MT4> && IsLower_v<MT5> ) || LOW ) {
               for( size_t i=0UL; i<ibegin; ++i ) {
                  reset( C(i,j) );
               }
            }
            else if( IsStrictlyLower_v<MT4> ) {
               reset( C(0UL,j) );
            }
            for( size_t i=ibegin; i<iend; ++i ) {
               C(i,j) = A(i,kbegin) * B(kbegin,j);
            }
            if( ( IsUpper_v<MT4> && IsUpper_v<MT5> ) || UPP ) {
               for( size_t i=iend; i<M; ++i ) {
                  reset( C(i,j) );
               }
            }
            else if( IsStrictlyUpper_v<MT4> ) {
               reset( C(M-1UL,j) );
            }
         }

         for( size_t k=kbegin+1UL; k<kend; ++k )
         {
            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( SYM || HERM || LOW ? max( j, k+1UL ) : k+1UL )
                                    :( SYM || HERM || LOW ? max( j, k ) : k ) )
                                 :( SYM || HERM || LOW ? j : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? min(j+1UL,k-1UL) : k-1UL )
                                  :( UPP ? min(j+1UL,k) : k ) )
                               :( UPP ? j+1UL : M ) );

            if( ( SYM || HERM || LOW || UPP ) && ( ibegin > iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            for( size_t i=ibegin; i<iend; ++i ) {
               C(i,j) += A(i,k) * B(k,j);
            }
            if( IsUpper_v<MT4> ) {
               C(iend,j) = A(iend,k) * B(k,j);
            }
         }

         {
            const size_t ibegin( ( IsLower_v<MT4> && IsLower_v<MT5> )
                                 ?( IsStrictlyLower_v<MT4> || IsStrictlyLower_v<MT5> ? j+1UL : j )
                                 :( ( SYM || HERM || LOW )?( j ):( 0UL ) ) );
            const size_t iend( ( IsUpper_v<MT4> && IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT4> || IsStrictlyUpper_v<MT5> ? j : j+1UL )
                               :( UPP ? j+1UL : M ) );

            if( ( SYM || HERM || LOW || UPP ) && ( ibegin > iend ) ) continue;
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            for( size_t i=ibegin; i<iend; ++i ) {
               C(i,j) *= scalar;
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=1UL; j<N; ++j ) {
            for( size_t i=0UL; i<j; ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (general/diagonal)*************************************
   /*!\brief Default assignment of a scaled general transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default assignment of a scaled general transpose dense matrix-
   // diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         if( IsLower_v<MT4> ) {
            for( size_t i=0UL; i<ibegin; ++i ) {
               reset( C(i,j) );
            }
         }
         for( size_t i=ibegin; i<iend; ++i ) {
            C(i,j) = A(i,j) * B(j,j) * scalar;
         }
         if( IsUpper_v<MT4> ) {
            for( size_t i=iend; i<M; ++i ) {
               reset( C(i,j) );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (diagonal/general)*************************************
   /*!\brief Default assignment of a scaled diagonal transpose dense matrix-general transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default assignment of a scaled diagonal transpose dense matrix-
   // general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         if( IsLower_v<MT4> ) {
            for( size_t i=0UL; i<ibegin; ++i ) {
               reset( C(i,j) );
            }
         }
         for( size_t i=ibegin; i<iend; ++i ) {
            C(i,j) = A(i,i) * B(i,j) * scalar;
         }
         if( IsUpper_v<MT4> ) {
            for( size_t i=iend; i<M; ++i ) {
               reset( C(i,j) );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (diagonal/diagonal)************************************
   /*!\brief Default assignment of a scaled diagonal transpose dense matrix-diagonal transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default assignment of a scaled diagonal transpose dense matrix-
   // diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      reset( C );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) = A(i,i) * B(i,i) * scalar;
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (small matrices)***************************************
   /*!\brief Default assignment of a small scaled transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to row-major dense matrices (small matrices)******************
   /*!\brief Vectorized default assignment of a small scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default assignment of a scaled transpose dense
   // matrix-transpose dense matrix multiplication expression to a row-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         assign( C, fwd( A * tmp ) * scalar );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         assign( C, fwd( tmp * B ) * scalar );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         assign( C, fwd( A * tmp ) * scalar );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         assign( C, fwd( tmp * B ) * scalar );
      }
   }
   //**********************************************************************************************

   //**Vectorized default assignment to column-major dense matrices (small matrices)***************
   /*!\brief Vectorized default assignment of a small scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default assignment of a scaled transpose dense
   // matrix-transpose dense matrix multiplication expression to a column-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( SYM || HERM || LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !SYM && !HERM && !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               size_t k( kbegin );

               if( k < kend )
               {
                  SIMDType b1( set( B(k,j) ) );
                  SIMDType xmm1( A.load(i             ,k) * b1 );
                  SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
                  SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
                  SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
                  SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );
                  SIMDType xmm6( A.load(i+SIMDSIZE*5UL,k) * b1 );
                  SIMDType xmm7( A.load(i+SIMDSIZE*6UL,k) * b1 );
                  SIMDType xmm8( A.load(i+SIMDSIZE*7UL,k) * b1 );

                  for( ++k; k<kend; ++k ) {
                     b1 = set( B(k,j) );
                     xmm1 += A.load(i             ,k) * b1;
                     xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                     xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                     xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                     xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
                     xmm6 += A.load(i+SIMDSIZE*5UL,k) * b1;
                     xmm7 += A.load(i+SIMDSIZE*6UL,k) * b1;
                     xmm8 += A.load(i+SIMDSIZE*7UL,k) * b1;
                  }

                  C.store( i             , j, xmm1 * factor );
                  C.store( i+SIMDSIZE    , j, xmm2 * factor );
                  C.store( i+SIMDSIZE*2UL, j, xmm3 * factor );
                  C.store( i+SIMDSIZE*3UL, j, xmm4 * factor );
                  C.store( i+SIMDSIZE*4UL, j, xmm5 * factor );
                  C.store( i+SIMDSIZE*5UL, j, xmm6 * factor );
                  C.store( i+SIMDSIZE*6UL, j, xmm7 * factor );
                  C.store( i+SIMDSIZE*7UL, j, xmm8 * factor );
               }
               else
               {
                  const SIMDType zero;
                  C.store( i             , j, zero );
                  C.store( i+SIMDSIZE    , j, zero );
                  C.store( i+SIMDSIZE*2UL, j, zero );
                  C.store( i+SIMDSIZE*3UL, j, zero );
                  C.store( i+SIMDSIZE*4UL, j, zero );
                  C.store( i+SIMDSIZE*5UL, j, zero );
                  C.store( i+SIMDSIZE*6UL, j, zero );
                  C.store( i+SIMDSIZE*7UL, j, zero );
               }
            }
         }
      }

      for( ; !SYM && !HERM && !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1 ( a1 * b1 );
               SIMDType xmm2 ( a2 * b1 );
               SIMDType xmm3 ( a3 * b1 );
               SIMDType xmm4 ( a4 * b1 );
               SIMDType xmm5 ( a5 * b1 );
               SIMDType xmm6 ( a1 * b2 );
               SIMDType xmm7 ( a2 * b2 );
               SIMDType xmm8 ( a3 * b2 );
               SIMDType xmm9 ( a4 * b2 );
               SIMDType xmm10( a5 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  a5 = A.load(i+SIMDSIZE*4UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1  += a1 * b1;
                  xmm2  += a2 * b1;
                  xmm3  += a3 * b1;
                  xmm4  += a4 * b1;
                  xmm5  += a5 * b1;
                  xmm6  += a1 * b2;
                  xmm7  += a2 * b2;
                  xmm8  += a3 * b2;
                  xmm9  += a4 * b2;
                  xmm10 += a5 * b2;
               }

               C.store( i             , j    , xmm1  * factor );
               C.store( i+SIMDSIZE    , j    , xmm2  * factor );
               C.store( i+SIMDSIZE*2UL, j    , xmm3  * factor );
               C.store( i+SIMDSIZE*3UL, j    , xmm4  * factor );
               C.store( i+SIMDSIZE*4UL, j    , xmm5  * factor );
               C.store( i             , j+1UL, xmm6  * factor );
               C.store( i+SIMDSIZE    , j+1UL, xmm7  * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm8  * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, xmm9  * factor );
               C.store( i+SIMDSIZE*4UL, j+1UL, xmm10 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i+SIMDSIZE*3UL, j    , zero );
               C.store( i+SIMDSIZE*4UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
               C.store( i+SIMDSIZE*3UL, j+1UL, zero );
               C.store( i+SIMDSIZE*4UL, j+1UL, zero );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
               }

               C.store( i             , j, xmm1 * factor );
               C.store( i+SIMDSIZE    , j, xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, xmm4 * factor );
               C.store( i+SIMDSIZE*4UL, j, xmm5 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
               C.store( i+SIMDSIZE*3UL, j, zero );
               C.store( i+SIMDSIZE*4UL, j, zero );
            }
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*4UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a4 * b1 );
               SIMDType xmm5( a1 * b2 );
               SIMDType xmm6( a2 * b2 );
               SIMDType xmm7( a3 * b2 );
               SIMDType xmm8( a4 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a4 * b1;
                  xmm5 += a1 * b2;
                  xmm6 += a2 * b2;
                  xmm7 += a3 * b2;
                  xmm8 += a4 * b2;
               }

               C.store( i             , j    , xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j    , xmm4 * factor );
               C.store( i             , j+1UL, xmm5 * factor );
               C.store( i+SIMDSIZE    , j+1UL, xmm6 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm7 * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, xmm8 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i+SIMDSIZE*3UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
               C.store( i+SIMDSIZE*3UL, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
               }

               C.store( i             , j, xmm1 * factor );
               C.store( i+SIMDSIZE    , j, xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, xmm4 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
               C.store( i+SIMDSIZE*3UL, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*4UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*3UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a1 * b2 );
               SIMDType xmm5( a2 * b2 );
               SIMDType xmm6( a3 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a1 * b2;
                  xmm5 += a2 * b2;
                  xmm6 += a3 * b2;
               }

               C.store( i             , j    , xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , xmm3 * factor );
               C.store( i             , j+1UL, xmm4 * factor );
               C.store( i+SIMDSIZE    , j+1UL, xmm5 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, xmm6 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j    , zero );
               C.store( i+SIMDSIZE    , j    , zero );
               C.store( i+SIMDSIZE*2UL, j    , zero );
               C.store( i             , j+1UL, zero );
               C.store( i+SIMDSIZE    , j+1UL, zero );
               C.store( i+SIMDSIZE*2UL, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               }

               C.store( i             , j, xmm1 * factor );
               C.store( i+SIMDSIZE    , j, xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, xmm3 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i             , j, zero );
               C.store( i+SIMDSIZE    , j, zero );
               C.store( i+SIMDSIZE*2UL, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*3UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType b4( set( B(k,j+3UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );
               SIMDType xmm7( a1 * b4 );
               SIMDType xmm8( a2 * b4 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  b4 = set( B(k,j+3UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
                  xmm7 += a1 * b4;
                  xmm8 += a2 * b4;
               }

               C.store( i         , j    , xmm1 * factor );
               C.store( i+SIMDSIZE, j    , xmm2 * factor );
               C.store( i         , j+1UL, xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, xmm4 * factor );
               C.store( i         , j+2UL, xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, xmm6 * factor );
               C.store( i         , j+3UL, xmm7 * factor );
               C.store( i+SIMDSIZE, j+3UL, xmm8 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
               C.store( i         , j+2UL, zero );
               C.store( i+SIMDSIZE, j+2UL, zero );
               C.store( i         , j+3UL, zero );
               C.store( i+SIMDSIZE, j+3UL, zero );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
               }

               C.store( i         , j    , xmm1 * factor );
               C.store( i+SIMDSIZE, j    , xmm2 * factor );
               C.store( i         , j+1UL, xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, xmm4 * factor );
               C.store( i         , j+2UL, xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, xmm6 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
               C.store( i         , j+2UL, zero );
               C.store( i+SIMDSIZE, j+2UL, zero );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
               }

               C.store( i         , j    , xmm1 * factor );
               C.store( i+SIMDSIZE, j    , xmm2 * factor );
               C.store( i         , j+1UL, xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, xmm4 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j    , zero );
               C.store( i+SIMDSIZE, j    , zero );
               C.store( i         , j+1UL, zero );
               C.store( i+SIMDSIZE, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i         ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i         ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE,k) * b1;
               }

               C.store( i         , j, xmm1 * factor );
               C.store( i+SIMDSIZE, j, xmm2 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i         , j, zero );
               C.store( i+SIMDSIZE, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE*2UL,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE,N) : N );
         size_t j( 0UL );

         if( SYM || HERM ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  C(ii,j) = HERM ? conj( C(j,ii) ) : C(j,ii);
               }
            }
         }
         else if( UPP ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<i; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );
               SIMDType xmm4( a1 * set( B(k,j+3UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
                  xmm4 += a1 * set( B(k,j+3UL) );
               }

               C.store( i, j    , xmm1 * factor );
               C.store( i, j+1UL, xmm2 * factor );
               C.store( i, j+2UL, xmm3 * factor );
               C.store( i, j+3UL, xmm4 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
               C.store( i, j+2UL, zero );
               C.store( i, j+3UL, zero );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
               }

               C.store( i, j    , xmm1 * factor );
               C.store( i, j+1UL, xmm2 * factor );
               C.store( i, j+2UL, xmm3 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
               C.store( i, j+2UL, zero );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
               }

               C.store( i, j    , xmm1 * factor );
               C.store( i, j+1UL, xmm2 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j    , zero );
               C.store( i, j+1UL, zero );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               SIMDType xmm1( A.load(i,k) * set( B(k,j) ) );

               for( ++k; k<K; ++k ) {
                  xmm1 += A.load(i,k) * set( B(k,j) );
               }

               C.store( i, j, xmm1 * factor );
            }
            else
            {
               const SIMDType zero;
               C.store( i, j, zero );
            }

            if( LOW ) ++j;
         }

         if( LOW ) {
            const size_t iiend( min(i+SIMDSIZE,M) );
            for( ; j<N; ++j ) {
               for( size_t ii=i; ii<iiend; ++ii ) {
                  reset( C(ii,j) );
               }
            }
         }
      }

      for( ; remainder && i<M; ++i )
      {
         size_t j( 0UL );

         if( SYM || HERM ) {
            for( ; j<i; ++j ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
         else if( UPP ) {
            for( ; j<i; ++j ) {
               reset( C(i,j) );
            }
         }

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               ElementType value1( A(i,k) * B(k,j    ) );
               ElementType value2( A(i,k) * B(k,j+1UL) );

               for( ++k; k<kend; ++k ) {
                  value1 += A(i,k) * B(k,j    );
                  value2 += A(i,k) * B(k,j+1UL);
               }

               C(i,j    ) = value1 * scalar;
               C(i,j+1UL) = value2 * scalar;
            }
            else
            {
               reset( C(i,j    ) );
               reset( C(i,j+1UL) );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               ElementType value( A(i,k) * B(k,j) );

               for( ++k; k<K; ++k ) {
                  value += A(i,k) * B(k,j);
               }

               C(i,j) = value * scalar;
            }
            else
            {
               reset( C(i,j) );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default assignment to dense matrices (large matrices)***************************************
   /*!\brief Default assignment of a large scaled transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default assignment to dense matrices (large matrices)****************************
   /*!\brief Vectorized default assignment of a large scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default assignment of a scaled transpose dense
   // matrix-transpose dense matrix multiplication expression to a dense matrix. This kernel
   // is optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      if( SYM )
         smmm( C, A, B, scalar );
      else if( HERM )
         hmmm( C, A, B, scalar );
      else if( LOW )
         lmmm( C, A, B, scalar, ST2(0) );
      else if( UPP )
         ummm( C, A, B, scalar, ST2(0) );
      else
         mmm( C, A, B, scalar, ST2(0) );
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense matrices (default)*******************************************
   /*!\brief Default assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the assignment of a large scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectLargeAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based assignment to dense matrices*****************************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*!\brief BLAS-based assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         assign( C, B );
         trmm( C, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(scalar) );
      }
      else if( IsTriangular_v<MT5> ) {
         assign( C, A );
         trmm( C, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(scalar) );
      }
      else {
         gemm( C, A, B, ET(scalar), ET(0) );
      }
   }
#endif
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*!\brief Assignment of a scaled transpose dense matrix-transpose dense matrix multiplication
   //        to a sparse matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a scaled transpose dense
   // matrix-transpose dense matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( serial( rhs ) );
      assign( *lhs, fwd( tmp ) );
   }
   //**********************************************************************************************

   //**Restructuring assignment to row-major matrices**********************************************
   /*!\brief Restructuring assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a row-major matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a row-major matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto assign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      assign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*!\brief Addition assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto addAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( right ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      DMatScalarMultExpr::selectAddAssignKernel( *lhs, A, B, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Addition assignment to dense matrices (kernel selection)************************************
   /*!\brief Selection of the kernel for an addition assignment of a scaled transpose dense
   //        matrix-transpose dense matrix multiplication to a dense matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallAddAssignKernel( C, A, B, scalar );
      else
         selectBlasAddAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (general/general)*****************************
   /*!\brief Default addition assignment of a scaled general transpose dense matrix-general
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment of a scaled general transpose
   // dense matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const ResultType tmp( serial( A * B * scalar ) );
      addAssign( C, tmp );
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (general/diagonal)****************************
   /*!\brief Default addition assignment of a scaled general transpose dense matrix-diagonal
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment of a scaled general transpose
   // dense matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) += A(i    ,j) * B(j,j) * scalar;
            C(i+1UL,j) += A(i+1UL,j) * B(j,j) * scalar;
         }
         if( ipos < iend ) {
            C(ipos,j) += A(ipos,j) * B(j,j) * scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (diagonal/general)****************************
   /*!\brief Default addition assignment of a scaled diagonal transpose dense matrix-general
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment of a scaled diagonal transpose
   // dense matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) += A(i    ,i    ) * B(i    ,j) * scalar;
            C(i+1UL,j) += A(i+1UL,i+1UL) * B(i+1UL,j) * scalar;
         }
         if( ipos < iend ) {
            C(ipos,j) += A(ipos,ipos) * B(ipos,j) * scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (diagonal/diagonal)***************************
   /*!\brief Default addition assignment of a scaled diagonal transpose dense matrix-diagonal
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default addition assignment of a scaled diagonal transpose
   // dense matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) += A(i,i) * B(i,i) * scalar;
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (small matrices)******************************
   /*!\brief Default addition assignment of a small scaled transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultAddAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to row-major dense matrices (small matrices)*********
   /*!\brief Vectorized default addition assignment of a small scaled transpose dense matrix-
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a row-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         addAssign( C, fwd( A * tmp ) * scalar );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         addAssign( C, fwd( tmp * B ) * scalar );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         addAssign( C, fwd( A * tmp ) * scalar );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         addAssign( C, fwd( tmp * B ) * scalar );
      }
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to column-major dense matrices (small matrices)******
   /*!\brief Vectorized default addition assignment of a small scaled transpose dense matrix-
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a column-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               size_t k( kbegin );

               if( k < kend )
               {
                  SIMDType b1( set( B(k,j) ) );
                  SIMDType xmm1( A.load(i             ,k) * b1 );
                  SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
                  SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
                  SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
                  SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );
                  SIMDType xmm6( A.load(i+SIMDSIZE*5UL,k) * b1 );
                  SIMDType xmm7( A.load(i+SIMDSIZE*6UL,k) * b1 );
                  SIMDType xmm8( A.load(i+SIMDSIZE*7UL,k) * b1 );

                  for( ++k; k<kend; ++k ) {
                     b1 = set( B(k,j) );
                     xmm1 += A.load(i             ,k) * b1;
                     xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                     xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                     xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                     xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
                     xmm6 += A.load(i+SIMDSIZE*5UL,k) * b1;
                     xmm7 += A.load(i+SIMDSIZE*6UL,k) * b1;
                     xmm8 += A.load(i+SIMDSIZE*7UL,k) * b1;
                  }

                  C.store( i             , j, C.load(i             ,j) + xmm1 * factor );
                  C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) + xmm2 * factor );
                  C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) + xmm3 * factor );
                  C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) + xmm4 * factor );
                  C.store( i+SIMDSIZE*4UL, j, C.load(i+SIMDSIZE*4UL,j) + xmm5 * factor );
                  C.store( i+SIMDSIZE*5UL, j, C.load(i+SIMDSIZE*5UL,j) + xmm6 * factor );
                  C.store( i+SIMDSIZE*6UL, j, C.load(i+SIMDSIZE*6UL,j) + xmm7 * factor );
                  C.store( i+SIMDSIZE*7UL, j, C.load(i+SIMDSIZE*7UL,j) + xmm8 * factor );
               }
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1 ( a1 * b1 );
               SIMDType xmm2 ( a2 * b1 );
               SIMDType xmm3 ( a3 * b1 );
               SIMDType xmm4 ( a4 * b1 );
               SIMDType xmm5 ( a5 * b1 );
               SIMDType xmm6 ( a1 * b2 );
               SIMDType xmm7 ( a2 * b2 );
               SIMDType xmm8 ( a3 * b2 );
               SIMDType xmm9 ( a4 * b2 );
               SIMDType xmm10( a5 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  a5 = A.load(i+SIMDSIZE*4UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1  += a1 * b1;
                  xmm2  += a2 * b1;
                  xmm3  += a3 * b1;
                  xmm4  += a4 * b1;
                  xmm5  += a5 * b1;
                  xmm6  += a1 * b2;
                  xmm7  += a2 * b2;
                  xmm8  += a3 * b2;
                  xmm9  += a4 * b2;
                  xmm10 += a5 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) + xmm1  * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) + xmm2  * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) + xmm3  * factor );
               C.store( i+SIMDSIZE*3UL, j    , C.load(i+SIMDSIZE*3UL,j    ) + xmm4  * factor );
               C.store( i+SIMDSIZE*4UL, j    , C.load(i+SIMDSIZE*4UL,j    ) + xmm5  * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) + xmm6  * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) + xmm7  * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) + xmm8  * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, C.load(i+SIMDSIZE*3UL,j+1UL) + xmm9  * factor );
               C.store( i+SIMDSIZE*4UL, j+1UL, C.load(i+SIMDSIZE*4UL,j+1UL) + xmm10 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) + xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) + xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) + xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) + xmm4 * factor );
               C.store( i+SIMDSIZE*4UL, j, C.load(i+SIMDSIZE*4UL,j) + xmm5 * factor );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a4 * b1 );
               SIMDType xmm5( a1 * b2 );
               SIMDType xmm6( a2 * b2 );
               SIMDType xmm7( a3 * b2 );
               SIMDType xmm8( a4 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a4 * b1;
                  xmm5 += a1 * b2;
                  xmm6 += a2 * b2;
                  xmm7 += a3 * b2;
                  xmm8 += a4 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) + xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) + xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) + xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j    , C.load(i+SIMDSIZE*3UL,j    ) + xmm4 * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) + xmm5 * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) + xmm6 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) + xmm7 * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, C.load(i+SIMDSIZE*3UL,j+1UL) + xmm8 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) + xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) + xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) + xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) + xmm4 * factor );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a1 * b2 );
               SIMDType xmm5( a2 * b2 );
               SIMDType xmm6( a3 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a1 * b2;
                  xmm5 += a2 * b2;
                  xmm6 += a3 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) + xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) + xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) + xmm3 * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) + xmm4 * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) + xmm5 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) + xmm6 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) + xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) + xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) + xmm3 * factor );
            }
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType b4( set( B(k,j+3UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );
               SIMDType xmm7( a1 * b4 );
               SIMDType xmm8( a2 * b4 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  b4 = set( B(k,j+3UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
                  xmm7 += a1 * b4;
                  xmm8 += a2 * b4;
               }

               C.store( i         , j    , C.load(i         ,j    ) + xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) + xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) + xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) + xmm4 * factor );
               C.store( i         , j+2UL, C.load(i         ,j+2UL) + xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, C.load(i+SIMDSIZE,j+2UL) + xmm6 * factor );
               C.store( i         , j+3UL, C.load(i         ,j+3UL) + xmm7 * factor );
               C.store( i+SIMDSIZE, j+3UL, C.load(i+SIMDSIZE,j+3UL) + xmm8 * factor );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
               }

               C.store( i         , j    , C.load(i         ,j    ) + xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) + xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) + xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) + xmm4 * factor );
               C.store( i         , j+2UL, C.load(i         ,j+2UL) + xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, C.load(i+SIMDSIZE,j+2UL) + xmm6 * factor );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
               }

               C.store( i         , j    , C.load(i         ,j    ) + xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) + xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) + xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) + xmm4 * factor );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i         ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i         ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE,k) * b1;
               }

               C.store( i         , j, C.load(i         ,j) + xmm1 * factor );
               C.store( i+SIMDSIZE, j, C.load(i+SIMDSIZE,j) + xmm2 * factor );
            }
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW && UPP ? min(i+SIMDSIZE,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );
               SIMDType xmm4( a1 * set( B(k,j+3UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
                  xmm4 += a1 * set( B(k,j+3UL) );
               }

               C.store( i, j    , C.load(i,j    ) + xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) + xmm2 * factor );
               C.store( i, j+2UL, C.load(i,j+2UL) + xmm3 * factor );
               C.store( i, j+3UL, C.load(i,j+3UL) + xmm4 * factor );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
               }

               C.store( i, j    , C.load(i,j    ) + xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) + xmm2 * factor );
               C.store( i, j+2UL, C.load(i,j+2UL) + xmm3 * factor );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
               }

               C.store( i, j    , C.load(i,j    ) + xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) + xmm2 * factor );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               SIMDType xmm1( A.load(i,k) * set( B(k,j) ) );

               for( ++k; k<K; ++k ) {
                  xmm1 += A.load(i,k) * set( B(k,j) );
               }

               C.store( i, j, C.load(i,j) + xmm1 * factor );
            }
         }
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jend( LOW ? i+1UL : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               ElementType value1( A(i,k) * B(k,j    ) );
               ElementType value2( A(i,k) * B(k,j+1UL) );

               for( ++k; k<kend; ++k ) {
                  value1 += A(i,k) * B(k,j    );
                  value2 += A(i,k) * B(k,j+1UL);
               }

               C(i,j    ) += value1 * scalar;
               C(i,j+1UL) += value2 * scalar;
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               ElementType value( A(i,k) * B(k,j) );

               for( ++k; k<K; ++k ) {
                  value += A(i,k) * B(k,j);
               }

               C(i,j) += value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default addition assignment to dense matrices (large matrices)******************************
   /*!\brief Default addition assignment of a large scaled transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultAddAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default addition assignment to dense matrices (large matrices)*******************
   /*!\brief Vectorized default addition assignment of a large scaled transpose dense matrix-
   //        transpose dense matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default addition assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix. This
   // kernel is optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      if( LOW )
         lmmm( C, A, B, scalar, ST2(1) );
      else if( UPP )
         ummm( C, A, B, scalar, ST2(1) );
      else
         mmm( C, A, B, scalar, ST2(1) );
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense matrices (default)**********************************
   /*!\brief Default addition assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the addition assignment of a large
   // scaled transpose dense matrix-transpose dense matrix multiplication expression to a dense
   // matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectLargeAddAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based addition assignment to dense matrices********************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*!\brief BLAS-based addition assignment of a scaled transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasAddAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         ResultType_t<MT3> tmp( serial( B ) );
         trmm( tmp, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(scalar) );
         addAssign( C, tmp );
      }
      else if( IsTriangular_v<MT5> ) {
         ResultType_t<MT3> tmp( serial( A ) );
         trmm( tmp, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(scalar) );
         addAssign( C, tmp );
      }
      else {
         gemm( C, A, B, ET(scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Restructuring addition assignment to row-major matrices*************************************
   /*!\brief Restructuring addition assignment of a scaled transpose dense matrix-transpose dense
   //        matrix multiplication to a row-major matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring addition assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only be
   // selected by the compiler in case the symmetry of either of the two matrix operands can be
   // exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto addAssign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      addAssign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*!\brief Subtraction assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto subAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( serial( left  ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( right ) );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      DMatScalarMultExpr::selectSubAssignKernel( *lhs, A, B, rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices (kernel selection)*********************************
   /*!\brief Selection of the kernel for a subtraction assignment of a scaled transpose dense
   //        matrix-transpose dense matrix multiplication to a dense matrix (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline void selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
   {
      if( ( IsDiagonal_v<MT4> ) ||
          ( !BLAZE_DEBUG_MODE && A.rows() <= SIMDSIZE*10UL ) ||
          ( C.rows() * C.columns() < TDMATTDMATMULT_THRESHOLD ) )
         selectSmallSubAssignKernel( C, A, B, scalar );
      else
         selectBlasSubAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (general/general)**************************
   /*!\brief Default subtraction assignment of a scaled general transpose dense matrix-general
   //        transpose dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment of a scaled general transpose
   // dense matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      const ResultType tmp( serial( A * B * scalar ) );
      subAssign( C, tmp );
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (general/diagonal)*************************
   /*!\brief Default subtraction assignment of a scaled general transpose dense matrix-diagonal
   //        transpose dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment of a scaled general transpose
   // dense matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< !IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT4> )
                              ?( IsStrictlyLower_v<MT4> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT4> )
                            ?( IsStrictlyUpper_v<MT4> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) -= A(i    ,j) * B(j,j) * scalar;
            C(i+1UL,j) -= A(i+1UL,j) * B(j,j) * scalar;
         }
         if( ipos < iend ) {
            C(ipos,j) -= A(ipos,j) * B(j,j) * scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (diagonal/general)*************************
   /*!\brief Default subtraction assignment of a scaled diagonal transpose dense matrix-general
   //        transpose dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment of a scaled diagonal transpose
   // dense matrix-general transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && !IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );

      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( ( IsLower_v<MT5> )
                              ?( IsStrictlyLower_v<MT5> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend( ( IsUpper_v<MT5> )
                            ?( IsStrictlyUpper_v<MT5> ? j : j+1UL )
                            :( M ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         const size_t inum( iend - ibegin );
         const size_t ipos( ibegin + prevMultiple( inum, 2UL ) );
         BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

         for( size_t i=ibegin; i<ipos; i+=2UL ) {
            C(i    ,j) -= A(i    ,i    ) * B(i    ,j) * scalar;
            C(i+1UL,j) -= A(i+1UL,i+1UL) * B(i+1UL,j) * scalar;
         }
         if( ipos < iend ) {
            C(ipos,j) -= A(ipos,ipos) * B(ipos,j) * scalar;
         }
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (diagonal/diagonal)************************
   /*!\brief Default subtraction assignment of a scaled diagonal transpose dense matrix-diagonal
   //        transpose dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the default subtraction assignment of a scaled diagonal transpose
   // dense matrix-diagonal transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectDefaultSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsDiagonal_v<MT4> && IsDiagonal_v<MT5> >
   {
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT3 );

      for( size_t i=0UL; i<A.rows(); ++i ) {
         C(i,i) -= A(i,i) * B(i,i) * scalar;
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (small matrices)***************************
   /*!\brief Default subtraction assignment of a small scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultSubAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to row-major dense matrices (small matrices)******
   /*!\brief Default subtraction assignment of a small scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a row-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsRowMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT4 );
      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT5 );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT4> );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType_t<MT5> );

      const ForwardFunctor fwd;

      if( IsResizable_v<MT4> && !IsResizable_v<MT5> ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         subAssign( C, fwd( A * tmp ) * scalar );
      }
      else if( !IsResizable_v<MT4> && IsResizable_v<MT5> ) {
         const OppositeType_t<MT4> tmp( serial( A ) );
         subAssign( C, fwd( tmp * B ) * scalar );
      }
      else if( B.rows() * B.columns() <= A.rows() * A.columns() ) {
         const OppositeType_t<MT5> tmp( serial( B ) );
         subAssign( C, fwd( A * tmp ) * scalar );
      }
      else {
         const OppositeType_t<MT4> tmp( serial( A ) );
         subAssign( C, fwd( tmp * B ) * scalar );
      }
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to column-major dense matrices (small matrices)***
   /*!\brief Default subtraction assignment of a small scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a column-major dense matrix.
   // This kernel is optimized for small matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectSmallSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< IsColumnMajorMatrix_v<MT3> && UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      const size_t M( A.rows()    );
      const size_t N( B.columns() );
      const size_t K( A.columns() );

      BLAZE_INTERNAL_ASSERT( !( LOW || UPP ) || ( M == N ), "Broken invariant detected" );

      const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
      BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

      const SIMDType factor( set( scalar ) );

      size_t i( 0UL );

      if( IsIntegral_v<ElementType> )
      {
         for( ; !LOW && !UPP && (i+SIMDSIZE*7UL) < ipos; i+=SIMDSIZE*8UL ) {
            for( size_t j=0UL; j<N; ++j )
            {
               const size_t kbegin( ( IsLower_v<MT5> )
                                    ?( ( IsUpper_v<MT4> )
                                       ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                       :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                    :( IsUpper_v<MT4> ? i : 0UL ) );
               const size_t kend( ( IsUpper_v<MT5> )
                                  ?( ( IsLower_v<MT4> )
                                     ?( min( i+SIMDSIZE*8UL, K, ( IsStrictlyUpper_v<MT5> ? j : j+1UL ) ) )
                                     :( IsStrictlyUpper_v<MT5> ? j : j+1UL ) )
                                  :( IsLower_v<MT4> ? min( i+SIMDSIZE*8UL, K ) : K ) );

               size_t k( kbegin );

               if( k < kend )
               {
                  SIMDType b1( set( B(k,j) ) );
                  SIMDType xmm1( A.load(i             ,k) * b1 );
                  SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
                  SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
                  SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
                  SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );
                  SIMDType xmm6( A.load(i+SIMDSIZE*5UL,k) * b1 );
                  SIMDType xmm7( A.load(i+SIMDSIZE*6UL,k) * b1 );
                  SIMDType xmm8( A.load(i+SIMDSIZE*7UL,k) * b1 );

                  for( ++k; k<kend; ++k ) {
                     b1 = set( B(k,j) );
                     xmm1 += A.load(i             ,k) * b1;
                     xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                     xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                     xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                     xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
                     xmm6 += A.load(i+SIMDSIZE*5UL,k) * b1;
                     xmm7 += A.load(i+SIMDSIZE*6UL,k) * b1;
                     xmm8 += A.load(i+SIMDSIZE*7UL,k) * b1;
                  }

                  C.store( i             , j, C.load(i             ,j) - xmm1 * factor );
                  C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) - xmm2 * factor );
                  C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) - xmm3 * factor );
                  C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) - xmm4 * factor );
                  C.store( i+SIMDSIZE*4UL, j, C.load(i+SIMDSIZE*4UL,j) - xmm5 * factor );
                  C.store( i+SIMDSIZE*5UL, j, C.load(i+SIMDSIZE*5UL,j) - xmm6 * factor );
                  C.store( i+SIMDSIZE*6UL, j, C.load(i+SIMDSIZE*6UL,j) - xmm7 * factor );
                  C.store( i+SIMDSIZE*7UL, j, C.load(i+SIMDSIZE*7UL,j) - xmm8 * factor );
               }
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*4UL) < ipos; i+=SIMDSIZE*5UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*5UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*5UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType a5( A.load(i+SIMDSIZE*4UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1 ( a1 * b1 );
               SIMDType xmm2 ( a2 * b1 );
               SIMDType xmm3 ( a3 * b1 );
               SIMDType xmm4 ( a4 * b1 );
               SIMDType xmm5 ( a5 * b1 );
               SIMDType xmm6 ( a1 * b2 );
               SIMDType xmm7 ( a2 * b2 );
               SIMDType xmm8 ( a3 * b2 );
               SIMDType xmm9 ( a4 * b2 );
               SIMDType xmm10( a5 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  a5 = A.load(i+SIMDSIZE*4UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1  += a1 * b1;
                  xmm2  += a2 * b1;
                  xmm3  += a3 * b1;
                  xmm4  += a4 * b1;
                  xmm5  += a5 * b1;
                  xmm6  += a1 * b2;
                  xmm7  += a2 * b2;
                  xmm8  += a3 * b2;
                  xmm9  += a4 * b2;
                  xmm10 += a5 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) - xmm1  * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) - xmm2  * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) - xmm3  * factor );
               C.store( i+SIMDSIZE*3UL, j    , C.load(i+SIMDSIZE*3UL,j    ) - xmm4  * factor );
               C.store( i+SIMDSIZE*4UL, j    , C.load(i+SIMDSIZE*4UL,j    ) - xmm5  * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) - xmm6  * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) - xmm7  * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) - xmm8  * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, C.load(i+SIMDSIZE*3UL,j+1UL) - xmm9  * factor );
               C.store( i+SIMDSIZE*4UL, j+1UL, C.load(i+SIMDSIZE*4UL,j+1UL) - xmm10 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*5UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );
               SIMDType xmm5( A.load(i+SIMDSIZE*4UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
                  xmm5 += A.load(i+SIMDSIZE*4UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) - xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) - xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) - xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) - xmm4 * factor );
               C.store( i+SIMDSIZE*4UL, j, C.load(i+SIMDSIZE*4UL,j) - xmm5 * factor );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*4UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*4UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType a4( A.load(i+SIMDSIZE*3UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a4 * b1 );
               SIMDType xmm5( a1 * b2 );
               SIMDType xmm6( a2 * b2 );
               SIMDType xmm7( a3 * b2 );
               SIMDType xmm8( a4 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  a4 = A.load(i+SIMDSIZE*3UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a4 * b1;
                  xmm5 += a1 * b2;
                  xmm6 += a2 * b2;
                  xmm7 += a3 * b2;
                  xmm8 += a4 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) - xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) - xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) - xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j    , C.load(i+SIMDSIZE*3UL,j    ) - xmm4 * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) - xmm5 * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) - xmm6 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) - xmm7 * factor );
               C.store( i+SIMDSIZE*3UL, j+1UL, C.load(i+SIMDSIZE*3UL,j+1UL) - xmm8 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*4UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );
               SIMDType xmm4( A.load(i+SIMDSIZE*3UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
                  xmm4 += A.load(i+SIMDSIZE*3UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) - xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) - xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) - xmm3 * factor );
               C.store( i+SIMDSIZE*3UL, j, C.load(i+SIMDSIZE*3UL,j) - xmm4 * factor );
            }
         }
      }

      for( ; !LOW && !UPP && (i+SIMDSIZE*2UL) < ipos; i+=SIMDSIZE*3UL )
      {
         size_t j( 0UL );

         for( ; (j+2UL) <= N; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*3UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*3UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i             ,k) );
               SIMDType a2( A.load(i+SIMDSIZE    ,k) );
               SIMDType a3( A.load(i+SIMDSIZE*2UL,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a3 * b1 );
               SIMDType xmm4( a1 * b2 );
               SIMDType xmm5( a2 * b2 );
               SIMDType xmm6( a3 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i             ,k);
                  a2 = A.load(i+SIMDSIZE    ,k);
                  a3 = A.load(i+SIMDSIZE*2UL,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a3 * b1;
                  xmm4 += a1 * b2;
                  xmm5 += a2 * b2;
                  xmm6 += a3 * b2;
               }

               C.store( i             , j    , C.load(i             ,j    ) - xmm1 * factor );
               C.store( i+SIMDSIZE    , j    , C.load(i+SIMDSIZE    ,j    ) - xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j    , C.load(i+SIMDSIZE*2UL,j    ) - xmm3 * factor );
               C.store( i             , j+1UL, C.load(i             ,j+1UL) - xmm4 * factor );
               C.store( i+SIMDSIZE    , j+1UL, C.load(i+SIMDSIZE    ,j+1UL) - xmm5 * factor );
               C.store( i+SIMDSIZE*2UL, j+1UL, C.load(i+SIMDSIZE*2UL,j+1UL) - xmm6 * factor );
            }
         }

         if( j < N )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*3UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i             ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE    ,k) * b1 );
               SIMDType xmm3( A.load(i+SIMDSIZE*2UL,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i             ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE    ,k) * b1;
                  xmm3 += A.load(i+SIMDSIZE*2UL,k) * b1;
               }

               C.store( i             , j, C.load(i             ,j) - xmm1 * factor );
               C.store( i+SIMDSIZE    , j, C.load(i+SIMDSIZE    ,j) - xmm2 * factor );
               C.store( i+SIMDSIZE*2UL, j, C.load(i+SIMDSIZE*2UL,j) - xmm3 * factor );
            }
         }
      }

      for( ; !( LOW && UPP ) && (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL )
      {
         const size_t jend( LOW ? min(i+SIMDSIZE*2UL,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType b4( set( B(k,j+3UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );
               SIMDType xmm7( a1 * b4 );
               SIMDType xmm8( a2 * b4 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  b4 = set( B(k,j+3UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
                  xmm7 += a1 * b4;
                  xmm8 += a2 * b4;
               }

               C.store( i         , j    , C.load(i         ,j    ) - xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) - xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) - xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) - xmm4 * factor );
               C.store( i         , j+2UL, C.load(i         ,j+2UL) - xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, C.load(i+SIMDSIZE,j+2UL) - xmm6 * factor );
               C.store( i         , j+3UL, C.load(i         ,j+3UL) - xmm7 * factor );
               C.store( i+SIMDSIZE, j+3UL, C.load(i+SIMDSIZE,j+3UL) - xmm8 * factor );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType b3( set( B(k,j+2UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );
               SIMDType xmm5( a1 * b3 );
               SIMDType xmm6( a2 * b3 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  b3 = set( B(k,j+2UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
                  xmm5 += a1 * b3;
                  xmm6 += a2 * b3;
               }

               C.store( i         , j    , C.load(i         ,j    ) - xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) - xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) - xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) - xmm4 * factor );
               C.store( i         , j+2UL, C.load(i         ,j+2UL) - xmm5 * factor );
               C.store( i+SIMDSIZE, j+2UL, C.load(i+SIMDSIZE,j+2UL) - xmm6 * factor );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( ( IsLower_v<MT4> )
                                  ?( min( i+SIMDSIZE*2UL, K, ( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) ) )
                                  :( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL ) )
                               :( IsLower_v<MT4> ? min( i+SIMDSIZE*2UL, K ) : K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i         ,k) );
               SIMDType a2( A.load(i+SIMDSIZE,k) );
               SIMDType b1( set( B(k,j    ) ) );
               SIMDType b2( set( B(k,j+1UL) ) );
               SIMDType xmm1( a1 * b1 );
               SIMDType xmm2( a2 * b1 );
               SIMDType xmm3( a1 * b2 );
               SIMDType xmm4( a2 * b2 );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i         ,k);
                  a2 = A.load(i+SIMDSIZE,k);
                  b1 = set( B(k,j    ) );
                  b2 = set( B(k,j+1UL) );
                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
                  xmm3 += a1 * b2;
                  xmm4 += a2 * b2;
               }

               C.store( i         , j    , C.load(i         ,j    ) - xmm1 * factor );
               C.store( i+SIMDSIZE, j    , C.load(i+SIMDSIZE,j    ) - xmm2 * factor );
               C.store( i         , j+1UL, C.load(i         ,j+1UL) - xmm3 * factor );
               C.store( i+SIMDSIZE, j+1UL, C.load(i+SIMDSIZE,j+1UL) - xmm4 * factor );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsLower_v<MT4> )?( min( i+SIMDSIZE*2UL, K ) ):( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType b1( set( B(k,j) ) );
               SIMDType xmm1( A.load(i         ,k) * b1 );
               SIMDType xmm2( A.load(i+SIMDSIZE,k) * b1 );

               for( ++k; k<kend; ++k ) {
                  b1 = set( B(k,j) );
                  xmm1 += A.load(i         ,k) * b1;
                  xmm2 += A.load(i+SIMDSIZE,k) * b1;
               }

               C.store( i         , j, C.load(i         ,j) - xmm1 * factor );
               C.store( i+SIMDSIZE, j, C.load(i+SIMDSIZE,j) - xmm2 * factor );
            }
         }
      }

      for( ; i<ipos; i+=SIMDSIZE )
      {
         const size_t jend( LOW && UPP ? min(i+SIMDSIZE,N) : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+4UL) <= jend; j+=4UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+3UL : j+4UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );
               SIMDType xmm4( a1 * set( B(k,j+3UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
                  xmm4 += a1 * set( B(k,j+3UL) );
               }

               C.store( i, j    , C.load(i,j    ) - xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) - xmm2 * factor );
               C.store( i, j+2UL, C.load(i,j+2UL) - xmm3 * factor );
               C.store( i, j+3UL, C.load(i,j+3UL) - xmm4 * factor );
            }
         }

         for( ; (j+3UL) <= jend; j+=3UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+2UL : j+3UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );
               SIMDType xmm3( a1 * set( B(k,j+2UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
                  xmm3 += a1 * set( B(k,j+2UL) );
               }

               C.store( i, j    , C.load(i,j    ) - xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) - xmm2 * factor );
               C.store( i, j+2UL, C.load(i,j+2UL) - xmm3 * factor );
            }
         }

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               SIMDType a1( A.load(i,k) );
               SIMDType xmm1( a1 * set( B(k,j    ) ) );
               SIMDType xmm2( a1 * set( B(k,j+1UL) ) );

               for( ++k; k<kend; ++k ) {
                  a1 = A.load(i,k);
                  xmm1 += a1 * set( B(k,j    ) );
                  xmm2 += a1 * set( B(k,j+1UL) );
               }

               C.store( i, j    , C.load(i,j    ) - xmm1 * factor );
               C.store( i, j+1UL, C.load(i,j+1UL) - xmm2 * factor );
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               SIMDType xmm1( A.load(i,k) * set( B(k,j) ) );

               for( ++k; k<K; ++k ) {
                  xmm1 += A.load(i,k) * set( B(k,j) );
               }

               C.store( i, j, C.load(i,j) - xmm1 * factor );
            }
         }
      }

      for( ; remainder && i<M; ++i )
      {
         const size_t jend( LOW ? i+1UL : N );
         size_t j( UPP ? i : 0UL );

         for( ; (j+2UL) <= jend; j+=2UL )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );
            const size_t kend( ( IsUpper_v<MT5> )
                               ?( IsStrictlyUpper_v<MT5> ? j+1UL : j+2UL )
                               :( K ) );

            size_t k( kbegin );

            if( k < kend )
            {
               ElementType value1( A(i,k) * B(k,j    ) );
               ElementType value2( A(i,k) * B(k,j+1UL) );

               for( ++k; k<kend; ++k ) {
                  value1 += A(i,k) * B(k,j    );
                  value2 += A(i,k) * B(k,j+1UL);
               }

               C(i,j    ) -= value1 * scalar;
               C(i,j+1UL) -= value2 * scalar;
            }
         }

         if( j < jend )
         {
            const size_t kbegin( ( IsLower_v<MT5> )
                                 ?( ( IsUpper_v<MT4> )
                                    ?( max( i, ( IsStrictlyLower_v<MT5> ? j+1UL : j ) ) )
                                    :( IsStrictlyLower_v<MT5> ? j+1UL : j ) )
                                 :( IsUpper_v<MT4> ? i : 0UL ) );

            size_t k( kbegin );

            if( k < K )
            {
               ElementType value( A(i,k) * B(k,j) );

               for( ++k; k<K; ++k ) {
                  value += A(i,k) * B(k,j);
               }

               C(i,j) -= value * scalar;
            }
         }
      }
   }
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices (large matrices)***************************
   /*!\brief Default subtraction assignment of a large scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectDefaultSubAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**Vectorized default subtraction assignment to dense matrices (large matrices)****************
   /*!\brief Default subtraction assignment of a large scaled transpose dense matrix-transpose
   //        dense matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function implements the vectorized default subtraction assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix. This kernel
   // is optimized for large matrices.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectLargeSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseVectorizedDefaultKernel_v<MT3,MT4,MT5,ST2> >
   {
      if( LOW )
         lmmm( C, A, B, -scalar, ST2(1) );
      else if( UPP )
         ummm( C, A, B, -scalar, ST2(1) );
      else
         mmm( C, A, B, -scalar, ST2(1) );
   }
   //**********************************************************************************************

   //**BLAS-based subtraction assignment to dense matrices (default)*******************************
   /*!\brief Default subtraction assignment of a scaled transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function relays to the default implementation of the subtraction assignment of a large
   // scaled transpose dense matrix-transpose dense matrix multiplication expression to a dense
   // matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> DisableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      selectLargeSubAssignKernel( C, A, B, scalar );
   }
   //**********************************************************************************************

   //**BLAS-based subraction assignment to dense matrices******************************************
#if BLAZE_BLAS_MODE && BLAZE_USE_BLAS_MATRIX_MATRIX_MULTIPLICATION
   /*!\brief BLAS-based subraction assignment of a scaled transpose dense matrix-transpose dense
   //        matrix multiplication (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side multiplication operand.
   // \param B The right-hand side multiplication operand.
   // \param scalar The scaling factor.
   // \return void
   //
   // This function performs the scaled transpose dense matrix-transpose dense matrix multiplication
   // based on the according BLAS functionality.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5    // Type of the right-hand side matrix operand
           , typename ST2 >  // Type of the scalar value
   static inline auto selectBlasSubAssignKernel( MT3& C, const MT4& A, const MT5& B, ST2 scalar )
      -> EnableIf_t< UseBlasKernel_v<MT3,MT4,MT5,ST2> >
   {
      using ET = ElementType_t<MT3>;

      if( IsTriangular_v<MT4> ) {
         ResultType_t<MT3> tmp( serial( B ) );
         trmm( tmp, A, CblasLeft, ( IsLower_v<MT4> )?( CblasLower ):( CblasUpper ), ET(scalar) );
         subAssign( C, tmp );
      }
      else if( IsTriangular_v<MT5> ) {
         ResultType_t<MT3> tmp( serial( A ) );
         trmm( tmp, B, CblasRight, ( IsLower_v<MT5> )?( CblasLower ):( CblasUpper ), ET(scalar) );
         subAssign( C, tmp );
      }
      else {
         gemm( C, A, B, ET(-scalar), ET(1) );
      }
   }
#endif
   //**********************************************************************************************

   //**Restructuring subtraction assignment to row-major matrices**********************************
   /*!\brief Restructuring subtraction assignment of a scaled transpose dense matrix-transpose
   //        dense matrix multiplication to a row-major matrix (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring subtraction assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a row-major matrix.
   // Due to the explicit application of the SFINAE principle this function can only be selected
   // by the compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto subAssign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      subAssign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*!\brief Schur product assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C\circ=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      schurAssign( *lhs, tmp );
   }
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
   /*!\brief SMP assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL ) {
         return;
      }
      else if( left.columns() == 0UL ) {
         reset( *lhs );
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT B( right );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      smpAssign( *lhs, A * B * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP assignment to sparse matrices***********************************************************
   /*!\brief SMP assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a sparse matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a scaled transpose
   // dense matrix-transpose dense matrix multiplication expression to a sparse matrix. Due to
   // the explicit application of the SFINAE principle, this function can only be selected by
   // the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      using TmpType = If_t< SO, ResultType, OppositeType >;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( OppositeType );
      BLAZE_CONSTRAINT_MATRICES_MUST_HAVE_SAME_STORAGE_ORDER( MT, TmpType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( TmpType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      const TmpType tmp( rhs );
      smpAssign( *lhs, fwd( tmp ) );
   }
   //**********************************************************************************************

   //**Restructuring SMP assignment to row-major matrices******************************************
   /*!\brief Restructuring SMP assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a row-major matrix (\f$ C=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP assignment of a scaled
   // transpose dense matrix-dense matrix multiplication expression to a row-major matrix. Due
   // to the explicit application of the SFINAE principle this function can only be selected by
   // the compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAssign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      smpAssign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*!\brief SMP addition assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   // Due to the explicit application of the SFINAE principle, this function can only be selected
   // by the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT B( right );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      smpAddAssign( *lhs, A * B * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Restructuring SMP addition assignment to row-major matrices*********************************
   /*!\brief Restructuring SMP addition assignment of a scaled transpose dense matrix-transpose
   //        dense matrix multiplication to a row-major matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP addition assignment of a
   // scaled transpose dense matrix-transpose dense matrix multiplication expression to a
   // row-major matrix. Due to the explicit application of the SFINAE principle this operator
   // can only be selected by the compiler in case the symmetry of either of the two matrix
   // operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAddAssign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      smpAddAssign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*!\brief SMP subtraction assignment of a scaled transpose dense matrix-transpose dense matrix
   //        multiplication to a dense matrix (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   // Due to the explicit application of the SFINAE principle, this function can only be selected
   // by the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LeftOperand_t<MMM>  left ( rhs.matrix_.leftOperand()  );
      RightOperand_t<MMM> right( rhs.matrix_.rightOperand() );

      if( (*lhs).rows() == 0UL || (*lhs).columns() == 0UL || left.columns() == 0UL ) {
         return;
      }

      LT A( left  );  // Evaluation of the left-hand side dense matrix operand
      RT B( right );  // Evaluation of the right-hand side dense matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == left.rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == left.columns()  , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == right.rows()    , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == right.columns() , "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns(), "Invalid number of columns" );

      smpSubAssign( *lhs, A * B * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**Restructuring SMP subtraction assignment to row-major matrices******************************
   /*!\brief Restructuring SMP subtraction assignment of a scaled transpose dense matrix-transpose
   //        dense matrix multiplication to a row-major matrix (\f$ C-=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side scaled multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP subtraction assignment of
   // a scaled transpose dense matrix-transpose dense matrix multiplication expression to a
   // row-major matrix. Due to the explicit application of the SFINAE principle this operator
   // can only be selected by the compiler in case the symmetry of either of the two matrix
   // operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpSubAssign( Matrix<MT,false>& lhs, const DMatScalarMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ForwardFunctor fwd;

      decltype(auto) A( transIf< IsSymmetric_v<MT1> >( rhs.matrix_.leftOperand()  ) );
      decltype(auto) B( transIf< IsSymmetric_v<MT2> >( rhs.matrix_.rightOperand() ) );

      smpSubAssign( *lhs, fwd( A * B ) * rhs.scalar_ );
   }
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*!\brief SMP Schur product assignment of a scaled transpose dense matrix-transpose dense
   //        matrix multiplication to a dense matrix (\f$ C+=s*A*B \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a scaled
   // transpose dense matrix-transpose dense matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const DMatScalarMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( rhs );
      smpSchurAssign( *lhs, tmp );
   }
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
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MMM );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MMM );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
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
/*!\brief Multiplication operator for the multiplication of two column-major dense matrices
//        (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the multiplication.
// \param rhs The right-hand side matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of two column-major dense matrices:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current number of columns of \a lhs and the current number of rows of \a rhs
// don't match, a \a std::invalid_argument is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
inline decltype(auto)
   operator*( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).columns() != (*rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,false,false,false,false>;
   return ReturnType( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-symmetric matrix multiplication expression as symmetric.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid symmetric matrix specification.
//
// The \a declsym function declares the given non-symmetric matrix multiplication expression
// \a dm as symmetric. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declsym function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declsym( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declsym( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid symmetric matrix specification" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,true,HF,LF,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-Hermitian matrix multiplication expression as Hermitian.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid Hermitian matrix specification.
//
// The \a declherm function declares the given non-Hermitian matrix multiplication expression
// \a dm as Hermitian. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declherm function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declherm( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declherm( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid Hermitian matrix specification" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,SF,true,LF,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-lower matrix multiplication expression as lower.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid lower matrix specification.
//
// The \a decllow function declares the given non-lower matrix multiplication expression
// \a dm as lower. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a decllow function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = decllow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) decllow( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,SF,HF,true,UF>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-unilower matrix multiplication expression as unilower.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid unilower matrix specification.
//
// The \a declunilow function declares the given non-unilower matrix multiplication expression
// \a dm as unilower. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a declunilow function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declunilow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool UF >     // Upper flag
inline decltype(auto) declunilow( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,false,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   return declunilow( decllow( *dm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-strictly-lower matrix multiplication expression as strictly lower.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid strictly lower matrix specification.
//
// The \a declstrlow function declares the given non-strictly-lower matrix multiplication
// expression \a dm as strictly lower. The function returns an expression representing the
// operation. In case the given expression does not represent a square matrix, a
// \a std::invalid_argument exception is thrown.\n
// The following example demonstrates the use of the \a declstrlow function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declstrlow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool UF >     // Upper flag
inline decltype(auto) declstrlow( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,false,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   return declstrlow( decllow( *dm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-upper matrix multiplication expression as upper.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid upper matrix specification.
//
// The \a declupp function declares the given non-upper matrix multiplication expression
// \a dm as upper. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a declupp function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) declupp( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,true>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-uniupper matrix multiplication expression as uniupper.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid uniupper matrix specification.
//
// The \a decluniupp function declares the given non-uniupper matrix multiplication expression
// \a dm as uniupper. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument
// exception is thrown.\n
// The following example demonstrates the use of the \a decluniupp function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = decluniupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF >     // Lower flag
inline decltype(auto) decluniupp( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,false>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   return decluniupp( declupp( *dm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-strictly-upper matrix multiplication expression as strictly upper.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid strictly upper matrix specification.
//
// The \a declstrupp function declares the given non-strictly-upper matrix multiplication
// expression \a dm as strictly upper. The function returns an expression representing the
// operation. In case the given expression does not represent a square matrix, a
// \a std::invalid_argument exception is thrown.\n
// The following example demonstrates the use of the \a declstrupp function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = declstrupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF >     // Lower flag
inline decltype(auto) declstrupp( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,false>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   return declstrupp( declupp( *dm ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-diagonal matrix multiplication expression as diagonal.
// \ingroup dense_matrix
//
// \param dm The input matrix multiplication expression.
// \return The redeclared dense matrix multiplication expression.
// \exception std::invalid_argument Invalid diagonal matrix specification.
//
// The \a decldiag function declares the given non-diagonal matrix multiplication expression
// \a dm as diagonal. The function returns an expression representing the operation. In case
// the given expression does not represent a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a decldiag function:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, B, C;
   // ... Resizing and initialization
   C = decldiag( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
inline decltype(auto) decldiag( const TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid diagonal matrix specification" );
   }

   using ReturnType = const TDMatTDMatMultExpr<MT1,MT2,SF,HF,true,true>;
   return ReturnType( dm.leftOperand(), dm.rightOperand() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct Size< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 0UL >
   : public Size<MT1,0UL>
{};

template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct Size< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 1UL >
   : public Size<MT2,1UL>
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
template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct IsAligned< TDMatTDMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public BoolConstant< IsAligned_v<MT1> && IsAligned_v<MT2> >
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
