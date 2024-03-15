//=================================================================================================
/*!
//  \file blaze/math/expressions/TDMatTSMatMultExpr.h
//  \brief Header file for the transpose dense matrix/transpose sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TDMATTSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TDMATTSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/functors/DeclDiag.h>
#include <blaze/math/functors/DeclHerm.h>
#include <blaze/math/functors/DeclLow.h>
#include <blaze/math/functors/DeclSym.h>
#include <blaze/math/functors/DeclUpp.h>
#include <blaze/math/functors/Noop.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
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
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TDMATTSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose dense matrix-transpose sparse matrix multiplications.
// \ingroup dense_matrix_expression
//
// The TDMatTSMatMultExpr class represents the compile time expression for multiplications between
// a column-major dense matrix and a column-major sparse matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF       // Lower flag
        , bool UF >     // Upper flag
class TDMatTSMatMultExpr
   : public MatMatMultExpr< DenseMatrix< TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side dense matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side sparse matrix expression.
   using ET1 = ElementType_t<RT1>;    //!< Element type of the left-hand side dense matrix expression.
   using ET2 = ElementType_t<RT2>;    //!< Element type of the right-hand side sparse matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side dense matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side dense matrix expression.
   static constexpr bool evaluateLeft = ( IsComputation_v<MT1> || RequiresEvaluation_v<MT1> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
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
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either of the two matrix operands requires an intermediate evaluation, the variable
       will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool IsEvaluationRequired_v = ( evaluateLeft || evaluateRight );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case all three involved data types are suited for a vectorized computation of the
       matrix multiplication, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseVectorizedKernel_v =
      ( useOptimizedKernels &&
        !IsDiagonal_v<T2> &&
        T1::simdEnabled && T2::simdEnabled &&
        IsColumnMajorMatrix_v<T1> &&
        IsSIMDCombinable_v< ElementType_t<T1>
                          , ElementType_t<T2>
                          , ElementType_t<T3> > &&
        HasSIMDAdd_v< ElementType_t<T2>, ElementType_t<T3> > &&
        HasSIMDMult_v< ElementType_t<T2>, ElementType_t<T3> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case a vectorized computation of the matrix multiplication is not possible, but a
       loop-unrolled computation is feasible, the variable will be set to 1, otherwise it will
       be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseOptimizedKernel_v =
      ( useOptimizedKernels &&
        !UseVectorizedKernel_v<T1,T2,T3> &&
        !IsDiagonal_v<T2> &&
        !IsResizable_v< ElementType_t<T1> > &&
        !IsResizable_v<ET2> );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! In case neither a vectorized nor optimized computation is possible, the variable will be
       set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool UseDefaultKernel_v =
      ( !UseVectorizedKernel_v<T1,T2,T3> && !UseOptimizedKernel_v<T1,T2,T3> );
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
   //! Type of this TDMatTSMatMultExpr instance.
   using This = TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>;

   //! Base type of this TDMatTSMatMultExpr instance.
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

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;

   //! Type for the assignment of the left-hand side dense matrix operand.
   using LT = If_t< evaluateLeft, const RT1, CT1 >;

   //! Type for the assignment of the right-hand side sparse matrix operand.
   using RT = If_t< evaluateRight, const RT2, CT2 >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled =
      ( !IsDiagonal_v<MT1> &&
        MT1::simdEnabled &&
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
   /*!\brief Constructor for the TDMatTSMatMultExpr class.
   //
   // \param lhs The left-hand side dense matrix operand of the multiplication expression.
   // \param rhs The right-hand side sparse matrix operand of the multiplication expression.
   */
   inline TDMatTSMatMultExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side dense matrix of the multiplication expression
      , rhs_( rhs )  // Right-hand side sparse matrix of the multiplication expression
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
   /*!\brief Returns the right-hand side transpose sparse matrix operand.
   //
   // \return The right-hand side transpose sparse matrix operand.
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
      return lhs_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( rows() * columns() >= SMP_TDMATTSMATMULT_THRESHOLD ) && !IsDiagonal_v<MT1>;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side dense matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-transpose sparse matrix multiplication to a
   //        dense matrix (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTSMatMultExpr::selectAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose dense matrix-transpose sparse matrix multiplication
   //        to dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose dense matrix-
   // transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> || IsDiagonal_v<MT4> ? A.rows() : 64UL );

      reset( C );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            auto element( B.begin(j) );
            const auto end( B.end(j) );

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               if( IsDiagonal_v<MT4> )
               {
                  C(j1,j) = A(j1,j1) * element->value();
               }
               else
               {
                  const size_t ibegin( ( IsLower_v<MT4> )
                                       ?( ( LOW )
                                          ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                          :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                       :( LOW ? max(j,ii) : ii ) );
                  const size_t iend( ( IsUpper_v<MT4> )
                                     ?( ( SYM || HERM || UPP )
                                        ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                        :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                     :( SYM || HERM || UPP ? min(j+1UL,itmp) : itmp ) );

                  if( ( SYM || HERM || LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                     continue;

                  BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

                  for( size_t i=ibegin; i<iend; ++i ) {
                     if( isDefault( C(i,j) ) )
                        C(i,j) = A(i,j1) * element->value();
                     else
                        C(i,j) += A(i,j1) * element->value();
                  }
               }
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=0UL; j<B.columns(); ++j ) {
            for( size_t i=j+1UL; i<A.rows(); ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized assignment to dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized assignment of a transpose dense matrix-transpose sparse matrix multiplication
   //        to dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the optimized assignment kernel for the transpose dense matrix-
   // transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> ? A.rows() : 64UL );

      reset( C );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            const auto end( B.end(j) );
            auto element( B.begin(j) );

            const size_t nonzeros( B.nonZeros(j) );
            const size_t kpos( prevMultiple( nonzeros, 4UL ) );
            BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

            for( size_t k=0UL; k<kpos; k+=4UL )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );
               ++element;
               const size_t j2( element->index() );
               const ET2    v2( element->value() );
               ++element;
               const size_t j3( element->index() );
               const ET2    v3( element->value() );
               ++element;
               const size_t j4( element->index() );
               const ET2    v4( element->value() );
               ++element;

               BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( SYM || HERM || UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) ) )
                                  :( SYM || HERM || UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( SYM || HERM || LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) += A(i    ,j1) * v1 + A(i    ,j2) * v2 + A(i    ,j3) * v3 + A(i    ,j4) * v4;
                  C(i+1UL,j) += A(i+1UL,j1) * v1 + A(i+1UL,j2) * v2 + A(i+1UL,j3) * v3 + A(i+1UL,j4) * v4;
                  C(i+2UL,j) += A(i+2UL,j1) * v1 + A(i+2UL,j2) * v2 + A(i+2UL,j3) * v3 + A(i+2UL,j4) * v4;
                  C(i+3UL,j) += A(i+3UL,j1) * v1 + A(i+3UL,j2) * v2 + A(i+3UL,j3) * v3 + A(i+3UL,j4) * v4;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) += A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
               }
            }

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( SYM || HERM || UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                  :( SYM || HERM || UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( SYM || HERM || LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) += A(i    ,j1) * v1;
                  C(i+1UL,j) += A(i+1UL,j1) * v1;
                  C(i+2UL,j) += A(i+2UL,j1) * v1;
                  C(i+3UL,j) += A(i+3UL,j1) * v1;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) += A(i,j1) * v1;
               }
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=0UL; j<B.columns(); ++j ) {
            for( size_t i=j+1UL; i<A.rows(); ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized assignment to column-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to column-major dense matrices (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the vectorized column-major assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      reset( C );

      for( size_t j=0UL; j<B.columns(); ++j )
      {
         const auto end( B.end(j) );
         auto element( B.begin(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( prevMultiple( nonzeros, 4UL ) );
         BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const SIMDType xmm1( set( v1 ) );
            const SIMDType xmm2( set( v2 ) );
            const SIMDType xmm3( set( v3 ) );
            const SIMDType xmm4( set( v4 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( SYM || HERM || UPP ? max(j+1UL,j4) : j4 )
                                  :( SYM || HERM || UPP ? max(j,j4)+1UL : j4+1UL ) )
                               :( SYM || HERM || UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) + A.load(i,j1) * xmm1 + A.load(i,j2) * xmm2 + A.load(i,j3) * xmm3 + A.load(i,j4) * xmm4 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) += A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const SIMDType xmm1( set( v1 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( SYM || HERM || UPP ? max(j+1UL,j1) : j1 )
                                  :( SYM || HERM || UPP ? max(j,j1)+1UL : j1+1UL ) )
                               :( SYM || HERM || UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) + A.load(i,j1) * xmm1 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) += A(i,j1) * v1;
            }
         }
      }

      if( SYM || HERM ) {
         for( size_t j=0UL; j<B.columns(); ++j ) {
            for( size_t i=j+1UL; i<A.rows(); ++i ) {
               C(i,j) = HERM ? conj( C(j,i) ) : C(j,i);
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose dense matrix-transpose sparse matrix multiplication to a
   //        sparse matrix (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose dense matrix-
   // transpose sparse matrix multiplication expression to a sparse matrix.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
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

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose dense matrix-transpose sparse matrix multiplication
   //        to a dense matrix (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose dense
   // matrix-transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTSMatMultExpr::selectAddAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices***********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> || IsDiagonal_v<MT4> ? A.rows() : 64UL );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            auto element( B.begin(j) );
            const auto end( B.end(j) );

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               if( IsDiagonal_v<MT4> )
               {
                  C(j1,j) += A(j1,j1) * element->value();
               }
               else
               {
                  const size_t ibegin( ( IsLower_v<MT4> )
                                       ?( ( LOW )
                                          ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                          :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                       :( LOW ? max(j,ii) : ii ) );
                  const size_t iend( ( IsUpper_v<MT4> )
                                     ?( ( UPP )
                                        ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                        :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                     :( UPP ? min(j+1UL,itmp) : itmp ) );

                  if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                     continue;

                  BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

                  const size_t inum( iend - ibegin );
                  const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
                  BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

                  for( size_t i=ibegin; i<ipos; i+=4UL ) {
                     C(i    ,j) += A(i    ,j1) * element->value();
                     C(i+1UL,j) += A(i+1UL,j1) * element->value();
                     C(i+2UL,j) += A(i+2UL,j1) * element->value();
                     C(i+3UL,j) += A(i+3UL,j1) * element->value();
                  }
                  for( size_t i=ipos; i<iend; ++i ) {
                     C(i,j) += A(i,j1) * element->value();
                  }
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized addition assignment to dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized addition assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the optimized addition assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> ? A.rows() : 64UL );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            const auto end( B.end(j) );
            auto element( B.begin(j) );

            const size_t nonzeros( B.nonZeros(j) );
            const size_t kpos( prevMultiple( nonzeros, 4UL ) );
            BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

            for( size_t k=0UL; k<kpos; k+=4UL )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );
               ++element;
               const size_t j2( element->index() );
               const ET2    v2( element->value() );
               ++element;
               const size_t j3( element->index() );
               const ET2    v3( element->value() );
               ++element;
               const size_t j4( element->index() );
               const ET2    v4( element->value() );
               ++element;

               BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) ) )
                                  :( UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) += A(i    ,j1) * v1 + A(i    ,j2) * v2 + A(i    ,j3) * v3 + A(i    ,j4) * v4;
                  C(i+1UL,j) += A(i+1UL,j1) * v1 + A(i+1UL,j2) * v2 + A(i+1UL,j3) * v3 + A(i+1UL,j4) * v4;
                  C(i+2UL,j) += A(i+2UL,j1) * v1 + A(i+2UL,j2) * v2 + A(i+2UL,j3) * v3 + A(i+2UL,j4) * v4;
                  C(i+3UL,j) += A(i+3UL,j1) * v1 + A(i+3UL,j2) * v2 + A(i+3UL,j3) * v3 + A(i+3UL,j4) * v4;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) += A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
               }
            }

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                  :( UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) += A(i    ,j1) * v1;
                  C(i+1UL,j) += A(i+1UL,j1) * v1;
                  C(i+2UL,j) += A(i+2UL,j1) * v1;
                  C(i+3UL,j) += A(i+3UL,j1) * v1;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) += A(i,j1) * v1;
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized addition assignment to column-major dense matrices*******************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized addition assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to column-major dense matrices (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the vectorized addition assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      for( size_t j=0UL; j<B.columns(); ++j )
      {
         const auto end( B.end(j) );
         auto element( B.begin(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( prevMultiple( nonzeros, 4UL ) );
         BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const SIMDType xmm1( set( v1 ) );
            const SIMDType xmm2( set( v2 ) );
            const SIMDType xmm3( set( v3 ) );
            const SIMDType xmm4( set( v4 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? max(j+1UL,j4) : j4 )
                                  :( UPP ? max(j,j4)+1UL : j4+1UL ) )
                               :( UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) + A.load(i,j1) * xmm1 + A.load(i,j2) * xmm2 + A.load(i,j3) * xmm3 + A.load(i,j4) * xmm4 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) += A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const SIMDType xmm1( set( v1 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? max(j+1UL,j1) : j1 )
                                  :( UPP ? max(j,j1)+1UL : j1+1UL ) )
                               :( UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) + A.load(i,j1) * xmm1 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) += A(i,j1) * v1;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // dense matrix-transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side dense matrix operand
      RT B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TDMatTSMatMultExpr::selectSubAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to matrices (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseDefaultKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> || IsDiagonal_v<MT4> ? A.rows() : 64UL );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            auto element( B.begin(j) );
            const auto end( B.end(j) );

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );

               if( IsDiagonal_v<MT4> )
               {
                  C(j1,j) -= A(j1,j1) * element->value();
               }
               else
               {
                  const size_t ibegin( ( IsLower_v<MT4> )
                                       ?( ( LOW )
                                          ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                          :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                       :( LOW ? max(j,ii) : ii ) );
                  const size_t iend( ( IsUpper_v<MT4> )
                                     ?( ( UPP )
                                        ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                        :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                     :( UPP ? min(j+1UL,itmp) : itmp ) );

                  if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                     continue;

                  BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

                  const size_t inum( iend - ibegin );
                  const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
                  BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

                  for( size_t i=ibegin; i<ipos; i+=4UL ) {
                     C(i    ,j) -= A(i    ,j1) * element->value();
                     C(i+1UL,j) -= A(i+1UL,j1) * element->value();
                     C(i+2UL,j) -= A(i+2UL,j1) * element->value();
                     C(i+3UL,j) -= A(i+3UL,j1) * element->value();
                  }
                  for( size_t i=ipos; i<iend; ++i ) {
                     C(i,j) -= A(i,j1) * element->value();
                  }
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Optimized subtraction assignment to dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Optimized subtraction assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to matrices (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the optimized subtraction assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseOptimizedKernel_v<MT3,MT4,MT5> >
   {
      const size_t block( IsColumnMajorMatrix_v<MT3> ? A.rows() : 64UL );

      for( size_t ii=0UL; ii<A.rows(); ii+=block )
      {
         const size_t itmp( min( ii+block, A.rows() ) );

         for( size_t j=0UL; j<B.columns(); ++j )
         {
            const auto end( B.end(j) );
            auto element( B.begin(j) );

            const size_t nonzeros( B.nonZeros(j) );
            const size_t kpos( prevMultiple( nonzeros, 4UL ) );
            BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

            for( size_t k=0UL; k<kpos; k+=4UL )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );
               ++element;
               const size_t j2( element->index() );
               const ET2    v2( element->value() );
               ++element;
               const size_t j3( element->index() );
               const ET2    v3( element->value() );
               ++element;
               const size_t j4( element->index() );
               const ET2    v4( element->value() );
               ++element;

               BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j4 : j4+1UL ) ) ) )
                                  :( UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) -= A(i    ,j1) * v1 + A(i    ,j2) * v2 + A(i    ,j3) * v3 + A(i    ,j4) * v4;
                  C(i+1UL,j) -= A(i+1UL,j1) * v1 + A(i+1UL,j2) * v2 + A(i+1UL,j3) * v3 + A(i+1UL,j4) * v4;
                  C(i+2UL,j) -= A(i+2UL,j1) * v1 + A(i+2UL,j2) * v2 + A(i+2UL,j3) * v3 + A(i+2UL,j4) * v4;
                  C(i+3UL,j) -= A(i+3UL,j1) * v1 + A(i+3UL,j2) * v2 + A(i+3UL,j3) * v3 + A(i+3UL,j4) * v4;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) -= A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
               }
            }

            for( ; element!=end; ++element )
            {
               const size_t j1( element->index() );
               const ET2    v1( element->value() );

               const size_t ibegin( ( IsLower_v<MT4> )
                                    ?( ( LOW )
                                       ?( max( j, ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) )
                                       :( max( ii, ( IsStrictlyLower_v<MT4> ? j1+1UL : j1 ) ) ) )
                                    :( LOW ? max(j,ii) : ii ) );
               const size_t iend( ( IsUpper_v<MT4> )
                                  ?( ( UPP )
                                     ?( min( j+1UL, itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) )
                                     :( min( itmp, ( IsStrictlyUpper_v<MT4> ? j1 : j1+1UL ) ) ) )
                                  :( UPP ? min(j+1UL,itmp) : itmp ) );

               if( ( LOW || UPP || IsTriangular_v<MT4> ) && ( ibegin >= iend ) )
                  continue;

               BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

               const size_t inum( iend - ibegin );
               const size_t ipos( ibegin + prevMultiple( inum, 4UL ) );
               BLAZE_INTERNAL_ASSERT( ipos <= ibegin+inum, "Invalid end calculation" );

               for( size_t i=ibegin; i<ipos; i+=4UL ) {
                  C(i    ,j) -= A(i    ,j1) * v1;
                  C(i+1UL,j) -= A(i+1UL,j1) * v1;
                  C(i+2UL,j) -= A(i+2UL,j1) * v1;
                  C(i+3UL,j) -= A(i+3UL,j1) * v1;
               }
               for( size_t i=ipos; i<iend; ++i ) {
                  C(i,j) -= A(i,j1) * v1;
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Vectorized subtraction assignment to column-major dense matrices****************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Vectorized subtraction assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to column-major matrices (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side dense matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the vectorized subtraction assignment kernel for the transpose dense
   // matrix-transpose sparse matrix multiplication to a column-major dense matrix.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline auto selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
      -> EnableIf_t< UseVectorizedKernel_v<MT3,MT4,MT5> >
   {
      constexpr bool remainder( !IsPadded_v<MT3> || !IsPadded_v<MT4> );

      for( size_t j=0UL; j<B.columns(); ++j )
      {
         const auto end( B.end(j) );
         auto element( B.begin(j) );

         const size_t nonzeros( B.nonZeros(j) );
         const size_t kpos( prevMultiple( nonzeros, 4UL ) );
         BLAZE_INTERNAL_ASSERT( kpos <= nonzeros, "Invalid end calculation" );

         for( size_t k=0UL; k<kpos; k+=4UL )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );
            ++element;
            const size_t j2( element->index() );
            const ET2    v2( element->value() );
            ++element;
            const size_t j3( element->index() );
            const ET2    v3( element->value() );
            ++element;
            const size_t j4( element->index() );
            const ET2    v4( element->value() );
            ++element;

            BLAZE_INTERNAL_ASSERT( j1 < j2 && j2 < j3 && j3 < j4, "Invalid sparse matrix index detected" );

            const SIMDType xmm1( set( v1 ) );
            const SIMDType xmm2( set( v2 ) );
            const SIMDType xmm3( set( v3 ) );
            const SIMDType xmm4( set( v4 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? max(j+1UL,j4) : j4 )
                                  :( UPP ? max(j,j4)+1UL : j4+1UL ) )
                               :( UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) - A.load(i,j1) * xmm1 - A.load(i,j2) * xmm2 - A.load(i,j3) * xmm3 - A.load(i,j4) * xmm4 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) -= A(i,j1) * v1 + A(i,j2) * v2 + A(i,j3) * v3 + A(i,j4) * v4;
            }
         }

         for( ; element!=end; ++element )
         {
            const size_t j1( element->index() );
            const ET2    v1( element->value() );

            const SIMDType xmm1( set( v1 ) );

            const size_t ibegin( ( IsLower_v<MT4> )
                                 ?( ( IsStrictlyLower_v<MT4> )
                                    ?( prevMultiple( ( LOW ? max(j,j1+1UL) : j1+1UL ), SIMDSIZE ) )
                                    :( prevMultiple( ( LOW ? max(j,j1) : j1 ), SIMDSIZE ) ) )
                                 :( LOW ? prevMultiple( j, SIMDSIZE ) : 0UL ) );
            const size_t iend( ( IsUpper_v<MT4> )
                               ?( ( IsStrictlyUpper_v<MT4> )
                                  ?( UPP ? max(j+1UL,j1) : j1 )
                                  :( UPP ? max(j,j1)+1UL : j1+1UL ) )
                               :( UPP ? j+1UL : A.rows() ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
            BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

            size_t i( ibegin );

            for( ; i<ipos; i+=SIMDSIZE ) {
               C.store( i, j, C.load(i,j) - A.load(i,j1) * xmm1 );
            }
            for( ; remainder && i<iend; ++i ) {
               C(i,j) -= A(i,j1) * v1;
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ A\circ=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // dense matrix-transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
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
   /*!\brief SMP assignment of a transpose dense matrix-transpose sparse matrix multiplication
   //        to a dense matrix (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-transpose sparse matrix multiplication expression to a dense matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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
   /*!\brief SMP assignment of a transpose dense matrix-transpose sparse matrix multiplication
   //        to a sparse matrix (\f$ A=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose dense
   // matrix-transpose sparse matrix multiplication expression to a sparse matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target sparse matrix
           , bool SO >    // Storage order of the target sparse matrix
   friend inline auto smpAssign( SparseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
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

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ A+=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // dense matrix-transpose sparse matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ A-=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a transpose
   // dense matrix-transpose sparse matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      LT A( rhs.lhs_ );  // Evaluation of the left-hand side dense matrix operand
      RT B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a transpose dense matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ A\circ=B*C \f$).
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // transpose dense matrix-transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const TDMatTSMatMultExpr& rhs )
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT2 );
   BLAZE_CONSTRAINT_MUST_FORM_VALID_MATMATMULTEXPR( MT1, MT2 );
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
/*!\brief Backend implementation of the multiplication between a column-major dense matrix and a
//        column-major sparse matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The product of the two matrices.
//
// This function implements a performance optimized treatment of the multiplication between a
// column-major dense matrix and a column-major sparse matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , DisableIf_t< ( IsIdentity_v<MT2> &&
                         IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > ) ||
                       IsZero_v<MT2> >* = nullptr >
inline const TDMatTSMatMultExpr<MT1,MT2,false,false,false,false>
   tdmattsmatmult( const DenseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   return TDMatTSMatMultExpr<MT1,MT2,false,false,false,false>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication between a column-major dense matrix and a
//        column-major identity matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return Reference to the left-hand side dense matrix.
//
// This function implements a performance optimized treatment of the multiplication between
// a column-major dense matrix and a column-major identity matrix. It returns a reference to
// the left-hand side dense matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsIdentity_v<MT2> &&
                      IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > >* = nullptr >
inline const MT1&
   tdmattsmatmult( const DenseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   return (*lhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication between a column-major dense matrix and a
//        column-major zero matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side zero matrix for the multiplication.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the multiplication between a
// column-major dense matrix and a column-major zero matrix. It returns a zero matrix.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsZero_v<MT2> >* = nullptr >
inline decltype(auto)
   tdmattsmatmult( const DenseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   using ReturnType = const MultTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_ZERO_TYPE( ReturnType );

   return ReturnType( (*lhs).rows(), (*rhs).columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the multiplication of a column-major dense matrix and a
//        column-major sparse matrix (\f$ A=B*C \f$).
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The resulting matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of a column-major dense matrix and a column-major
// sparse matrix:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an expression representing a dense matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given matrices don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const DenseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).columns() != (*rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return tdmattsmatmult( *lhs, *rhs );
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
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
inline decltype(auto) declsym( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid symmetric matrix specification" );
   }

   using ReturnType = const TDMatTSMatMultExpr<MT1,MT2,true,HF,LF,UF>;
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
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
inline decltype(auto) declherm( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid Hermitian matrix specification" );
   }

   using ReturnType = const TDMatTSMatMultExpr<MT1,MT2,SF,true,LF,UF>;
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
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
inline decltype(auto) decllow( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   using ReturnType = const TDMatTSMatMultExpr<MT1,MT2,SF,HF,true,UF>;
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = declunilow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool UF >     // Upper flag
inline decltype(auto) declunilow( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,false,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid unilower matrix specification" );
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = declstrlow( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool UF >     // Upper flag
inline decltype(auto) declstrlow( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,false,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid strictly lower matrix specification" );
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
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
inline decltype(auto) declupp( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   using ReturnType = const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,true>;
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = decluniupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF >     // Lower flag
inline decltype(auto) decluniupp( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,false>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uniupper matrix specification" );
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
   // ... Resizing and initialization
   C = declstrupp( A * B );
   \endcode
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SF       // Symmetry flag
        , bool HF       // Hermitian flag
        , bool LF >     // Lower flag
inline decltype(auto) declstrupp( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,false>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid strictly upper matrix specification" );
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

   blaze::DynamicMatrix<double,columnMajor> A, C;
   blaze::CompressedMatrix<double,columnMajor> B;
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
inline decltype(auto) decldiag( const TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>& dm )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid diagonal matrix specification" );
   }

   using ReturnType = const TDMatTSMatMultExpr<MT1,MT2,SF,HF,true,true>;
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
struct Size< TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 0UL >
   : public Size<MT1,0UL>
{};

template< typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF >
struct Size< TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF>, 1UL >
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
struct IsAligned< TDMatTSMatMultExpr<MT1,MT2,SF,HF,LF,UF> >
   : public IsAligned<MT1>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
