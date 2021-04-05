//=================================================================================================
/*!
//  \file blaze/math/expressions/TSMatTSMatMultExpr.h
//  \brief Header file for the transpose sparse matrix/transpose sparse matrix multiplication expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_TSMATTSMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_TSMATTSMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Identity.h>
#include <blaze/math/constraints/MatMatMultExpr.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Zero.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Infinity.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TSMATTSMATMULTEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for transpose sparse matrix-transpose sparse matrix multiplications.
// \ingroup sparse_matrix_expression
//
// The TSMatTSMatMultExpr class represents the compile time expression for multiplications between
// two column-major sparse matrices.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
class TSMatTSMatMultExpr
   : public MatMatMultExpr< SparseMatrix< TSMatTSMatMultExpr<MT1,MT2>, true > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT1 = ResultType_t<MT1>;     //!< Result type of the left-hand side sparse matrix expression.
   using RT2 = ResultType_t<MT2>;     //!< Result type of the right-hand side sparse matrix expression.
   using CT1 = CompositeType_t<MT1>;  //!< Composite type of the left-hand side sparse matrix expression.
   using CT2 = CompositeType_t<MT2>;  //!< Composite type of the right-hand side sparse matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the left-hand side sparse matrix expression.
   static constexpr bool evaluateLeft = RequiresEvaluation_v<MT1>;
   //**********************************************************************************************

   //**********************************************************************************************
   //! Compilation switch for the composite type of the right-hand side sparse matrix expression.
   static constexpr bool evaluateRight = RequiresEvaluation_v<MT2>;
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the optimal evaluation strategy.
       In case the target matrix is row-major and both matrix operands are symmetric, the variable
       is set to 1 and an optimized evaluation strategy is selected. Otherwise the variable is set
       to 0 and the default strategy is chosen. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool CanExploitSymmetry_v =
      ( IsRowMajorMatrix_v<T1> && IsSymmetric_v<T2> && IsSymmetric_v<T3> );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   /*! This variable template is a helper for the selection of the parallel evaluation strategy.
       In case either of the two matrix operands requires an intermediate evaluation and in case
       no symmetry can be exploited, the variable will be set to 1, otherwise it will be 0. */
   template< typename T1, typename T2, typename T3 >
   static constexpr bool IsEvaluationRequired_v =
      ( ( evaluateLeft || evaluateRight ) && !CanExploitSymmetry_v<T1,T2,T3> );
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this TSMatTSMatMultExpr instance.
   using This = TSMatTSMatMultExpr<MT1,MT2>;

   //! Base type of this TSMatTSMatMultExpr instance.
   using BaseType = MatMatMultExpr< SparseMatrix<This,true> >;

   using ResultType    = MultTrait_t<RT1,RT2>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<ResultType>;    //!< Resulting element type.
   using ReturnType    = const ElementType;            //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite type of the left-hand side sparse matrix expression.
   using LeftOperand = If_t< IsExpression_v<MT1>, const MT1, const MT1& >;

   //! Composite type of the right-hand side sparse matrix expression.
   using RightOperand = If_t< IsExpression_v<MT2>, const MT2, const MT2& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable =
      ( !evaluateLeft  && MT1::smpAssignable && !evaluateRight && MT2::smpAssignable );
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the TSMatTSMatMultExpr class.
   //
   // \param lhs The left-hand side sparse matrix operand of the multiplication expression.
   // \param rhs The right-hand side sparse matrix operand of the multiplication expression.
   */
   inline TSMatTSMatMultExpr( const MT1& lhs, const MT2& rhs ) noexcept
      : lhs_( lhs )  // Left-hand side sparse matrix of the multiplication expression
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

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the sparse matrix.
   //
   // \return The number of non-zero elements in the sparse matrix.
   */
   constexpr size_t nonZeros() const noexcept {
      return 0UL;
   }
   //**********************************************************************************************

   //**NonZeros function***************************************************************************
   /*!\brief Returns the number of non-zero elements in the specified row.
   //
   // \param i The index of the row.
   // \return The number of non-zero elements of row \a i.
   */
   inline size_t nonZeros( size_t i ) const noexcept {
      MAYBE_UNUSED( i );
      return 0UL;
   }
   //**********************************************************************************************

   //**Left operand access*************************************************************************
   /*!\brief Returns the left-hand side transpose sparse matrix operand.
   //
   // \return The left-hand side transpose sparse matrix operand.
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
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return ( rows() * columns() >= SMP_TSMATTSMATMULT_THRESHOLD );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   LeftOperand  lhs_;  //!< Left-hand side sparse matrix of the multiplication expression.
   RightOperand rhs_;  //!< Right-hand side sparse matrix of the multiplication expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-transpose sparse matrix multiplication
   //        to a dense matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matrix
   friend inline auto assign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TSMatTSMatMultExpr::selectAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default assignment to dense matrices********************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default assignment of a transpose sparse matrix-transpose sparse matrix multiplication
   //        (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default assignment kernel for the transpose sparse matrix-
   // transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      for( size_t j=0UL; j<C.columns(); ++j ) {
         const auto rend( B.end(j) );
         for( auto relem=B.begin(j); relem!=rend; ++relem ) {
            const auto lend( A.end( relem->index() ) );
            for( auto lelem=A.begin( relem->index() ); lelem!=lend; ++lelem )
            {
               if( IsResizable_v< ElementType_t<MT3> > &&
                   isDefault( C(lelem->index(),j) ) ) {
                  C(lelem->index(),j) = lelem->value() * relem->value();
               }
               else {
                  C(lelem->index(),j) += lelem->value() * relem->value();
               }
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-transpose sparse matrix multiplication to a
   //        row-major sparse matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose sparse matrix multiplication expression to a row-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline auto assign( SparseMatrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ResultType );

      const ResultType tmp( serial( rhs ) );
      (*lhs).reserve( tmp.nonZeros() );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to column-major sparse matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a transpose sparse matrix-transpose sparse matrix multiplication to a
   //        column-major sparse matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a transpose sparse matrix-
   // transpose sparse matrix multiplication expression to a column-major sparse matrix.
   */
   template< typename MT >  // Type of the target sparse matrix
   friend inline void assign( SparseMatrix<MT,true>& lhs, const TSMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      // (Over-)Estimating the number of non-zero entries in the resulting matrix
      size_t nonzeros( 0UL );

      for( size_t j=0UL; j<(*lhs).columns(); ++j ) {
         const auto rend( B.end(j) );
         for( auto relem=B.begin(j); relem!=rend; ++relem ) {
            nonzeros += A.nonZeros( relem->index() );
         }
      }

      if( nonzeros > (*lhs).rows() * (*lhs).columns() ) {
         nonzeros = (*lhs).rows() * (*lhs).columns();
      }

      (*lhs).reserve( nonzeros );
      nonzeros = 0UL;

      // Performing the matrix-matrix multiplication
      SmallArray<ElementType,128UL> values ( (*lhs).rows(), ElementType() );
      SmallArray<bool,128UL>        valid  ( (*lhs).rows(), false );
      SmallArray<size_t,128UL>      indices( (*lhs).rows(), 0UL );
      size_t minIndex( inf ), maxIndex( 0UL );

      for( size_t j=0UL; j<(*lhs).columns(); ++j )
      {
         const auto rend( B.end(j) );
         for( auto relem=B.begin(j); relem!=rend; ++relem )
         {
            const auto lend( A.end( relem->index() ) );
            for( auto lelem=A.begin( relem->index() ); lelem!=lend; ++lelem )
            {
               if( !valid[lelem->index()] ) {
                  values[lelem->index()] = lelem->value() * relem->value();
                  valid [lelem->index()] = true;
                  indices[nonzeros] = lelem->index();
                  ++nonzeros;
                  if( lelem->index() < minIndex ) minIndex = lelem->index();
                  if( lelem->index() > maxIndex ) maxIndex = lelem->index();
               }
               else {
                  values[lelem->index()] += lelem->value() * relem->value();
               }
            }
         }

         BLAZE_INTERNAL_ASSERT( nonzeros <= (*lhs).rows(), "Invalid number of non-zero elements" );

         if( nonzeros > 0UL )
         {
            BLAZE_INTERNAL_ASSERT( minIndex <= maxIndex, "Invalid index detected" );

            if( ( nonzeros + nonzeros ) < ( maxIndex - minIndex ) )
            {
               std::sort( indices.begin(), indices.begin() + nonzeros );

               for( size_t i=0UL; i<nonzeros; ++i )
               {
                  const size_t index( indices[i] );
                  if( !isDefault( values[index] ) ) {
                     (*lhs).append( index, j, values[index] );
                     reset( values[index] );
                  }

                  reset( valid[index] );
               }
            }
            else {
               for( size_t i=minIndex; i<=maxIndex; ++i )
               {
                  if( !isDefault( values[i] ) ) {
                     (*lhs).append( i, j, values[i] );
                     reset( values[i] );
                  }

                  reset( valid[i] );
               }
            }

            nonzeros = 0UL;
            minIndex = inf;
            maxIndex = 0UL;
         }

         (*lhs).finalize( j );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring assignment to row-major matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a row-major matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring assignment of a transpose sparse
   // matrix-transpose sparse matrix multiplication expression to a row-major matrix. Due to the
   // explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto assign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      assign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a transpose sparse
   // matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline auto addAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TSMatTSMatMultExpr::selectAddAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default addition assignment to dense matrices***********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default addition assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication (\f$ C+=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default addition assignment kernel for the transpose sparse
   // matrix-transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectAddAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      for( size_t j=0UL; j<C.columns(); ++j ) {
         const auto rend( B.end(j) );
         for( auto relem=B.begin(j); relem!=rend; ++relem ) {
            const auto lend( A.end( relem->index() ) );
            for( auto lelem=A.begin( relem->index() ); lelem!=lend; ++lelem ) {
               C(lelem->index(),j) += lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring addition assignment to row-major matrices*************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring addition assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a row-major matrix (\f$ C+=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring addition assignment of a transpose
   // sparse matrix-transpose sparse matrix multiplication expression to a row-major matrix. Due
   // to the explicit application of the SFINAE principle this function can only be selected by
   // the compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto addAssign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      addAssign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a transpose
   // sparse matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline auto subAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> DisableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( serial( rhs.lhs_ ) );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( serial( rhs.rhs_ ) );  // Evaluation of the right-hand side sparse matrix operand

      BLAZE_INTERNAL_ASSERT( A.rows()    == rhs.lhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( A.columns() == rhs.lhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( B.rows()    == rhs.rhs_.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == rhs.rhs_.columns(), "Invalid number of columns" );
      BLAZE_INTERNAL_ASSERT( A.rows()    == (*lhs).rows()     , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( B.columns() == (*lhs).columns()  , "Invalid number of columns" );

      TSMatTSMatMultExpr::selectSubAssignKernel( *lhs, A, B );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Default subtraction assignment to dense matrices********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Default subtraction assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication (\f$ C-=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param C The target left-hand side dense matrix.
   // \param A The left-hand side sparse matrix operand.
   // \param B The right-hand side sparse matrix operand.
   // \return void
   //
   // This function implements the default subtraction assignment kernel for the transpose sparse
   // matrix-transpose sparse matrix multiplication.
   */
   template< typename MT3    // Type of the left-hand side target matrix
           , typename MT4    // Type of the left-hand side matrix operand
           , typename MT5 >  // Type of the right-hand side matrix operand
   static inline void selectSubAssignKernel( MT3& C, const MT4& A, const MT5& B )
   {
      for( size_t j=0UL; j<C.columns(); ++j ) {
         const auto rend( B.end(j) );
         for( auto relem=B.begin(j); relem!=rend; ++relem ) {
            const auto lend( A.end( relem->index() ) );
            for( auto lelem=A.begin( relem->index() ); lelem!=lend; ++lelem ) {
               C(lelem->index(),j) -= lelem->value() * relem->value();
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Restructuring subtraction assignment to row-major matrices**********************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring subtraction assignment of a transpose sparse matrix-transpose sparse
   //        matrix multiplication to a row-major matrix (\f$ C-=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring subtraction assignment of a
   // transpose sparse matrix-transpose sparse matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only
   // be selected by the compiler in case the symmetry of either of the two matrix operands
   // can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto subAssign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      subAssign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a transpose
   // sparse matrix-sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline void schurAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
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

   //**SMP assignment to matrices******************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP assignment of a transpose sparse matrix-transpose sparse matrix multiplication
   //        to a matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized SMP assignment of a transpose sparse
   // matrix-transpose sparse matrix multiplication expression to a matrix. Due to the explicit
   // application of the SFINAE principle this function can only be selected by the compiler in
   // case either of the two matrix operands requires an intermediate evaluation and no symmetry
   // can be exploited.
   */
   template< typename MT  // Type of the target matrix
           , bool SO >    // Storage order of the target matrix
   friend inline auto smpAssign( Matrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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

   //**Restructuring SMP assignment to row-major matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Restructuring SMP assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a row-major matrix (\f$ C=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be assigned.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP assignment of a transpose
   // sparse matrix-transpose sparse matrix multiplication expression to a row-major matrix. Due
   // to the explicit application of the SFINAE principle this function can only be selected by
   // the compiler in case the symmetry of either of the two matrix operands can be exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAssign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAssign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to dense matrices***************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP addition assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C+=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the performance optimized SMP addition assignment of a transpose
   // sparse matrix-transpose sparse matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by
   // the compiler in case either of the two matrix operands requires an intermediate evaluation
   // and no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline auto smpAddAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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
   /*!\brief Restructuring SMP addition assignment of a transpose sparse matrix-transpose sparse
   //        matrix multiplication to a row-major matrix (\f$ C+=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be added.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP addition assignment of a
   // transpose sparse matrix-transpose sparse matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only be
   // selected by the compiler in case the symmetry of either of the two matrix operands can be
   // exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpAddAssign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpAddAssign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP addition assignment to sparse matrices**************************************************
   // No special implementation for the SMP addition assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP subtraction assignment to dense matrices************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP subtraction assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C-=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized SMP subtraction assignment of a transpose
   // sparse matrix-transpose sparse matrix multiplication expression to a dense matrix. Due to
   // the explicit application of the SFINAE principle this function can only be selected by the
   // compiler in case either of the two matrix operands requires an intermediate evaluation and
   // no symmetry can be exploited.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline auto smpSubAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< IsEvaluationRequired_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CT1 A( rhs.lhs_ );  // Evaluation of the left-hand side sparse matrix operand
      CT2 B( rhs.rhs_ );  // Evaluation of the right-hand side sparse matrix operand

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
   /*!\brief Restructuring SMP subtraction assignment of a transpose sparse matrix-transpose sparse
   //        matrix multiplication to a row-major matrix (\f$ C-=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side matrix.
   // \param rhs The right-hand side multiplication expression to be subtracted.
   // \return void
   //
   // This function implements the symmetry-based restructuring SMP subtraction assignment of a
   // transpose sparse matrix-transpose sparse matrix multiplication expression to a row-major
   // matrix. Due to the explicit application of the SFINAE principle this function can only be
   // selected by the compiler in case the symmetry of either of the two matrix operands can be
   // exploited.
   */
   template< typename MT >  // Type of the target matrix
   friend inline auto smpSubAssign( Matrix<MT,false>& lhs, const TSMatTSMatMultExpr& rhs )
      -> EnableIf_t< CanExploitSymmetry_v<MT,MT1,MT2> >
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      smpSubAssign( *lhs, trans( rhs.lhs_ ) * trans( rhs.rhs_ ) );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**SMP subtraction assignment to sparse matrices***********************************************
   // No special implementation for the SMP subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**SMP Schur product assignment to dense matrices**********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief SMP Schur product assignment of a transpose sparse matrix-transpose sparse matrix
   //        multiplication to a dense matrix (\f$ C\circ=A*B \f$).
   // \ingroup sparse_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side multiplication expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized SMP Schur product assignment of a
   // transpose sparse matrix-transpose sparse matrix multiplication expression to a dense matrix.
   */
   template< typename MT  // Type of the target dense matrix
           , bool SO >    // Storage order of the target dense matarix
   friend inline void smpSchurAssign( DenseMatrix<MT,SO>& lhs, const TSMatTSMatMultExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ResultType );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ZERO_TYPE( MT1 );
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
/*!\brief Backend implementation of the multiplication between a column-major sparse matrix and a
//        column-major sparse matrix (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The product of the two matrices.
//
// This function implements a performance optimized treatment of the multiplication between a
// column-major sparse matrix and a column-major spares matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , DisableIf_t< ( ( IsIdentity_v<MT1> || IsIdentity_v<MT2> ) &&
                         IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > ) ||
                       ( IsIdentity_v<MT1> && IsIdentity_v<MT2> ) ||
                       ( IsZero_v<MT1> || IsZero_v<MT2> ) >* = nullptr >
inline const TSMatTSMatMultExpr<MT1,MT2>
   tsmattsmatmult( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   return TSMatTSMatMultExpr<MT1,MT2>( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication between a column-major sparse matrix and
//        a column-major identity matrix (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return Reference to the left-hand side sparse matrix.
//
// This function implements a performance optimized treatment of the multiplication between
// a column-major sparse matrix and a column-major identity matrix. It returns a reference to
// the left-hand side sparse matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< !IsIdentity_v<MT1> && IsIdentity_v<MT2> &&
                      IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > >* = nullptr >
inline const MT1&
   tsmattsmatmult( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
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
/*!\brief Backend implementation of the multiplication between a column-major identity matrix
//        and a column-major sparse matrix (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side identity matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return Reference to the right-hand side sparse matrix.
//
// This function implements a performance optimized treatment of the multiplication between a
// column-major identity matrix and a column-major sparse matrix. It returns a reference to the
// right-hand side sparse matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsIdentity_v<MT1> && !IsIdentity_v<MT2> &&
                      IsSame_v< ElementType_t<MT1>, ElementType_t<MT2> > >* = nullptr >
inline const MT2&
   tsmattsmatmult( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( lhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   return (*rhs);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication between two column-major identity matrices
//        (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side identity matrix for the multiplication.
// \param rhs The right-hand side identity matrix for the multiplication.
// \return The resulting identity matrix.
//
// This function implements a performance optimized treatment of the multiplication between
// two column-major identity matrices. It returns an identity matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsIdentity_v<MT1> && IsIdentity_v<MT2> >* = nullptr >
inline decltype(auto)
   tsmattsmatmult( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( rhs );

   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).rows(), "Invalid matrix sizes" );

   using ReturnType = const MultTrait_t< ResultType_t<MT1>, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( ReturnType );
   BLAZE_CONSTRAINT_MUST_BE_IDENTITY_MATRIX_TYPE( ReturnType );

   return ReturnType( (*lhs).rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend implementation of the multiplication between a column-major (zero) sparse matrix
//        and a column-major (zero) sparse matrix (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the multiplication.
// \param rhs The right-hand side sparse matrix for the multiplication.
// \return The resulting zero matrix.
//
// This function implements a performance optimized treatment of the multiplication between a
// column-major (zero) sparse matrix and a column-major (zero) sparse matrix. It returns a zero
// matrix.
*/
template< typename MT1  // Type of the left-hand side sparse matrix
        , typename MT2  // Type of the right-hand side sparse matrix
        , EnableIf_t< IsZero_v<MT1> || IsZero_v<MT2> >* = nullptr >
inline decltype(auto)
   tsmattsmatmult( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
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
/*!\brief Multiplication operator for the multiplication of two column-major sparse matrices
//        (\f$ A=B*C \f$).
// \ingroup sparse_matrix
//
// \param lhs The left-hand side sparse matrix for the matrix multiplication.
// \param rhs The right-hand side sparse matrix for the matrix multiplication.
// \return The product of the two matrices.
// \exception std::invalid_argument Matrix sizes do not match.
//
// This operator represents the multiplication of two column-major sparse matrices:

   \code
   using blaze::columnMajor;

   blaze::CompressedMatrix<double,rowMajor> A, B, C;
   // ... Resizing and initialization
   C = A * B;
   \endcode

// The operator returns an expression representing a sparse matrix of the higher-order element
// type of the two involved matrix element types \a MT1::ElementType and \a MT2::ElementType.
// Both matrix types \a MT1 and \a MT2 as well as the two element types \a MT1::ElementType
// and \a MT2::ElementType have to be supported by the MultTrait class template.\n
// In case the current number of rows and columns of the two given matrices don't match, a
// \a std::invalid_argument is thrown.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const SparseMatrix<MT1,true>& lhs, const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   if( (*lhs).columns() != (*rhs).rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   return tsmattsmatmult( *lhs, *rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
