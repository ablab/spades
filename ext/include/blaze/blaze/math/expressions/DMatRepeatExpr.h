//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatRepeatExpr.h
//  \brief Header file for the dense matrix repeat expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATREPEATEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATREPEATEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatRepeatExpr.h>
#include <blaze/math/expressions/RepeatExprData.h>
#include <blaze/math/expressions/Transformation.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/submatrix/Dense.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATREPEATEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense matrix repeat() function.
// \ingroup dense_matrix_expression
//
// The DMatRepeatExpr class represents the compile time expression for repeating a dense matrix
// via the repeat() function.
*/
template< typename MT       // Type of the dense matrix
        , bool SO           // Storage order
        , size_t... CRAs >  // Compile time repeater arguments
class DMatRepeatExpr
   : public MatRepeatExpr< DenseMatrix< DMatRepeatExpr<MT,SO,CRAs...>, SO >, CRAs... >
   , private If_t< IsComputation_v<MT>, Computation, Transformation >
   , private RepeatExprData<2UL,CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RepeatExprData<2UL,CRAs...>;  //!< The type of the RepeatExprData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatRepeatExpr instance.
   using This = DMatRepeatExpr<MT,SO,CRAs...>;

   //! Base type of this DMatRepeatExpr instance.
   using BaseType = MatRepeatExpr< DenseMatrix<This,SO>, CRAs... >;

   using ResultType    = RepeatTrait_t<MT,CRAs...>;    //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations.
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.

   //! Composite data type of the dense matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

 public:
   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatRepeatExpr class.
   //
   // \param dm The dense matrix operand of the repeater expression.
   // \param args The number of repetitions.
   */
   template< typename... RRAs >  // Runtime repeater arguments
   explicit inline DMatRepeatExpr( const MT& dm, RRAs... args ) noexcept
      : DataType( args... )  // Base class initialization
      , dm_     ( dm )       // Dense matrix of the repeater expression
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
      BLAZE_INTERNAL_ASSERT( i < rows()   , "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
      return dm_( i%dm_.rows(), j%dm_.columns() );
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
      if( i >= rows() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= columns() ) {
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
      return dm_.rows() * this->template repetitions<0UL>();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.columns() * this->template repetitions<1UL>();
   }
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the dense matrix operand.
   //
   // \return The dense matrix operand.
   */
   inline Operand operand() const noexcept {
      return dm_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::repetitions;
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return IsExpression_v<MT> && dm_.canAlias( alias );
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
      return dm_.isAliased( alias );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return dm_.isAligned();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return false;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   Operand dm_;    //!< Dense matrix of the repeater expression.
   size_t  reps_;  //!< The number of repetitions.
   //**********************************************************************************************

   //**Assignment to row-major dense matrices******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix repeater expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side repeater expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix repeater
   // expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT2,SO2>& lhs, const DMatRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CompositeType_t<MT> A( serial( rhs.dm_ ) );  // Evaluation of the dense matrix operand

      const size_t reps0( rhs.template repetitions<0UL>() );
      const size_t reps1( rhs.template repetitions<1UL>() );
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
         for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
            submatrix( *lhs, rep0*M, rep1*N, M, N, unchecked ) = serial( A );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to row-major sparse matrices*****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix repeater expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side repeater expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix repeater
   // expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT2,SO2>& lhs, const DMatRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CompositeType_t<MT> A( serial( rhs.dm_ ) );  // Evaluation of the dense matrix operand

      const size_t reps0( rhs.template repetitions<0UL>() );
      const size_t reps1( rhs.template repetitions<1UL>() );
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      (*lhs).reset();
      (*lhs).reserve( reps0*reps1*M*N );

      if( SO2 == rowMajor )
      {
         for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
            for( size_t i=0UL; i<M; ++i ) {
               for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
                  for( size_t j=0UL; j<N; ++j ) {
                     (*lhs).append( rep0*M+i, rep1*N+j, A(i,j), true );
                  }
               }
               (*lhs).finalize( rep0*M+i );
            }
         }
      }
      else
      {
         for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
            for( size_t j=0UL; j<N; ++j ) {
               for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
                  for( size_t i=0UL; i<M; ++i ) {
                     (*lhs).append( rep0*M+i, rep1*N+j, A(i,j), true );
                  }
               }
               (*lhs).finalize( rep1*N+j );
            }
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to row-major dense matrices*********************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix repeater expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side repeater expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix
   // repeater expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CompositeType_t<MT> A( serial( rhs.dm_ ) );  // Evaluation of the dense matrix operand

      const size_t reps0( rhs.template repetitions<0UL>() );
      const size_t reps1( rhs.template repetitions<1UL>() );
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
         for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
            submatrix( *lhs, rep0*M, rep1*N, M, N, unchecked ) += serial( A );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to row-major dense matrices******************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix repeater expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side repeater expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix
   // repeater expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CompositeType_t<MT> A( serial( rhs.dm_ ) );  // Evaluation of the dense matrix operand

      const size_t reps0( rhs.template repetitions<0UL>() );
      const size_t reps1( rhs.template repetitions<1UL>() );
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
         for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
            submatrix( *lhs, rep0*M, rep1*N, M, N, unchecked ) -= serial( A );
         }
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to row-major dense matrices****************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix repeater expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side repeater expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix repeater expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatRepeatExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      CompositeType_t<MT> A( serial( rhs.dm_ ) );  // Evaluation of the dense matrix operand

      const size_t reps0( rhs.template repetitions<0UL>() );
      const size_t reps1( rhs.template repetitions<1UL>() );
      const size_t M( A.rows() );
      const size_t N( A.columns() );

      for( size_t rep0=0UL; rep0<reps0; ++rep0 ) {
         for( size_t rep1=0UL; rep1<reps1; ++rep1 ) {
            submatrix( *lhs, rep0*M, rep1*N, M, N, unchecked ) %= serial( A );
         }
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

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
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
/*!\brief Repeats the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be repeated.
// \param m The number of row-wise repetitions.
// \param n The number of column-wise repetitions.
// \return The repeated dense matrix.
//
// This function returns an expression representing the repeated dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<int,rowMajor> A1{ { 1, 0, -2 }, { 0, 5, 0 } };
   blaze::DynamicMatrix<int,columnMajor> B1{ { 0, -1 }, { 0, 4 }, { 7, 0 } };

   blaze::DynamicMatrix<int,rowMajor> A2;
   blaze::DynamicMatrix<int,columnMajor> B2;

   // ... Resizing and initialization

   // Repeating the 2x3 matrix 'A1' results in
   //
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //
   A2 = repeat( A1, 2UL, 3UL );

   // Repeating the 3x2 matrix 'B1' results in
   //
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //
   B2 = repeat( B1, 2UL, 3UL );
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) repeat( const DenseMatrix<MT,SO>& dm, size_t m, size_t n )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatRepeatExpr<MT,SO>;
   return ReturnType( *dm, m, n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Repeats the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be repeated.
// \return The repeated dense matrix.
//
// This function returns an expression representing the repeated dense matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   blaze::DynamicMatrix<int,rowMajor> A1{ { 1, 0, -2 }, { 0, 5, 0 } };
   blaze::DynamicMatrix<int,columnMajor> B1{ { 0, -1 }, { 0, 4 }, { 7, 0 } };

   blaze::DynamicMatrix<int,rowMajor> A2;
   blaze::DynamicMatrix<int,columnMajor> B2;

   // ... Resizing and initialization

   // Repeating the 2x3 matrix 'A1' results in
   //
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //    (  1  0 -2  1  0 -2  1  0 -2 )
   //    (  0  5  0  0  5  0  0  5  0 )
   //
   A2 = repeat<2UL,3UL>( A1 );

   // Repeating the 3x2 matrix 'B1' results in
   //
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //    (  0 -1  0 -1  0 -1 )
   //    (  0  4  0  4  0  4 )
   //    (  7  0  7  0  7  0 )
   //
   B2 = repeat<2UL,3UL>( B1 );
   \endcode
*/
template< size_t R0    // Compile time row-wise repetitions
        , size_t R1    // Compile time column-wise repetitions
        , typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) repeat( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatRepeatExpr<MT,SO,R0,R1>;
   return ReturnType( *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Repeating the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be repeated.
// \param m The number of row-wise repetitions.
// \param n The number of column-wise repetitions.
// \return The repeated dense matrix.
//
// This auxiliary overload of the \c repeat() function accepts both a compile time and a runtime
// expansion. The runtime argument is discarded in favor of the compile time argument.
*/
template< size_t R0    // Compile time row-wise repetitions
        , size_t R1    // Compile time column-wise repetitions
        , typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) repeat( const DenseMatrix<MT,SO>& dm, size_t m, size_t n )
{
   MAYBE_UNUSED( m, n );

   BLAZE_FUNCTION_TRACE;

   using ReturnType = const DMatRepeatExpr<MT,SO,R0,R1>;
   return ReturnType( *dm );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
