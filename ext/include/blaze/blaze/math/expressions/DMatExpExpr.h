//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatExpExpr.h
//  \brief Header file for the dense matrix exponential expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATEXPEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATEXPEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatExpExpr.h>
#include <blaze/math/shims/Exp.h>
#include <blaze/math/shims/Frexp.h>
#include <blaze/math/shims/Pow.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATEXPEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for dense matrix exponential.
// \ingroup dense_matrix_expression
//
// The DMatExpExpr class represents the compile time expression for the dense matrices exponential
// operation (see https://en.wikipedia.org/wiki/Matrix_exponential).
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
class DMatExpExpr
   : public MatExpExpr< DenseMatrix< DMatExpExpr<MT,SO>, SO > >
   , private Computation
{
 private:
   //**Type definitions****************************************************************************
   using RT = ResultType_t<MT>;   //!< Result type of the dense matrix expression.
   using ET = ElementType_t<MT>;  //!< Element type of the dense matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr size_t K = 24;  //!< The approximation limit for the exponential computation.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatExpExpr instance.
   using This = DMatExpExpr<MT,SO>;

   //! Base type of this DMatExpExpr instance.
   using BaseType = MatExpExpr< DenseMatrix<This,SO> >;

   using ResultType    = RemoveAdaptor_t<RT>;  //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<MT>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<MT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;    //!< Resulting element type.
   using ReturnType    = ReturnType_t<MT>;     //!< Return type for expression template evaluations.

   //! Data type for composite expression templates.
   using CompositeType = const ResultType;

   //! Composite data type of the dense matrix expression.
   using Operand = If_t< IsExpression_v<MT>, const MT, const MT& >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatExpExpr class.
   //
   // \param dm The dense matrix operand of the exponential expression.
   */
   explicit inline DMatExpExpr( const MT& dm ) noexcept
      : dm_( dm )  // Dense matrix of the exponential expression
   {
      BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return dm_.columns();
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return dm_.rows();
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
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return dm_.isAliased( alias );
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

 private:
   //**Member variables****************************************************************************
   Operand dm_;  //!< Dense matrix of the exponential expression.
   //**********************************************************************************************

   //**Assignment to dense matrices****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix exponential expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side exponential expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix exponential
   // expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void assign( DenseMatrix<MT2,SO2>& lhs, const DMatExpExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const size_t N( rhs.rows() );

      if( IsDiagonal_v<MT> || N < 2UL )
      {
         assign( *lhs, rhs.dm_ );
         for( size_t i=0UL; i<N; ++i ) {
            (*lhs)(i,i) = exp( (*lhs)(i,i) );
         }
      }
      else
      {
         using BT = UnderlyingBuiltin_t<ET>;

         const BT norm    ( maxNorm( rhs.dm_ ) );
         const BT log2norm( norm > BT(0) ? log2( norm ) : BT(0) );

         int exponent( 0 );
         frexp( ceil( log2norm ), &exponent );
         exponent = max( 0, exponent+1 );

         ResultType R( rhs.dm_ / pow( 2.0, double(exponent) ) );
         ResultType A( R );
         ResultType B( R );

         for( size_t i=0UL; i<N; ++i ) {
            B(i,i) += BT(1);
         }

         BT factor( 1 );
         for( size_t k=2UL; k<K; ++k ) {
            factor *= BT( k );
            A *= R;
            addAssign( B, ( A / factor ) );
         }

         for( int i=0; i<exponent; ++i ) {
            B *= B;
         }

         assign( *lhs, B );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Assignment to sparse matrices***************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Assignment of a dense matrix exponential expression to a sparse matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side sparse matrix.
   // \param rhs The right-hand side exponential expression to be assigned.
   // \return void
   //
   // This function implements the performance optimized assignment of a dense matrix exponential
   // expression to a sparse matrix.
   */
   template< typename MT2  // Type of the target sparse matrix
           , bool SO2 >    // Storage order of the target sparse matrix
   friend inline void assign( SparseMatrix<MT2,SO2>& lhs, const DMatExpExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ResultType );
      BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( ResultType, SO );
      BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const ResultType tmp( serial( rhs ) );
      assign( *lhs, tmp );
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to dense matrices*******************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Addition assignment of a dense matrix exponential expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side exponential expression to be added.
   // \return void
   //
   // This function implements the performance optimized addition assignment of a dense matrix
   // exponential expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void addAssign( DenseMatrix<MT2,SO2>& lhs, const DMatExpExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const size_t N( rhs.rows() );

      if( IsDiagonal_v<MT> || N < 2UL )
      {
         CompositeType_t<MT> tmp( rhs.dm_ );
         for( size_t i=0UL; i<N; ++i ) {
            (*lhs)(i,i) += exp( tmp(i,i) );
         }
      }
      else
      {
         using BT = UnderlyingBuiltin_t<ET>;

         const BT norm    ( maxNorm( rhs.dm_ ) );
         const BT log2norm( norm > BT(0) ? log2( norm ) : BT(0) );

         int exponent( 0 );
         frexp( ceil( log2norm ), &exponent );
         exponent = max( 0, exponent+1 );

         ResultType R( rhs.dm_ / pow( 2.0, double(exponent) ) );
         ResultType A( R );
         ResultType B( R );

         for( size_t i=0UL; i<N; ++i ) {
            B(i,i) += BT(1);
         }

         BT factor( 1 );
         for( size_t k=2UL; k<K; ++k ) {
            factor *= BT( k );
            A *= R;
            addAssign( B, ( A / factor ) );
         }

         for( int i=0; i<exponent; ++i ) {
            B *= B;
         }

         addAssign( *lhs, B );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Addition assignment to sparse matrices******************************************************
   // No special implementation for the addition assignment to sparse matrices.
   //**********************************************************************************************

   //**Subtraction assignment to dense matrices****************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Subtraction assignment of a dense matrix exponential expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side exponential expression to be subtracted.
   // \return void
   //
   // This function implements the performance optimized subtraction assignment of a dense matrix
   // exponential expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void subAssign( DenseMatrix<MT2,SO2>& lhs, const DMatExpExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const size_t N( rhs.rows() );

      if( IsDiagonal_v<MT> || N < 2UL )
      {
         CompositeType_t<MT> tmp( rhs.dm_ );
         for( size_t i=0UL; i<N; ++i ) {
            (*lhs)(i,i) -= exp( tmp(i,i) );
         }
      }
      else
      {
         using BT = UnderlyingBuiltin_t<ET>;

         const BT norm    ( maxNorm( rhs.dm_ ) );
         const BT log2norm( norm > BT(0) ? log2( norm ) : BT(0) );

         int exponent( 0 );
         frexp( ceil( log2norm ), &exponent );
         exponent = max( 0, exponent+1 );

         ResultType R( rhs.dm_ / pow( 2.0, double(exponent) ) );
         ResultType A( R );
         ResultType B( R );

         for( size_t i=0UL; i<N; ++i ) {
            B(i,i) += BT(1);
         }

         BT factor( 1 );
         for( size_t k=2UL; k<K; ++k ) {
            factor *= BT( k );
            A *= R;
            addAssign( B, ( A / factor ) );
         }

         for( int i=0; i<exponent; ++i ) {
            B *= B;
         }

         subAssign( *lhs, B );
      }
   }
   /*! \endcond */
   //**********************************************************************************************

   //**Subtraction assignment to sparse matrices***************************************************
   // No special implementation for the subtraction assignment to sparse matrices.
   //**********************************************************************************************

   //**Schur product assignment to dense matrices**************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Schur product assignment of a dense matrix exponential expression to a dense matrix.
   // \ingroup dense_matrix
   //
   // \param lhs The target left-hand side dense matrix.
   // \param rhs The right-hand side exponential expression for the Schur product.
   // \return void
   //
   // This function implements the performance optimized Schur product assignment of a dense
   // matrix exponential expression to a dense matrix.
   */
   template< typename MT2  // Type of the target dense matrix
           , bool SO2 >    // Storage order of the target dense matrix
   friend inline void schurAssign( DenseMatrix<MT2,SO2>& lhs, const DMatExpExpr& rhs )
   {
      BLAZE_FUNCTION_TRACE;

      BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == rhs.rows()   , "Invalid number of rows"    );
      BLAZE_INTERNAL_ASSERT( (*lhs).columns() == rhs.columns(), "Invalid number of columns" );

      const size_t N( rhs.rows() );

      if( IsDiagonal_v<MT> || N < 2UL )
      {
         CompositeType_t<MT> tmp( rhs.dm_ );
         for( size_t i=0UL; i<N; ++i ) {
            (*lhs)(i,i) *= exp( tmp(i,i) );
         }
      }
      else
      {
         using BT = UnderlyingBuiltin_t<ET>;

         const BT norm    ( maxNorm( rhs.dm_ ) );
         const BT log2norm( norm > BT(0) ? log2( norm ) : BT(0) );

         int exponent( 0 );
         frexp( ceil( log2norm ), &exponent );
         exponent = max( 0, exponent+1 );

         ResultType R( rhs.dm_ / pow( 2.0, double(exponent) ) );
         ResultType A( R );
         ResultType B( R );

         for( size_t i=0UL; i<N; ++i ) {
            B(i,i) += BT(1);
         }

         BT factor( 1 );
         for( size_t k=2UL; k<K; ++k ) {
            factor *= BT( k );
            A *= R;
            addAssign( B, ( A / factor ) );
         }

         for( int i=0; i<exponent; ++i ) {
            B *= B;
         }

         schurAssign( *lhs, B );
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
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Calculation of the exponential of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix for the matrix exponential.
// \return The exponential of the given matrix.
//
// This function returns an expression representing the exponential of the given dense matrix:

                  \f[ e^X = \sum\limits_{k=0}^\infty \frac{1}{k!} X^k \f]

// Example:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A, B;
   // ... Resizing and initialization
   B = matexp( A );
   \endcode

// \note The matrix exponential can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note It is not possible to use any kind of view on the expression object returned by the
// \c matexp() function. Also, it is not possible to access individual elements via the function
// call operator on the expression object:

   \code
   row( matexp( A ), 2UL );  // Compilation error: Views cannot be used on an matexp() expression!
   matexp( A )(1,2);         // Compilation error: It is not possible to access individual elements!
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) matexp( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   if( !isSquare( *dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   using ReturnType = const DMatExpExpr<MT,SO>;
   return ReturnType( *dm );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given dense matrix exponential.
// \ingroup dense_matrix
//
// \param dm The given dense matrix exponential.
// \return The determinant of the given matrix exponential.
//
// This function computes the determinant of the given dense matrix exponential.
//
// \note The computation of the determinant is numerically unreliable since especially for large
// matrices the value can overflow during the computation. Please note that this function does
// not guarantee that it is possible to compute the determinant with the given matrix!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline decltype(auto) det( const DMatExpExpr<MT,SO>& dm )
{
   BLAZE_FUNCTION_TRACE;

   return det( evaluate( *dm ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
