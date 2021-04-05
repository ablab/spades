//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatGenExpr.h
//  \brief Header file for the dense matrix generator expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATGENEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATGENEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Expression.h>
#include <blaze/math/constraints/StorageOrder.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatGenExpr.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/views/Check.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/util/Assert.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/TypeList.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DMATGENEXPR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Expression object for the dense matrix generate() function.
// \ingroup dense_matrix_expression
//
// The DMatGenExpr class represents the compile time expression for the generation of a custom
// dense matrix via the generate() function.
*/
template< typename MT  // Type of the dense matrix
        , typename OP  // Type of the custom binary operation
        , bool SO >    // Storage order
class DMatGenExpr
   : public MatGenExpr< DenseMatrix< DMatGenExpr<MT,OP,SO>, SO > >
   , private Computation
{
 public:
   //**Type definitions****************************************************************************
   //! Type of this DMatGenExpr instance.
   using This = DMatGenExpr<MT,OP,SO>;

   //! Base type of this DMatGenExpr instance.
   using BaseType = MatGenExpr< DenseMatrix<This,SO> >;

   using ResultType    = RemoveCVRef_t<MT>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<MT>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;    //!< Resulting element type.

   //! Return type for expression template evaluations.
   using ReturnType = decltype( std::declval<OP>()( std::declval<size_t>(), std::declval<size_t>() ) );

   //! Data type for composite expression templates.
   using CompositeType = const DMatGenExpr&;

   //! Data type of the custom binary operation.
   using Operation = OP;
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the dense matrix generator expression.
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
      // \param m Row index to the initial matrix element.
      // \param n Column index to the initial matrix element.
      // \param op The custom binary operation.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator( size_t m, size_t n, OP op )
         : m_ ( m )              // Row index of the current matrix element
         , n_ ( n )              // Column index of the current matrix element
         , op_( std::move(op) )  // The custom binary operation
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator+=( size_t inc ) {
         if( SO )
            m_ += inc;
         else
            n_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator-=( size_t dec ) {
         if( SO )
            m_ += dec;
         else
            n_ += dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator++() {
         if( SO )
            ++m_;
         else
            ++n_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator++( int ) {
         if( SO )
            return ConstIterator( m_++, n_, op_ );
         else
            return ConstIterator( m_, n_++, op_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE ConstIterator& operator--() {
         if( SO )
            --m_;
         else
            --n_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const ConstIterator operator--( int ) {
         if( SO )
            return ConstIterator( m_--, n_, op_ );
         else
            return ConstIterator( m_, n_--, op_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline BLAZE_DEVICE_CALLABLE ReturnType operator*() const {
         return op_( m_, n_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ == rhs.m_;
         else
            return n_ == rhs.n_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ != rhs.m_;
         else
            return n_ != rhs.n_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ < rhs.m_;
         else
            return n_ < rhs.n_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ > rhs.m_;
         else
            return n_ > rhs.n_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ <= rhs.m_;
         else
            return n_ <= rhs.n_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ >= rhs.m_;
         else
            return n_ >= rhs.n_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const ConstIterator& rhs ) const noexcept {
         if( SO )
            return m_ - rhs.m_;
         else
            return n_ - rhs.n_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ConstIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator+( const ConstIterator& it, size_t inc ) {
         if( SO )
            return ConstIterator( it.m_ + inc, it.n_, it.op_ );
         else
            return ConstIterator( it.m_, it.n_ + inc, it.op_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ConstIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator+( size_t inc, const ConstIterator& it ) {
         if( SO )
            return ConstIterator( it.m_ + inc, it.n_, it.op_ );
         else
            return ConstIterator( it.m_, it.n_ + inc, it.op_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ConstIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const ConstIterator operator-( const ConstIterator& it, size_t dec ) {
         if( SO )
            return ConstIterator( it.m_ - dec, it.n_, it.op_ );
         else
            return ConstIterator( it.m_, it.n_ - dec, it.op_ );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      size_t m_;   //!< Row index of the current matrix element.
      size_t n_;   //!< Column index of the current matrix element.
      OP     op_;  //!< The custom binary operation.
      //*******************************************************************************************
   };
   //**********************************************************************************************

 public:
   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = true;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for the DMatGenExpr class.
   //
   // \param m The number of rows of the dense matrix generator.
   // \param n The number of columns of the dense matrix generator.
   // \param op The custom binary operation.
   */
   inline DMatGenExpr( size_t m, size_t n, OP&& op ) noexcept
      : m_ ( m )                // The number of rows of the dense matrix generator
      , n_ ( n )                // The number of columns of the dense matrix generator
      , op_( std::move( op ) )  // The custom binary operation
   {}
   //**********************************************************************************************

   //**Access operator*****************************************************************************
   /*!\brief 2D-access to the matrix elements.
   //
   // \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
   // \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator()( size_t i, size_t j ) const noexcept {
      BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index"    );
      BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
      return op_(i,j);
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
      if( i >= m_ ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
      }
      if( j >= n_ ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
      }
      return (*this)(i,j);
   }
   //**********************************************************************************************

   //**Begin function******************************************************************************
   /*!\brief Returns an iterator to the first non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator to the first non-zero element of row/column \a i.
   */
   inline ConstIterator begin( size_t i ) const {
      return ConstIterator( ( SO ? 0UL : i ), ( SO ? i : 0UL ), op_ );
   }
   //**********************************************************************************************

   //**End function********************************************************************************
   /*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
   //
   // \param i The row/column index.
   // \return Iterator just past the last non-zero element of row/column \a i.
   */
   inline ConstIterator end( size_t i ) const {
      return ConstIterator( ( SO ? m_ : i ), ( SO ? i : n_ ), op_ );
   }
   //**********************************************************************************************

   //**Rows function*******************************************************************************
   /*!\brief Returns the current number of rows of the matrix.
   //
   // \return The number of rows of the matrix.
   */
   inline size_t rows() const noexcept {
      return m_;
   }
   //**********************************************************************************************

   //**Columns function****************************************************************************
   /*!\brief Returns the current number of columns of the matrix.
   //
   // \return The number of columns of the matrix.
   */
   inline size_t columns() const noexcept {
      return n_;
   }
   //**********************************************************************************************

   //**Operation access****************************************************************************
   /*!\brief Returns a copy of the custom binary operation.
   //
   // \return A copy of the custom binary operation.
   */
   inline Operation operation() const {
      return op_;
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
      MAYBE_UNUSED( alias );
      return false;
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
      MAYBE_UNUSED( alias );
      return false;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operands of the expression are properly aligned in memory.
   //
   // \return \a true in case the operands are aligned, \a false if not.
   */
   inline bool isAligned() const noexcept {
      return true;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can be used in SMP assignments.
   //
   // \return \a true in case the expression can be used in SMP assignments, \a false if not.
   */
   inline bool canSMPAssign() const noexcept {
      return true;
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   size_t m_;      //!< The number of rows of the dense matrix generator.
   size_t n_;      //!< The number of columns of the dense matrix generator.
   Operation op_;  //!< The custom binary operation.
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_WITH_STORAGE_ORDER( MT, SO );
   BLAZE_CONSTRAINT_MUST_NOT_BE_EXPRESSION_TYPE( MT );
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
/*!\brief Generates a new dense matrix filled via the given custom binary operation.
// \ingroup dense_matrix
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param op The custom binary operation.
// \return The newly generated dense matrix.
// \exception std::invalid_argument Invalid size specification.
//
// The \a generate() function returns a dense matrix filled elementwise via the given custom
// binary operation. The function returns an expression representing an accordingly filled matrix
// of the specified type \a MT.\n
// The following example demonstrates the use of the \a generate() function:

   \code
   using blaze::generate;
   using blaze::rowMajor;
   using blaze::columnMajor>

   // Generates the uniform integer matrix ( ( 2, 2, 2 ), ( 2, 2, 2 ) )
   using A = blaze::DynamicMatrix<int,rowMajor>;
   const auto a = generate<A>( 2UL, 3UL, []( size_t i, size_t j ){ return 2; } );

   // Generates the linearly spaced float matrix ( ( 2.1, 3.2, 4.3 ), ( 5.4, 6.5, 7.6 ) )
   using B = blaze::DynamicMatrix<float,rowMajor>;
   const auto b = generate<B>( 2UL, 3UL, []( size_t i, size_t j ){ return 2.1F + 1.1F*(i*3UL+j); } );

   // Generates the logarithmically spaced double vector ( ( 1.0, 10.0 ), ( 100.0, 1000.0 ) )
   using C = blaze::DynamicMatrix<double,rowMajor>;
   const auto c = generate<C>( 2UL, 2UL, []( size_t i, size_t j ) { return blaze::exp10( 1.0 + 1.0*(i*2UL+j) ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = StaticVector<int,2UL>;
   using D = blaze::DynamicMatrix<VT,columnMajor>;
   const auto d = generate<D>( 2UL, 2UL, []( size_t i, size_t j ) { return evaluate( VT{ 1, 2 } + (i*2UL+j) ); } );
   \endcode

// In case the specified size does not match the size of the given matrix type \a MT, a
// \a std::invalid_argument is thrown.
*/
template< typename MT    // Type of the dense matrix
        , typename OP >  // Type of the custom binary operation
inline decltype(auto) generate( size_t m, size_t n, OP op )
{
   BLAZE_FUNCTION_TRACE;

   if( Size_v<MT,0UL> != DefaultSize_v && size_t( Size_v<MT,0UL> ) != m &&
       Size_v<MT,1UL> != DefaultSize_v && size_t( Size_v<MT,1UL> ) != n ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid size specification" );
   }

   using ReturnType = const DMatGenExpr< MT, OP, StorageOrder_v<MT> >;
   return ReturnType( m, n, std::move( op ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Generates a new dense matrix filled via the given custom binary operation.
// \ingroup dense_matrix
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param op The custom binary operation.
// \return The newly generated dense matrix.
//
// The \c generate() function returns a dense matrix filled elementwise via the given custom
// binary operation. By default, the returned matrix is a row-major matrix, but this setting
// can be changed via the \c BLAZE_DEFAULT_STORAGE_ORDER switch. Alternatively it is possible
// to specify the storage order explicitly.\n
// The following example demonstrates the use of the \a generate() function:

   \code
   using blaze::generate;
   using blaze::rowMajor;
   using blaze::columnMajor>

   // Generates the uniform integer matrix ( ( 2, 2, 2 ), ( 2, 2, 2 ) )
   blaze::DynamicMatrix<int,rowMajor> A;
   A = generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2; } );

   // Generates the linearly spaced float matrix ( ( 2.1, 3.2, 4.3 ), ( 5.4, 6.5, 7.6 ) )
   blaze::DynamicMatrix<float,rowMajor> B;
   B = generate( 2UL, 3UL, []( size_t i, size_t j ){ return 2.1F + 1.1F*(i*3UL+j); } );

   // Generates the logarithmically spaced double vector ( ( 1.0, 10.0 ), ( 100.0, 1000.0 ) )
   blaze::DynamicMatrix<double,rowMajor> C;
   C = generate<rowMajor>( 2UL, 2UL, []( size_t i, size_t j ) { return blaze::exp10( 1.0 + 1.0*(i*2UL+j) ); } );

   // Generates the vector of integer vectors ( ( 1, 2 ), ( 2, 3 ), ( 3, 4 ), ( 4, 5 ) )
   using VT = StaticVector<int,2UL>;
   blaze::DynamicMatrix<VT,columnMajor> D;
   D = generate<columnMajor>( 2UL, 2UL, []( size_t i, size_t j ) { return evaluate( VT{ 1, 2 } + (i*2UL+j) ); } );
   \endcode
*/
template< bool SO = defaultStorageOrder  // Storage order
        , typename OP >                  // Type of the custom binary operation
inline decltype(auto) generate( size_t m, size_t n, OP op )
{
   BLAZE_FUNCTION_TRACE;

   using ET = RemoveCVRef_t< decltype( std::declval<OP>()( std::declval<size_t>(), std::declval<size_t>() ) ) >;
   using MT = DynamicMatrix<ET,SO>;

   return generate<MT>( m, n, std::move( op ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBMATRIX)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submatrix of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the submatrix.
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the dense matrix generator expression.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given dense
// matrix generator expression.
*/
template< AlignmentFlag AF    // Alignment flag
        , size_t I            // Index of the first row
        , size_t J            // Index of the first column
        , size_t M            // Number of rows
        , size_t N            // Number of columns
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RSAs >  // Optional arguments
inline decltype(auto) submatrix( const DMatGenExpr<MT,OP,SO>& expr, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RSAs...>, Unchecked > );

   if( isChecked ) {
      if( ( I + M > expr.rows() ) || ( J + N > expr.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I + M <= expr.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( J + N <= expr.columns(), "Invalid submatrix specification" );
   }

   return generate( M, N, [op=expr.operation()]( size_t i, size_t j ) {
      return op( i+I, j+J );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific submarix of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the submatrix.
// \param row The index of the first row of the submatrix.
// \param column The index of the first column of the submatrix.
// \param m The number of rows of the submatrix.
// \param n The number of columns of the submatrix.
// \param args The optional submatrix arguments.
// \return View on the specified submatrix of the dense matrix generator expression.
// \exception std::invalid_argument Invalid submatrix specification.
//
// This function returns an expression representing the specified submatrix of the given dense
// matrix generator expression.
*/
template< AlignmentFlag AF    // Alignment flag
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RSAs >  // Optional arguments
inline decltype(auto)
   submatrix( const DMatGenExpr<MT,OP,SO>& expr, size_t row, size_t column, size_t m, size_t n, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RSAs...>, Unchecked > );

   if( isChecked ) {
      if( ( row + m > expr.rows() ) || ( column + n > expr.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row    + m <= expr.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column + n <= expr.columns(), "Invalid submatrix specification" );
   }

   return generate( m, n, [op=expr.operation(),row,column]( size_t i, size_t j ) {
      return op( i+row, j+column );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ROW)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the row.
// \param args The optional row arguments.
// \return View on the specified row of the dense matrix generator expression.
//
// This function returns an expression representing the specified row of the given dense matrix
// generator expression.
*/
template< size_t I            // Row index
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const DMatGenExpr<MT,OP,SO>& expr, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      if( expr.rows() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I < expr.rows(), "Invalid row access index" );
   }

   return generate<rowVector>( expr.columns(), [op=expr.operation()]( size_t i ){
      return op( I, i );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the row.
// \param index The index of the row.
// \param args The optional row arguments.
// \return View on the specified row of the dense matrix generator expression.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given dense matrix
// generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const DMatGenExpr<MT,OP,SO>& expr, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      if( expr.rows() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( index < expr.rows(), "Invalid row access index" );
   }

   return generate<rowVector>( expr.columns(), [op=expr.operation(),index]( size_t i ){
      return op( index, i );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ROWS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param args The optional row arguments.
// \return View on the specified rows of the dense matrix generator expression.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified rows of the given dense matrix
// generator expression.
*/
template< size_t I            // First required row index
        , size_t... Is        // Remaining required row indices
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const DMatGenExpr<MT,OP,SO>& expr, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( expr.rows() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid elements specification" );
         }
      }
   }

   return generate( sizeof...(Is)+1UL, expr.columns(), [op=expr.operation()]( size_t i, size_t j ) {
      static constexpr size_t indices[] = { I, Is... };
      return op( indices[i], j );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args The optional row arguments.
// \return View on the specified rows of the dense matrix generator expression.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given dense matrix
// generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename T          // Type of the row indices
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const DMatGenExpr<MT,OP,SO> expr, T* indices, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.rows() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices( indices, indices+n );

   return generate( n, expr.columns(), [op=expr.operation(),newIndices]( size_t i, size_t j ){
      return op( newIndices[i], j );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific rows of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The optional row arguments.
// \return View on the specified rows of the dense matrix generator expression.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given dense matrix
// generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename P          // Type of the index producer
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) rows( const DMatGenExpr<MT,OP,SO>& expr, P p, size_t n, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RRAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.rows() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
         }
      }
   }

   return generate( n, expr.columns(), [op=expr.operation(),p]( size_t i, size_t j ){
      return op( p(i), j );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMN)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the dense matrix generator expression.
//
// This function returns an expression representing the specified column of the given dense matrix
// generator expression.
*/
template< size_t I            // Column index
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... CRAs >  // Optional column arguments
inline decltype(auto) column( const DMatGenExpr<MT,OP,SO>& expr, CRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<CRAs...>, Unchecked > );

   if( isChecked ) {
      if( expr.columns() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( I < expr.columns(), "Invalid column access index" );
   }

   return generate<columnVector>( expr.rows(), [op=expr.operation()]( size_t i ){
      return op( i, I );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the column.
// \param index The index of the column.
// \param args The optional column arguments.
// \return View on the specified column of the dense matrix generator expression.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given dense matrix
// generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... CRAs >  // Runtime column arguments
inline decltype(auto) column( const DMatGenExpr<MT,OP,SO>& expr, size_t index, CRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<CRAs...>, Unchecked > );

   if( isChecked ) {
      if( expr.columns() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( index < expr.columns(), "Invalid column access index" );
   }

   return generate<columnVector>( expr.rows(), [op=expr.operation(),index]( size_t i ){
      return op( i, index );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMNS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param args The optional column arguments.
// \return View on the specified columns of the dense matrix generator expression.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified columns of the given dense
// matrix generator expression.
*/
template< size_t I            // First required column index
        , size_t... Is        // Remaining required column indices
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const DMatGenExpr<MT,OP,SO>& expr, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      static constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( expr.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid elements specification" );
         }
      }
   }

   return generate( expr.rows(), sizeof...(Is)+1UL, [op=expr.operation()]( size_t i, size_t j ) {
      static constexpr size_t indices[] = { I, Is... };
      return op( i, indices[j] );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the dense matrix generator expression.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given dense
// matrix generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const DMatGenExpr<MT,OP,SO> expr, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.columns() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices( indices, indices+n );

   return generate( expr.rows(), n, [op=expr.operation(),newIndices]( size_t i, size_t j ){
      return op( i, newIndices[j] );
   } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The given dense matrix generator expression.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the dense matrix generator expression.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given dense
// matrix generator expression.
*/
template< typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const DMatGenExpr<MT,OP,SO>& expr, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( args... );

   constexpr bool isChecked( !Contains_v< TypeList<RCAs...>, Unchecked > );

   if( isChecked ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( expr.columns() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   return generate( expr.rows(), n, [op=expr.operation(),p]( size_t i, size_t j ){
      return op( i, p(j) );
   } );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (BAND)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given dense matrix generator expression.
// \ingroup dense_matrix
//
// \param expr The dense matrix generator expression containing the band.
// \param args The runtime band arguments.
// \return View on the specified band of the dense matrix generator expression.
//
// This function returns an expression representing the specified band of the given dense matrix
// generator expression.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the dense matrix
        , typename OP         // Type of the custom binary operation
        , bool SO             // Storage order
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const DMatGenExpr<MT,OP,SO>& expr, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const BandData<CBAs...> bd( args... );

   const size_t row   ( bd.band() >= 0L ? 0UL : -bd.band() );
   const size_t column( bd.band() >= 0L ? bd.band() :  0UL );
   const size_t size  ( min( expr.rows() - row, expr.columns() - column ) );

   const bool isChecked( !Contains_v< TypeList<RBAs...>, Unchecked > );

   if( isChecked ) {
      if( ( bd.band() > 0L && column >= expr.columns() ) ||
          ( bd.band() < 0L && row >= expr.rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( bd.band() <= 0L || column < expr.columns(), "Invalid band access index" );
      BLAZE_USER_ASSERT( bd.band() >= 0L || row < expr.rows(), "Invalid band access index" );
   }

   return generate( size, [op=expr.operation(),row,column]( size_t i ){
      return op( i+row, i+column );
   } );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
