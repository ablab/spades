//=================================================================================================
/*!
//  \file blaze/math/views/column/Sparse.h
//  \brief Column specialization for sparse matrices
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_SPARSE_H_
#define _BLAZE_MATH_VIEWS_COLUMN_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/column/ColumnData.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Column for columns on column-major sparse matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of column-major
// sparse matrices.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,true,false,SF,CCAs...>
   : public View< SparseVector< Column<MT,true,false,SF,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,true,false,SF,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< SparseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator= ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator+=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator-=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>&       rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>&       rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   reserve( size_t n );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const ElementType& value );
   inline Iterator insert( size_t index, const ElementType& value );
   inline void     append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t index );
   inline Iterator erase( Iterator pos );
   inline Iterator erase( Iterator first, Iterator last );

   template< typename Pred, typename = DisableIf_t< IsIntegral_v<Pred> > >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for columns on column-major sparse matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,true,false,SF,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Reference
   Column<MT,true,false,SF,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstReference
   Column<MT,true,false,SF,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Reference
   Column<MT,true,false,SF,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstReference
   Column<MT,true,false,SF,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::begin()
{
   return matrix_.begin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::begin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::cbegin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::end()
{
   return matrix_.end( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::end() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::cend() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   const InitializerVector<ElementType,false> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Column.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator=( const Column& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && column() == rhs.column() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( rhs.nonZeros() );
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for sparse vectors.
//
// \param rhs Sparse vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( right.nonZeros() );
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator+=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a sparse vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator+=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a dense vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator-=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a sparse vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator-=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the sparse column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,true,false,SF,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,true,false,SF,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the sparse column.
//
// \return The size of the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,false,SF,CCAs...>::size() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,false,SF,CCAs...>::capacity() const noexcept
{
   return matrix_.capacity( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,false,SF,CCAs...>::nonZeros() const
{
   return matrix_.nonZeros( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,true,false,SF,CCAs...>::reset()
{
   matrix_.reset( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
void Column<MT,true,false,SF,CCAs...>::reserve( size_t n )
{
   matrix_.reserve( column(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new sparse column capacity.
//
// \return The new sparse column capacity.
//
// This function calculates a new column capacity based on the current capacity of the sparse
// column. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,true,false,SF,CCAs...>::extendCapacity() const noexcept
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = max( nonzeros, 7UL    );
   nonzeros = min( nonzeros, size() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::set( size_t index, const ElementType& value )
{
   return matrix_.set( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,true,false,SF,CCAs...>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( index, column(), value, check );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,true,false,SF,CCAs...>::erase( size_t index )
{
   matrix_.erase( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::erase( Iterator pos )
{
   return matrix_.erase( column(), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( column(), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Column<MT,true,false,SF,CCAs...>::erase( Pred predicate )
{
   matrix_.erase( column(), begin(), end(), predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse column. The
// elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements and to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( col2.begin(), col2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred >   // Type of the unary predicate
inline void Column<MT,true,false,SF,CCAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   matrix_.erase( column(), first, last, predicate );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::find( size_t index )
{
   return matrix_.find( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::find( size_t index ) const
{
   return matrix_.find( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::lowerBound( size_t index )
{
   return matrix_.lowerBound( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::Iterator
   Column<MT,true,false,SF,CCAs...>::upperBound( size_t index )
{
   return matrix_.upperBound( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,true,false,SF,CCAs...>::ConstIterator
   Column<MT,true,false,SF,CCAs...>::upperBound( size_t index ) const
{
   return matrix_.upperBound( index, column() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,true,false,SF,CCAs...>&
   Column<MT,true,false,SF,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,true,false,SF,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,true,false,SF,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,true,false,SF,CCAs...>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( column() ) == matrix_.capacity( column() ) )
         matrix_.reserve( column(), extendCapacity() );

      matrix_.append( i, column(), (*rhs)[i], true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,false,SF,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      matrix_.append( element->index(), column(), element->value(), true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,true,false,SF,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( column() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,false,SF,CCAs...>::addAssign( const SparseVector<VT,false>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( column() );
   matrix_.reserve( column(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,true,false,SF,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( column() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,true,false,SF,CCAs...>::subAssign( const SparseVector<VT,false>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( column() );
   matrix_.reserve( column(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL ROW-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Column for general row-major sparse matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of general
// row-major sparse matrices.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,false,false,false,CCAs...>
   : public View< SparseVector< Column<MT,false,false,false,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,false,false,false,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< SparseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;
   //**********************************************************************************************

   //**ColumnElement class definition**************************************************************
   /*!\brief Access proxy for a specific element of the sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnElement class.
      //
      // \param pos Iterator to the current position within the sparse column.
      // \param row The row index.
      */
      inline ColumnElement( IteratorType pos, size_t row )
         : pos_( pos )  // Iterator to the current position within the sparse column
         , row_( row )  // Index of the according row
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse column element.
      //
      // \param value The new value of the sparse column element.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse column element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const ColumnElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse column element.
      //
      // \return The current value of the sparse column element.
      */
      inline decltype(auto) value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return row_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse column.
      size_t       row_;  //!< Index of the according row.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**ColumnIterator class definition*************************************************************
   /*!\brief Iterator over the elements of the sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;               //!< The iterator category.
      using ValueType        = ColumnElement<MatrixType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                               //!< Pointer return type.
      using ReferenceType    = ValueType;                               //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                               //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the ColumnIterator class.
      */
      inline ColumnIterator()
         : matrix_( nullptr )  // The sparse matrix containing the column
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the column
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {
         for( ; row_<matrix_->rows(); ++row_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( row_ ) ) break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator.
      */
      inline ColumnIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the column
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ( pos     )  // Iterator to the current sparse element
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline ColumnIterator( const ColumnIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the column
         , row_   ( it.row_    )  // The current row index
         , column_( it.column_ )  // The current column index
         , pos_   ( it.pos_    )  // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ColumnIterator& operator++() {
         ++row_;
         for( ; row_<matrix_->rows(); ++row_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( row_ ) ) break;
         }

         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ColumnIterator operator++( int ) {
         const ColumnIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, row_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, row_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
         return row_ == rhs.row_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnIterator<MatrixType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two column iterators.
      //
      // \param rhs The right-hand side column iterator.
      // \return The number of elements between the two column iterators.
      */
      inline DifferenceType operator-( const ColumnIterator& rhs ) const {
         size_t counter( 0UL );
         for( size_t i=rhs.row_; i<row_; ++i ) {
            if( matrix_->find( i, column_ ) != matrix_->end( i ) )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The sparse matrix containing the column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnIterator;
      template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CCAs2 > friend class Column;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = ColumnIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, ColumnIterator< MT, Iterator_t<MT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>& rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const;
   inline size_t capacity() const;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   reserve( size_t n );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const ElementType& value );
   inline Iterator insert( size_t index, const ElementType& value );
   inline void     append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t index );
   inline Iterator erase( Iterator pos );
   inline Iterator erase( Iterator first, Iterator last );

   template< typename Pred, typename = DisableIf_t< IsIntegral_v<Pred> > >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const;
   template< typename Other > inline bool isAliased( const Other* alias ) const;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const Vector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const Vector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for columns on row-major sparse matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,false,false,false,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Reference
   Column<MT,false,false,false,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstReference
   Column<MT,false,false,false,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(index,column());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Reference
   Column<MT,false,false,false,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstReference
   Column<MT,false,false,false,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::begin()
{
   return Iterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::begin() const
{
   return ConstIterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::cbegin() const
{
   return ConstIterator( matrix_, 0UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::end()
{
   return Iterator( matrix_, size(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::end() const
{
   return ConstIterator( matrix_, size(), column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::cend() const
{
   return ConstIterator( matrix_, size(), column() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   const InitializerVector<ElementType,false> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Column.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator=( const Column& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && column() == rhs.column() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      assign( left, tmp );
   }
   else {
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const CompositeType_t<VT> tmp( *rhs );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator+=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator-=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the sparse column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,false,false,false,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,false,false,false,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the column.
//
// \return The size of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,false,CCAs...>::size() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,false,CCAs...>::capacity() const
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,false,CCAs...>::nonZeros() const
{
   size_t counter( 0UL );
   for( ConstIterator element=begin(); element!=end(); ++element ) {
      ++counter;
   }
   return counter;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,false,CCAs...>::reset()
{
   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( column()+1UL )
                           :( column() ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( column() )
                           :( column()+1UL ) )
                        :( size() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_.erase( i, column() );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
void Column<MT,false,false,false,CCAs...>::reserve( size_t n )
{
   MAYBE_UNUSED( n );

   return;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::set( size_t index, const ElementType& value )
{
   return Iterator( matrix_, index, column(), matrix_.set( index, column(), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::insert( size_t index, const ElementType& value )
{
   return Iterator( matrix_, index, column(), matrix_.insert( index, column(), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,false,CCAs...>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( index, column(), value );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,false,CCAs...>::erase( size_t index )
{
   matrix_.erase( index, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::erase( Iterator pos )
{
   const size_t row( pos.row_ );

   if( row == size() )
      return pos;

   matrix_.erase( row, pos.pos_ );
   return Iterator( matrix_, row+1UL, column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      matrix_.erase( first.row_, first.pos_ );
   }
   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Column<MT,false,false,false,CCAs...>::erase( Pred predicate )
{
   for( Iterator element=begin(); element!=end(); ++element ) {
      if( predicate( element->value() ) )
         matrix_.erase( element.row_, element.pos_ );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse column.
// The elements are selected by the given unary predicate \a predicate, which is expected
// to accept a single argument of the type of the elements and to be pure. The following
// example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( col2.begin(), col2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred >   // Type of the unary predicate
inline void Column<MT,false,false,false,CCAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   for( ; first!=last; ++first ) {
      if( predicate( first->value() ) )
         matrix_.erase( first.row_, first.pos_ );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::find( size_t index )
{
   const Iterator_t<MT> pos( matrix_.find( index, column() ) );

   if( pos != matrix_.end( index ) )
      return Iterator( matrix_, index, column(), pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::find( size_t index ) const
{
   const ConstIterator_t<MT> pos( matrix_.find( index, column() ) );

   if( pos != matrix_.end( index ) )
      return ConstIterator( matrix_, index, column(), pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::lowerBound( size_t index )
{
   for( size_t i=index; i<size(); ++i )
   {
      const Iterator_t<MT> pos( matrix_.find( i, column() ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, i, column(), pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::lowerBound( size_t index ) const
{
   for( size_t i=index; i<size(); ++i )
   {
      const ConstIterator_t<MT> pos( matrix_.find( i, column() ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, i, column(), pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::Iterator
   Column<MT,false,false,false,CCAs...>::upperBound( size_t index )
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const Iterator_t<MT> pos( matrix_.find( i, column() ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, i, column(), pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,false,CCAs...>::ConstIterator
   Column<MT,false,false,false,CCAs...>::upperBound( size_t index ) const
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const ConstIterator_t<MT> pos( matrix_.find( i, column() ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, i, column(), pos );
   }

   return end();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,false,false,false,CCAs...>&
   Column<MT,false,false,false,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,false,false,CCAs...>::canAlias( const Other* alias ) const
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column, \a false if not.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,false,false,CCAs...>::isAliased( const Other* alias ) const
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,false,false,CCAs...>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      matrix_(i,column()) = (*rhs)[i];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,false,false,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      for( ; i<element->index(); ++i )
         matrix_.erase( i, column() );
      matrix_(i++,column()) = element->value();
   }
   for( ; i<size(); ++i ) {
      matrix_.erase( i, column() );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a vector.
//
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline void Column<MT,false,false,false,CCAs...>::addAssign( const Vector<VT,false>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a vector.
//
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline void Column<MT,false,false,false,CCAs...>::subAssign( const Vector<VT,false>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC ROW-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Column for symmetric row-major sparse matrices.
// \ingroup column
//
// This specialization of Column adapts the class template to the requirements of symmetric
// row-major matrices.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
class Column<MT,false,false,true,CCAs...>
   : public View< SparseVector< Column<MT,false,false,true,CCAs...>, false > >
   , private ColumnData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnData<CCAs...>;                  //!< The type of the ColumnData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Column instance.
   using This = Column<MT,false,false,true,CCAs...>;

   //! Base type of this Column instance.
   using BaseType = View< SparseVector<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Column instance.
   using ResultType    = ColumnTrait_t<MT,CCAs...>;    //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Column&;                //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Column( MT& matrix, RCAs... args );

   Column( const Column& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Column() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Column& operator=( initializer_list<ElementType> list );
   inline Column& operator=( const Column& rhs );

   template< typename VT > inline Column& operator= ( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator= ( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator+=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator+=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator-=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator-=( const SparseVector<VT,false>& rhs );
   template< typename VT > inline Column& operator*=( const Vector<VT,false>&       rhs );
   template< typename VT > inline Column& operator/=( const DenseVector<VT,false>&  rhs );
   template< typename VT > inline Column& operator%=( const Vector<VT,false>&       rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::column;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t size() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   reserve( size_t n );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const ElementType& value );
   inline Iterator insert( size_t index, const ElementType& value );
   inline void     append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t index );
   inline Iterator erase( Iterator pos );
   inline Iterator erase( Iterator first, Iterator last );

   template< typename Pred, typename = DisableIf_t< IsIntegral_v<Pred> > >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Column& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,false>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,false>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the column.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for columns on row-major symmetric sparse matrices.
//
// \param matrix The matrix containing the column.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , size_t... CCAs >    // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Column<MT,false,false,true,CCAs...>::Column( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the column
{
   if( isChecked( args... ) ) {
      if( matrix_.columns() <= column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( column() < matrix_.columns(), "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Reference
   Column<MT,false,false,true,CCAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return matrix_(column(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstReference
   Column<MT,false,false,true,CCAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid column access index" );
   return const_cast<const MT&>( matrix_ )(column(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Reference
   Column<MT,false,false,true,CCAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the column elements.
//
// \param index Access index. The index must be smaller than the number of matrix rows.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid column access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstReference
   Column<MT,false,false,true,CCAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::begin()
{
   return matrix_.begin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::begin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the column.
//
// \return Iterator to the first element of the column.
//
// This function returns an iterator to the first element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::cbegin() const
{
   return matrix_.cbegin( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::end()
{
   return matrix_.end( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::end() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the column.
//
// \return Iterator just past the last element of the column.
//
// This function returns an iterator just past the last element of the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::cend() const
{
   return matrix_.cend( column() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all column elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// column by means of an initializer list. The column elements are assigned the values from the
// given initializer list. Missing values are reset to their default state. Note that in case
// the size of the initializer list exceeds the size of the column, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column" );
   }

   const InitializerVector<ElementType,false> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Column.
//
// \param rhs Sparse column to be copied.
// \return Reference to the assigned column.
// \exception std::invalid_argument Column sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two columns don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator=( const Column& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && column() == rhs.column() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Column sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( rhs.nonZeros() );
      assign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right );
      left.reset();
      assign( left, tmp );
   }
   else {
      left.reset();
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for sparse vectors.
//
// \param rhs Sparse vector to be assigned.
// \return Reference to the assigned column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT> tmp( right);
      left.reset();
      left.reserve( tmp.nonZeros() );
      assign( left, tmp );
   }
   else {
      left.reset();
      left.reserve( right.nonZeros() );
      assign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a dense vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator+=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a sparse vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be added to the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator+=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a dense vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator-=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a sparse vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be subtracted from the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator-=( const SparseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   left.reserve( tmp.nonZeros() );
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the sparse column.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator*=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the sparse column.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator/=( const DenseVector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the sparse column.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side vector
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::operator%=( const Vector<VT,false>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, 0UL, column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline MT& Column<MT,false,false,true,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the column.
//
// \return The matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline const MT& Column<MT,false,false,true,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the sparse column.
//
// \return The size of the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,true,CCAs...>::size() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column.
//
// \return The capacity of the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,true,CCAs...>::capacity() const noexcept
{
   return matrix_.capacity( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the column.
//
// \return The number of non-zero elements in the column.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of rows of the matrix containing the column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,true,CCAs...>::nonZeros() const
{
   return matrix_.nonZeros( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,true,CCAs...>::reset()
{
   matrix_.reset( column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column.
//
// \param n The new minimum capacity of the sparse column.
// \return void
//
// This function increases the capacity of the sparse column to at least \a n elements. The
// current values of the column elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
void Column<MT,false,false,true,CCAs...>::reserve( size_t n )
{
   matrix_.reserve( column(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new sparse column capacity.
//
// \return The new sparse column capacity.
//
// This function calculates a new column capacity based on the current capacity of the sparse
// column. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline size_t Column<MT,false,false,true,CCAs...>::extendCapacity() const noexcept
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = max( nonzeros, 7UL    );
   nonzeros = min( nonzeros, size() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the sparse column.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse column. In case the sparse column
// already contains an element with index \a index its value is modified, else a new element
// with the given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::set( size_t index, const ElementType& value )
{
   return matrix_.set( column(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse column access index.
//
// This function inserts a new element into the sparse column. However, duplicate elements
// are not allowed. In case the sparse column already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( column(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse column.
//
// \param index The index of the new element. The index must be smaller than the number of matrix rows.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column with elements. It appends
// a new element to the end of the sparse column without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse column
//  - the current number of non-zero elements must be smaller than the capacity of the column
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,true,CCAs...>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( column(), index, value, check );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline void Column<MT,false,false,true,CCAs...>::erase( size_t index )
{
   matrix_.erase( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::erase( Iterator pos )
{
   return matrix_.erase( column(), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse column.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( column(), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Column<MT,false,false,true,CCAs...>::erase( Pred predicate )
{
   matrix_.erase( column(), begin(), end(), predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse column. The
// elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements and to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto col2 = column( A, 2UL );
   col2.erase( col2.begin(), col2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Pred >   // Type of the unary predicate
inline void Column<MT,false,false,true,CCAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   matrix_.erase( column(), first, last, predicate );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::find( size_t index )
{
   return matrix_.find( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific column element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// column. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse column (the end() iterator) is returned. Note that
// the returned sparse column iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::find( size_t index ) const
{
   return matrix_.find( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::lowerBound( size_t index )
{
   return matrix_.lowerBound( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::Iterator
   Column<MT,false,false,true,CCAs...>::upperBound( size_t index )
{
   return matrix_.upperBound( column(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse column iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
inline typename Column<MT,false,false,true,CCAs...>::ConstIterator
   Column<MT,false,false,true,CCAs...>::upperBound( size_t index ) const
{
   return matrix_.upperBound( column(), index );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the column scaling.
// \return Reference to the sparse column.
//
// This function scales the column by applying the given scalar value \a scalar to each element
// of the column. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a column
// on a lower or upper unitriangular matrix. The attempt to scale such a column results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the scalar value
inline Column<MT,false,false,true,CCAs...>&
   Column<MT,false,false,true,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address can alias with the sparse column. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,false,true,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse column is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse column, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename Other >  // Data type of the foreign expression
inline bool Column<MT,false,false,true,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,false,true,CCAs...>::assign( const DenseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( column() ) == matrix_.capacity( column() ) )
         matrix_.reserve( column(), extendCapacity() );

      matrix_.append( column(), i, (*rhs)[i], true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,false,true,CCAs...>::assign( const SparseVector<VT,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      matrix_.append( column(), element->index(), element->value(), true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,false,true,CCAs...>::addAssign( const DenseVector<VT,false>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( column() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,false,true,CCAs...>::addAssign( const SparseVector<VT,false>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( column() );
   matrix_.reserve( column(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Column<MT,false,false,true,CCAs...>::subAssign( const DenseVector<VT,false>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( column() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CCAs >  // Compile time column arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Column<MT,false,false,true,CCAs...>::subAssign( const SparseVector<VT,false>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( column() );
   matrix_.reserve( column(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
