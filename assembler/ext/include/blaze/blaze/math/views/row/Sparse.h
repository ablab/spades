//=================================================================================================
/*!
//  \file blaze/math/views/row/Sparse.h
//  \brief Row specialization for sparse matrices
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

#ifndef _BLAZE_MATH_VIEWS_ROW_SPARSE_H_
#define _BLAZE_MATH_VIEWS_ROW_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/RowVector.h>
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
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RowTrait.h>
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
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/row/RowData.h>
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
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Row for rows on row-major sparse matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of row-major
// sparse matrices.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,true,false,SF,CRAs...>
   : public View< SparseVector< Row<MT,true,false,SF,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,true,false,SF,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< SparseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
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
   template< typename... RRAs >
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
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
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator= ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator+=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator-=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

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
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,true>& rhs );
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
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
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
/*!\brief Constructor for rows on row-major sparse matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,true,false,SF,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Reference
   Row<MT,true,false,SF,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstReference
   Row<MT,true,false,SF,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Reference
   Row<MT,true,false,SF,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstReference
   Row<MT,true,false,SF,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::begin()
{
   return matrix_.begin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::begin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::cbegin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::end()
{
   return matrix_.end( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::end() const
{
   return matrix_.cend( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::cend() const
{
   return matrix_.cend( row() );
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
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   const InitializerVector<ElementType,true> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Copy assignment operator for Row.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator=( const Row& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
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
// \return Reference to the assigned row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
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
// \return Reference to the assigned row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be added to the sparse row.
// \return Reference to the sparse row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator+=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be added to the sparse row.
// \return Reference to the sparse row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator+=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator-=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator-=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,true,false,SF,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,true,false,SF,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the sparse row.
//
// \return The size of the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,false,SF,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,false,SF,CRAs...>::capacity() const noexcept
{
   return matrix_.capacity( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,false,SF,CRAs...>::nonZeros() const
{
   return matrix_.nonZeros( row() );
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,true,false,SF,CRAs...>::reset()
{
   matrix_.reset( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
void Row<MT,true,false,SF,CRAs...>::reserve( size_t n )
{
   matrix_.reserve( row(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new sparse row capacity.
//
// \return The new sparse row capacity.
//
// This function calculates a new row capacity based on the current capacity of the sparse
// row. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,true,false,SF,CRAs...>::extendCapacity() const noexcept
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
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::set( size_t index, const ElementType& value )
{
   return matrix_.set( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,true,false,SF,CRAs...>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( row(), index, value, check );
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
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,true,false,SF,CRAs...>::erase( size_t index )
{
   matrix_.erase( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::erase( Iterator pos )
{
   return matrix_.erase( row(), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( row(), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Row<MT,true,false,SF,CRAs...>::erase( Pred predicate )
{
   matrix_.erase( row(), begin(), end(), predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse row. The elements
// are selected by the given unary predicate \a predicate, which is expected to accept a single
// argument of the type of the elements and to be pure. The following example demonstrates how
// to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( row2.begin(), row2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred >   // Type of the unary predicate
inline void Row<MT,true,false,SF,CRAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   matrix_.erase( row(), first, last, predicate );
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::find( size_t index )
{
   return matrix_.find( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::find( size_t index ) const
{
   return matrix_.find( row(), index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::lowerBound( size_t index )
{
   return matrix_.lowerBound( row(), index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( row(), index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::Iterator
   Row<MT,true,false,SF,CRAs...>::upperBound( size_t index )
{
   return matrix_.upperBound( row(), index );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,true,false,SF,CRAs...>::ConstIterator
   Row<MT,true,false,SF,CRAs...>::upperBound( size_t index ) const
{
   return matrix_.upperBound( row(), index );
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
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,true,false,SF,CRAs...>&
   Row<MT,true,false,SF,CRAs...>::scale( const Other& scalar )
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
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,true,false,SF,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,true,false,SF,CRAs...>::isAliased( const Other* alias ) const noexcept
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,true,false,SF,CRAs...>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<size(); ++j )
   {
      if( matrix_.nonZeros( row() ) == matrix_.capacity( row() ) )
         matrix_.reserve( row(), extendCapacity() );

      matrix_.append( row(), j, (*rhs)[j], true );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,false,SF,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      matrix_.append( row(), element->index(), element->value(), true );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,true,false,SF,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( row() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,false,SF,CRAs...>::addAssign( const SparseVector<VT,true>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( row() );
   matrix_.reserve( row(), tmp.nonZeros() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,true,false,SF,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( row() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,true,false,SF,CRAs...>::subAssign( const SparseVector<VT,true>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( row() );
   matrix_.reserve( row(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Row for general column-major sparse matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of general
// column-major sparse matrices.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,false,false,false,CRAs...>
   : public View< SparseVector< Row<MT,false,false,false,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,false,false,false,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< SparseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;
   //**********************************************************************************************

   //**RowElement class definition*****************************************************************
   /*!\brief Access proxy for a specific element of the sparse row.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the RowElement class.
      //
      // \param pos Iterator to the current position within the sparse row.
      // \param column The column index.
      */
      inline RowElement( IteratorType pos, size_t column )
         : pos_   ( pos    )  // Iterator to the current position within the sparse row
         , column_( column )  // Index of the according column
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse row element.
      //
      // \param value The new value of the sparse row element.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the addition.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the subtraction.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the multiplication.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse row element.
      //
      // \param value The right-hand side value for the division.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const RowElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse row element.
      //
      // \return The current value of the sparse row element.
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
         return column_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse row.
      size_t column_;     //!< Index of the according column.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**RowIterator class definition****************************************************************
   /*!\brief Iterator over the elements of the sparse row.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;            //!< The iterator category.
      using ValueType        = RowElement<MatrixType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                            //!< Pointer return type.
      using ReferenceType    = ValueType;                            //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                            //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the RowIterator class.
      */
      inline RowIterator()
         : matrix_( nullptr )  // The sparse matrix containing the row
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the RowIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      */
      inline RowIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the row
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {
         for( ; column_<matrix_->columns(); ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( column_ ) ) break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the RowIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator.
      */
      inline RowIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the row
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ( pos     )  // Iterator to the current sparse element
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline RowIterator( const RowIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the row
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
      inline RowIterator& operator++() {
         ++column_;
         for( ; column_<matrix_->columns(); ++column_ ) {
            pos_ = matrix_->find( row_, column_ );
            if( pos_ != matrix_->end( column_ ) ) break;
         }

         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const RowIterator operator++( int ) {
         const RowIterator tmp( *this );
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
         return ReferenceType( pos_, column_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, column_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ == rhs.column_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const RowIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowIterator& rhs ) const {
         size_t counter( 0UL );
         for( size_t j=rhs.column_; j<column_; ++j ) {
            if( matrix_->find( row_, j ) != matrix_->end( j ) )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      MatrixType*  matrix_;  //!< The sparse matrix containing the row.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class RowIterator;
      template< typename MT2, bool SO2, bool DF2, bool SF2, size_t... CRAs2 > friend class Row;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = RowIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, RowIterator< MT, Iterator_t<MT> > >;
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
   template< typename... RRAs >
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
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
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>& rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

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
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const Vector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const Vector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE ( MT );
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
/*!\brief Constructor for rows on column-major sparse matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
*/
template< typename MT         // Type of the sparse matrix
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,false,false,false,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Reference
   Row<MT,false,false,false,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstReference
   Row<MT,false,false,false,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(row(),index);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Reference
   Row<MT,false,false,false,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstReference
   Row<MT,false,false,false,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::begin()
{
   return Iterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::begin() const
{
   return ConstIterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::cbegin() const
{
   return ConstIterator( matrix_, row(), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::end()
{
   return Iterator( matrix_, row(), size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::end() const
{
   return ConstIterator( matrix_, row(), size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::cend() const
{
   return ConstIterator( matrix_, row(), size() );
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
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   const InitializerVector<ElementType,true> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Copy assignment operator for Row.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator=( const Row& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
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
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const CompositeType_t<VT> tmp( *rhs );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side vector to be added to the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator+=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator-=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,false,false,false,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,false,false,false,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the row.
//
// \return The size of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,false,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,false,CRAs...>::capacity() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,false,CRAs...>::nonZeros() const
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,false,CRAs...>::reset()
{
   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( row()+1UL )
                           :( row() ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( row() )
                           :( row()+1UL ) )
                        :( size() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_.erase( row(), j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
void Row<MT,false,false,false,CRAs...>::reserve( size_t n )
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
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::set( size_t index, const ElementType& value )
{
   return Iterator( matrix_, row(), index, matrix_.set( row(), index, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::insert( size_t index, const ElementType& value )
{
   return Iterator( matrix_, row(), index, matrix_.insert( row(), index, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,false,CRAs...>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( row(), index, value );
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
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,false,CRAs...>::erase( size_t index )
{
   matrix_.erase( row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::erase( Iterator pos )
{
   const size_t column( pos.column_ );

   if( column == size() )
      return pos;

   matrix_.erase( column, pos.pos_ );
   return Iterator( matrix_, row(), column+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      matrix_.erase( first.column_, first.pos_ );
   }
   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Row<MT,false,false,false,CRAs...>::erase( Pred predicate )
{
   for( Iterator element=begin(); element!=end(); ++element ) {
      if( predicate( element->value() ) )
         matrix_.erase( element.column_, element.pos_ );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse row. The
// elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( row2.begin(), row2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred >   // Type of the unary predicate
inline void Row<MT,false,false,false,CRAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   for( ; first!=last; ++first ) {
      if( predicate( first->value() ) )
         matrix_.erase( first.column_, first.pos_ );
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::find( size_t index )
{
   const Iterator_t<MT> pos( matrix_.find( row(), index ) );

   if( pos != matrix_.end( index ) )
      return Iterator( matrix_, row(), index, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::find( size_t index ) const
{
   const ConstIterator_t<MT> pos( matrix_.find( row(), index ) );

   if( pos != matrix_.end( index ) )
      return ConstIterator( matrix_, row(), index, pos );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::lowerBound( size_t index )
{
   for( size_t i=index; i<size(); ++i )
   {
      const Iterator_t<MT> pos( matrix_.find( row(), i ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, row(), i, pos );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::lowerBound( size_t index ) const
{
   for( size_t i=index; i<size(); ++i )
   {
      const ConstIterator_t<MT> pos( matrix_.find( row(), i ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, row(), i, pos );
   }

   return end();
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::Iterator
   Row<MT,false,false,false,CRAs...>::upperBound( size_t index )
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const Iterator_t<MT> pos( matrix_.find( row(), i ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, row(), i, pos );
   }

   return end();
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,false,CRAs...>::ConstIterator
   Row<MT,false,false,false,CRAs...>::upperBound( size_t index ) const
{
   for( size_t i=index+1UL; i<size(); ++i )
   {
      const ConstIterator_t<MT> pos( matrix_.find( row(), i ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, row(), i, pos );
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
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,false,false,false,CRAs...>&
   Row<MT,false,false,false,CRAs...>::scale( const Other& scalar )
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
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,false,false,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row, \a false if not.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,false,false,CRAs...>::isAliased( const Other* alias ) const noexcept
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,false,false,CRAs...>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( size_t j=0UL; j<(*rhs).size(); ++j ) {
      matrix_(row(),j) = (*rhs)[j];
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,false,false,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t j( 0UL );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      for( ; j<element->index(); ++j )
         matrix_.erase( row(), j );
      matrix_(row(),j++) = element->value();
   }
   for( ; j<size(); ++j ) {
      matrix_.erase( row(), j );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline void Row<MT,false,false,false,CRAs...>::addAssign( const Vector<VT,true>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline void Row<MT,false,false,false,CRAs...>::subAssign( const Vector<VT,true>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Row for symmetric column-major sparse matrices.
// \ingroup row
//
// This specialization of Row adapts the class template to the requirements of symmetric
// column-major sparse matrices.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
class Row<MT,false,false,true,CRAs...>
   : public View< SparseVector< Row<MT,false,false,true,CRAs...>, true > >
   , private RowData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowData<CRAs...>;                     //!< The type of the RowData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the dense matrix expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Row instance.
   using This = Row<MT,false,false,true,CRAs...>;

   //! Base type of this Row instance.
   using BaseType = View< SparseVector<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Row instance.
   using ResultType    = RowTrait_t<MT,CRAs...>;       //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Row&;                   //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
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
   template< typename... RRAs >
   explicit inline Row( MT& matrix, RRAs... args );

   Row( const Row& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Row() = default;
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
   inline Row& operator=( initializer_list<ElementType> list );
   inline Row& operator=( const Row& rhs );

   template< typename VT > inline Row& operator= ( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator= ( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator+=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator+=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator-=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator-=( const SparseVector<VT,true>& rhs );
   template< typename VT > inline Row& operator*=( const Vector<VT,true>&       rhs );
   template< typename VT > inline Row& operator/=( const DenseVector<VT,true>&  rhs );
   template< typename VT > inline Row& operator%=( const Vector<VT,true>&       rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;

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
   template< typename Other > inline Row& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT >    inline void assign   ( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void assign   ( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void addAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void addAssign( const SparseVector<VT,true>& rhs );
   template< typename VT >    inline void subAssign( const DenseVector <VT,true>& rhs );
   template< typename VT >    inline void subAssign( const SparseVector<VT,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the row.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE   ( MT );
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
/*!\brief Constructor for rows on column-major symmetric sparse matrices.
//
// \param matrix The matrix containing the row.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , size_t... CRAs >    // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Row<MT,false,false,true,CRAs...>::Row( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the row
{
   if( isChecked( args... ) ) {
      if( matrix_.rows() <= row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row() < matrix_.rows(), "Invalid row access index" );
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
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Reference
   Row<MT,false,false,true,CRAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return matrix_(index,row());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstReference
   Row<MT,false,false,true,CRAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid row access index" );
   return const_cast<const MT&>( matrix_ )(index,row());
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Reference
   Row<MT,false,false,true,CRAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the row elements.
//
// \param index Access index. The index must be smaller than the number of matrix columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid row access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstReference
   Row<MT,false,false,true,CRAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::begin()
{
   return matrix_.begin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::begin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the row.
//
// \return Iterator to the first element of the row.
//
// This function returns an iterator to the first element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::cbegin() const
{
   return matrix_.cbegin( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::end()
{
   return matrix_.end( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::end() const
{
   return matrix_.cend( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the row.
//
// \return Iterator just past the last element of the row.
//
// This function returns an iterator just past the last element of the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::cend() const
{
   return matrix_.cend( row() );
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
/*!\brief List assignment to all row elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to row.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the sparse
// row by means of an initializer list. The row elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the row, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is restricted and the assignment would violate
// an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row" );
   }

   const InitializerVector<ElementType,true> tmp( list, size() );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Copy assignment operator for Row.
//
// \param rhs Sparse row to be copied.
// \return Reference to the assigned row.
// \exception std::invalid_argument Row sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two rows don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator=( const Row& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Row sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), 0UL ) ) {
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
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
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
// \return Reference to the assigned row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<VT>, const VT& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be added to the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator+=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be added to the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator+=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side dense vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator-=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side sparse vector to be subtracted from the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower or upper triangular matrix and the
// assignment would violate its lower or upper property, respectively, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator-=( const SparseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \param rhs The right-hand side vector to be multiplied with the sparse row.
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator*=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( MultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator/=( const DenseVector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( DivType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( DivType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
// \return Reference to the sparse row.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side vector
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::operator%=( const Vector<VT,true>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( ResultType_t<VT> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( CrossType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), 0UL ) ) {
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
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline MT& Row<MT,false,false,true,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the row.
//
// \return The matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline const MT& Row<MT,false,false,true,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current size/dimension of the sparse row.
//
// \return The size of the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,true,CRAs...>::size() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row.
//
// \return The capacity of the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,true,CRAs...>::capacity() const noexcept
{
   return matrix_.capacity( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the row.
//
// \return The number of non-zero elements in the row.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of columns of the matrix containing the row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,true,CRAs...>::nonZeros() const
{
   return matrix_.nonZeros( row() );
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,true,CRAs...>::reset()
{
   matrix_.reset( row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row.
//
// \param n The new minimum capacity of the sparse row.
// \return void
//
// This function increases the capacity of the sparse row to at least \a n elements. The
// current values of the row elements are preserved.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
void Row<MT,false,false,true,CRAs...>::reserve( size_t n )
{
   matrix_.reserve( row(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new sparse row capacity.
//
// \return The new sparse row capacity.
//
// This function calculates a new row capacity based on the current capacity of the sparse
// row. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline size_t Row<MT,false,false,true,CRAs...>::extendCapacity() const
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
/*!\brief Setting an element of the sparse row.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse row. In case the sparse row already
// contains an element with index \a index its value is modified, else a new element with the
// given \a value is inserted.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::set( size_t index, const ElementType& value )
{
   return matrix_.set( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse row access index.
//
// This function inserts a new element into the sparse row. However, duplicate elements
// are not allowed. In case the sparse row already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::insert( size_t index, const ElementType& value )
{
   return matrix_.insert( index, row(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse row.
//
// \param index The index of the new element. The index must be smaller than the number of matrix columns.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row with elements. It appends
// a new element to the end of the sparse row without any memory allocation. Therefore it is
// strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse row
//  - the current number of non-zero elements must be smaller than the capacity of the row
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
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,true,CRAs...>::append( size_t index, const ElementType& value, bool check )
{
   matrix_.append( index, row(), value, check );
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
/*!\brief Erasing an element from the sparse row.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline void Row<MT,false,false,true,CRAs...>::erase( size_t index )
{
   matrix_.erase( index, row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::erase( Iterator pos )
{
   return matrix_.erase( row(), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse row.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::erase( Iterator first, Iterator last )
{
   return matrix_.erase( row(), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row. The elements are selected by the
// given unary predicate \a predicate, which is expected to accept a single argument of the type
// of the elements and to be pure. The following example demonstrates how to remove all elements
// that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Row<MT,false,false,true,CRAs...>::erase( Pred predicate )
{
   matrix_.erase( row(), begin(), end(), predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse row. The elements
// are selected by the given unary predicate \a predicate, which is expected to accept a single
// argument of the type of the elements and to be pure. The following example demonstrates how
// to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto row2 = row( A, 2UL );
   row2.erase( row2.begin(), row2.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Pred >   // Type of the unary predicate
inline void Row<MT,false,false,true,CRAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   matrix_.erase( row(), first, last, predicate );
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
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::find( size_t index )
{
   return matrix_.find( index, row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific row element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// row. It specifically searches for the element with index \a index. In case the element is
// found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse row (the end() iterator) is returned. Note that
// the returned sparse row iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::find( size_t index ) const
{
   return matrix_.find( index, row() );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::lowerBound( size_t index )
{
   return matrix_.lowerBound( index, row() );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::lowerBound( size_t index ) const
{
   return matrix_.lowerBound( index, row() );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::Iterator
   Row<MT,false,false,true,CRAs...>::upperBound( size_t index )
{
   return matrix_.upperBound( index, row() );
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
// pair of iterators specifying a range of indices. Note that the returned sparse row iterator
// is subject to invalidation due to inserting operations via the subscript operator, the set()
// function or the insert() function!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
inline typename Row<MT,false,false,true,CRAs...>::ConstIterator
   Row<MT,false,false,true,CRAs...>::upperBound( size_t index ) const
{
   return matrix_.upperBound( index, row() );
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
/*!\brief Scaling of the sparse row by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the row scaling.
// \return Reference to the sparse row.
//
// This function scales the row by applying the given scalar value \a scalar to each element
// of the row. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator. Note that the function cannot be used to scale a row
// on a lower or upper unitriangular matrix. The attempt to scale such a row results in a
// compile time error!
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the scalar value
inline Row<MT,false,false,true,CRAs...>&
   Row<MT,false,false,true,CRAs...>::scale( const Other& scalar )
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
/*!\brief Returns whether the sparse row can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address can alias with the sparse row. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,false,true,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse row is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse row, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the sparse matrix
        , size_t... CRAs >  // Compile time row arguments
template< typename Other >  // Data type of the foreign expression
inline bool Row<MT,false,false,true,CRAs...>::isAliased( const Other* alias ) const noexcept
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,false,true,CRAs...>::assign( const DenseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<size(); ++i )
   {
      if( matrix_.nonZeros( row() ) == matrix_.capacity( row() ) )
         matrix_.reserve( row(), extendCapacity() );

      matrix_.append( i, row(), (*rhs)[i], true );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,false,true,CRAs...>::assign( const SparseVector<VT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( ConstIterator_t<VT> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      matrix_.append( element->index(), row(), element->value(), true );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,false,true,CRAs...>::addAssign( const DenseVector<VT,true>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( row() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,false,true,CRAs...>::addAssign( const SparseVector<VT,true>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( AddType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   matrix_.reset( row() );
   matrix_.reserve( row(), tmp.nonZeros() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side dense vector
inline void Row<MT,false,false,true,CRAs...>::subAssign( const DenseVector<VT,true>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( row() );
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
        , size_t... CRAs >  // Compile time row arguments
template< typename VT >     // Type of the right-hand side sparse vector
inline void Row<MT,false,false,true,CRAs...>::subAssign( const SparseVector<VT,true>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( SubType );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE    ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   matrix_.reset( row() );
   matrix_.reserve( row(), tmp.nonZeros() );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
