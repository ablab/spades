//=================================================================================================
/*!
//  \file blaze/math/views/columns/Sparse.h
//  \brief Columns specialization for sparse matrices
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_SPARSE_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Columns.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/columns/BaseTemplate.h>
#include <blaze/math/views/columns/ColumnsData.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Columns for column selections on column-major sparse matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of column-major
// sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,true,false,SF,CCAs...>
   : public View< SparseMatrix< Columns<MT,true,false,SF,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,true,false,SF,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< SparseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

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
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Columns& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Columns& operator=( const Columns& rhs );

   template< typename MT2, bool SO > inline Columns& operator= ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator+=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator-=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator%=( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t j, size_t nonzeros );
   inline void   trim();
   inline void   trim( size_t j );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const ElementType& value );
   inline Iterator insert  ( size_t i, size_t j, const ElementType& value );
   inline void     append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void     finalize( size_t j );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t j, Iterator pos );
   inline Iterator erase( size_t j, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t j, Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t i, size_t j );
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline Iterator      lowerBound( size_t i, size_t j );
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline Iterator      upperBound( size_t i, size_t j );
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Columns& transpose();
   inline Columns& ctranspose();

   template< typename Other > inline Columns& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename MT2, bool SO > inline void assign     ( const DenseMatrix<MT2,SO>& rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2, bool SO > inline void addAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void subAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void schurAssign( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity( size_t j ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE        ( MT );
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
/*!\brief Constructor for column selections on column-major sparse matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not properly
// specified (i.e. if any specified index is greater than the number of columns of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,true,false,SF,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
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
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Reference
   Columns<MT,true,false,SF,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstReference
   Columns<MT,true,false,SF,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Reference
   Columns<MT,true,false,SF,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstReference
   Columns<MT,true,false,SF,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.begin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.end( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid row access index" );

   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid row access index" );

   return matrix_.cend( idx(j) );
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
/*!\brief List assignment to all elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the column
// selection by means of an initializer list. The column selection elements are assigned the
// values from the given initializer list. Missing values are initialized as default. Note that
// if the size of the top-level initializer list does not match the number of rows of the column
// selection or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column selection" );
   }

   const InitializerMatrix<ElementType> tmp( list, columns() );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, j ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Columns.
//
// \param rhs Sparse column selection to be copied.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse column selection is initialized as a copy of the given sparse column selection.
// In case the current sizes of the two column selections don't match, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular,
// or symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator=( const Columns& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && compareIndices( *this, rhs ) ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( rhs, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      left.reset();
      smpAssign( left, tmp );
   }
   else {
      left.reset();
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse column selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = CompositeType_t<MT2>;
   Right right( *rhs );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( right, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      left.reset();
      smpAssign( left, tmp );
   }
   else {
      left.reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side matrix to be for the Schur product.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::operator%=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

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
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,true,false,SF,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,true,false,SF,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column selection.
//
// \return The capacity of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::capacity() const noexcept
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::capacity( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.capacity( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse column selection.
//
// \return The number of non-zero elements in the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j )
      nonzeros += nonZeros( j );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.nonZeros( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.reset( idx(j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value. Note that the
// capacity of the column remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::reset( size_t j )
{
   matrix_.reset( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column selection.
//
// \param nonzeros The new minimum capacity of the sparse column selection.
// \return void
//
// This function increases the capacity of the sparse column selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// columns are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::reserve( size_t nonzeros )
{
   const size_t current( capacity() );

   if( nonzeros > current ) {
      matrix_.reserve( matrix_.capacity() + nonzeros - current );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific column of the column selection.
//
// \param j The column index of the new element \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified column.
// \return void
//
// This function increases the capacity of column \a j of the sparse column selection to at
// least \a nonzeros elements, but not beyond the current number of columns, respectively. The
// current values of the sparse column selection and all other individual column capacities
// are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,true,false,SF,CCAs...>::reserve( size_t j, size_t nonzeros )
{
   matrix_.reserve( idx(j), nonzeros );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all column-specific reserve() calls.
// The function removes all excessive capacity from all columns. Note that this function does not
// remove the overall capacity but only reduces the capacity per column.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,true,false,SF,CCAs...>::trim()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      trim( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific column of the sparse matrix.
//
// \param j The index of the column to be trimmed \f$[0..N-1]\f$.
// \return void
//
// This function can be used to reverse the effect of a column-specific reserve() call. It removes
// all excessive capacity from the specified column. The excessive capacity is assigned to the
// subsequent column.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,true,false,SF,CCAs...>::trim( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.trim( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new capacity for the specified sparse column.
//
// \param j The column index.
// \return The new sparse column capacity.
//
// This function calculates a new column capacity based on the current capacity of the specified
// sparse column. Note that the new capacity is restricted to the interval \f$[7..N]\f$.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,true,false,SF,CCAs...>::extendCapacity( size_t j ) const noexcept
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity( j )+1UL );
   nonzeros = max( nonzeros, 7UL );
   nonzeros = min( nonzeros, rows() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity( j ), "Invalid capacity value" );

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
/*!\brief Setting an element of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse column selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return matrix_.set( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid column access index.
//
// This function inserts a new element into the sparse column selection. However, duplicate
// elements are not allowed. In case the sparse column selection already contains an element
// with row index \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return matrix_.insert( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified column of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column selection with elements.
// It appends a new element to the end of the specified column without any additional memory
// allocation. Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified column of the sparse column selection
//  - the current number of non-zero elements in the column selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse column selection:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A( 42, 54 );
   auto B = columns( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in column 0 with column index 1
   B.finalize( 0 );        // Finalizing column 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in column 1 with column index 1
   B.finalize( 1 );        // Finalizing column 1
   B.finalize( 2 );        // Finalizing the empty column 2 to prepare column 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in column 3 with column index 0
   B.finalize( 3 );        // Finalizing column 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a column.
//
// \param j The index of the column to be finalized \f$[0..N-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a column selection with
// elements. After completion of column \a j via the append() function, this function can be
// called to finalize column \a j and prepare the next column for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::finalize( size_t j )
{
   MAYBE_UNUSED( j );

   return;
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
/*!\brief Erasing an element from the sparse column selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,true,false,SF,CCAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::erase( size_t j, Iterator pos )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.erase( idx(j), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from column \a j of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::erase( size_t j, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.erase( idx(j), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column selection. The elements are
// selected by the given unary predicate \a predicate, which is expected to accept a single
// argument of the type of the elements and to be pure. The following example demonstrates
// how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,true,false,SF,CCAs...>::erase( Pred predicate )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.erase( idx(j), begin(j), end(j), predicate );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column selection.
//
// \param j The column index of the elements to be erased. The index has to be in the range \f$[0..n-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in column \a j of the sparse
// column selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( 2UL, C.begin(2UL), C.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,true,false,SF,CCAs...>::erase( size_t j, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( idx(j), first, last, predicate );
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
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::find( size_t i, size_t j )
{
   return matrix_.find( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::find( size_t i, size_t j ) const
{
   return matrix_.find( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::lowerBound( size_t i, size_t j )
{
   return matrix_.lowerBound( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::lowerBound( size_t i, size_t j ) const
{
   return matrix_.lowerBound( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::Iterator
   Columns<MT,true,false,SF,CCAs...>::upperBound( size_t i, size_t j )
{
   return matrix_.upperBound( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,true,false,SF,CCAs...>::ConstIterator
   Columns<MT,true,false,SF,CCAs...>::upperBound( size_t i, size_t j ) const
{
   return matrix_.upperBound( i, idx(j) );
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
/*!\brief In-place transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::transpose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::ctranspose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the sparse column selection.
//
// This function scales the column selection by applying the given scalar value \a scalar to each
// element of the column selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to scale
// a column selection on a lower or upper unitriangular matrix. The attempt to scale such a column
// selection results in a compile time error!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the scalar value
inline Columns<MT,true,false,SF,CCAs...>&
   Columns<MT,true,false,SF,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t j=0UL; j<columns(); ++j ) {
      const Iterator last( end(j) );
      for( Iterator element=begin(j); element!=last; ++element )
         element->value() *= scalar;
   }

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
/*!\brief Returns whether the column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,true,false,SF,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,true,false,SF,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection can be used in SMP assignments.
//
// \return \a true in case the column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the column selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,true,false,SF,CCAs...>::canSMPAssign() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side dense matrix
        , bool SO >           // Storage order of the right-hand side dense matrix
inline void Columns<MT,true,false,SF,CCAs...>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      size_t remaining( matrix_.capacity( index ) );

      for( size_t i=0UL; i<rows(); ++i )
      {
         if( remaining == 0UL ) {
            matrix_.reserve( index, extendCapacity( j ) );
            remaining = matrix_.capacity( index ) - matrix_.nonZeros( index );
         }

         matrix_.append( i, index, (*rhs)(i,j), true );
         --remaining;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,false,SF,CCAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t index( idx(j) );
      size_t remaining( matrix_.capacity( index ) );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
      {
         if( remaining == 0UL ) {
            matrix_.reserve( index, extendCapacity( j ) );
            remaining = matrix_.capacity( index ) - matrix_.nonZeros( index );
         }

         matrix_.append( element->index(), index, element->value(), true );
         --remaining;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,true,false,SF,CCAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   // Counting the number of elements per column
   std::vector<size_t> columnLengths( columns(), 0UL );
   for( size_t i=0UL; i<rows(); ++i ) {
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         ++columnLengths[element->index()];
   }

   // Resizing the sparse matrix
   for( size_t j=0UL; j<columns(); ++j ) {
      reserve( j, columnLengths[j] );
   }

   // Appending the elements to the columns of the sparse column selection
   for( size_t i=0UL; i<rows(); ++i ) {
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         append( i, element->index(), element->value(), true );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a matrix.
//
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,true,false,SF,CCAs...>::addAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (*rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a matrix.
//
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,true,false,SF,CCAs...>::subAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (*rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a matrix.
//
// \param rhs The right-hand side matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,true,false,SF,CCAs...>::schurAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( SchurType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const SchurType tmp( serial( *this % (*rhs) ) );
   reset();
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
/*!\brief Specialization of Columns for column selections on general row-major sparse matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of general
// row-major sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,false,false,false,CCAs...>
   : public View< SparseMatrix< Columns<MT,false,false,false,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the sparse matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,false,false,false,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< SparseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

   //! Reference to a constant column value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant column value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;
   //**********************************************************************************************

   //**ColumnsElement class definition*************************************************************
   /*!\brief Access proxy for a specific element of the sparse column selection.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnsElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the ColumnsElement class.
      //
      // \param pos Iterator to the current position of the sparse column element.
      // \param row The row index.
      */
      inline ColumnsElement( IteratorType pos, size_t row )
         : pos_( pos )  // Iterator to the current position of the sparse column element
         , row_( row )  // Index of the according row
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse column element.
      //
      // \param value The new value of the sparse column element.
      // \return Reference to the sparse column element.
      */
      template< typename T > inline ColumnsElement& operator=( const T& v ) {
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
      template< typename T > inline ColumnsElement& operator+=( const T& v ) {
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
      template< typename T > inline ColumnsElement& operator-=( const T& v ) {
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
      template< typename T > inline ColumnsElement& operator*=( const T& v ) {
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
      template< typename T > inline ColumnsElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const ColumnsElement* operator->() const {
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
      IteratorType pos_;  //!< Iterator to the current position of the sparse column element.
      size_t       row_;  //!< Index of the according row.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**ColumnsIterator class definition************************************************************
   /*!\brief Iterator over the elements of a selected sparse column.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class ColumnsIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;                //!< The iterator category.
      using ValueType        = ColumnsElement<MatrixType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                                //!< Pointer return type.
      using ReferenceType    = ValueType;                                //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                                //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the ColumnsIterator class.
      */
      inline ColumnsIterator()
         : matrix_( nullptr )  // The sparse matrix containing the selected column
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ColumnsIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      */
      inline ColumnsIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the selected column
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
      /*!\brief Constructor for the ColumnsIterator class.
      //
      // \param matrix The matrix containing the column.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator.
      */
      inline ColumnsIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the selected column
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ( pos     )  // Iterator to the current sparse element
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ColumnsIterator instances.
      //
      // \param it The column iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline ColumnsIterator( const ColumnsIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the selected column
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
      inline ColumnsIterator& operator++() {
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
      inline const ColumnsIterator operator++( int ) {
         const ColumnsIterator tmp( *this );
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
      /*!\brief Equality comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return row_ == rhs.row_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ColumnsIterator objects.
      //
      // \param rhs The right-hand side column iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const ColumnsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two column iterators.
      //
      // \param rhs The right-hand side column iterator.
      // \return The number of elements between the two column iterators.
      */
      inline DifferenceType operator-( const ColumnsIterator& rhs ) const {
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
      MatrixType*  matrix_;  //!< The sparse matrix containing the selected column.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class ColumnsIterator;
      template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CCAs2 > friend class Columns;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = ColumnsIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, ColumnsIterator< MT, Iterator_t<MT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Columns& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Columns& operator=( const Columns& rhs );

   template< typename MT2, bool SO > inline Columns& operator= ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator+=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator-=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Columns& operator%=( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t j, size_t nonzeros );
   inline void   trim();
   inline void   trim( size_t j );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const ElementType& value );
   inline Iterator insert  ( size_t i, size_t j, const ElementType& value );
   inline void     append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void     finalize( size_t j );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t j, Iterator pos );
   inline Iterator erase( size_t j, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t j, Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t i, size_t j );
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline Iterator      lowerBound( size_t i, size_t j );
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline Iterator      upperBound( size_t i, size_t j );
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Columns& transpose();
   inline Columns& ctranspose();

   template< typename Other > inline Columns& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename MT2, bool SO > inline void assign     ( const DenseMatrix<MT2,SO>& rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,true>& rhs );
   template< typename MT2, bool SO > inline void addAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void subAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void schurAssign( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE         ( MT );
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
/*!\brief Constructor for column selections on general row-major sparse matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,false,false,false,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
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
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Reference
   Columns<MT,false,false,false,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstReference
   Columns<MT,false,false,false,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(i,idx(j));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Reference
   Columns<MT,false,false,false,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstReference
   Columns<MT,false,false,false,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return Iterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return ConstIterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return ConstIterator( matrix_, 0UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid row access index" );

   return Iterator( matrix_, rows(), idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return ConstIterator( matrix_, rows(), idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return ConstIterator( matrix_, rows(), idx(j) );
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
/*!\brief List assignment to all elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to column selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the column
// selection by means of an initializer list. The column selection elements are assigned the
// values from the given initializer list. Missing values are initialized as default. Note that
// if the size of the top-level initializer list does not match the number of rows of the column
// selection or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to column selection" );
   }

   const InitializerMatrix<ElementType> tmp( list, columns() );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, j ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Columns.
//
// \param rhs Sparse column selection to be copied.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse column selection is initialized as a copy of the given sparse column selection.
// In case the current sizes of the two column selections don't match, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular,
// or symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator=( const Columns& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && compareIndices( *this, rhs ) ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( rhs, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      left.reset();
      smpAssign( left, tmp );
   }
   else {
      left.reset();
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse column selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = CompositeType_t<MT2>;
   Right right( *rhs );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( right, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      if( IsSparseMatrix_v< ResultType_t<MT2> > )
         left.reset();
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseMatrix_v<MT2> )
         left.reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the column selection.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator+=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<AddType> ) {
      left.reset();
   }

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the column selection.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator-=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SubType> ) {
      left.reset();
   }

   smpAssign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side matrix to be for the Schur product.
// \return Reference to the sparse column selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::operator%=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      left.reset();
   }

   smpAssign( left, tmp );

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
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,false,false,false,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,false,false,false,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,false,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column selection.
//
// \return The capacity of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,false,CCAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,false,CCAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse column selection.
//
// \return The number of non-zero elements in the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,false,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i ) {
      const auto end( matrix_.end( i ) );
      for( size_t j=0UL; j<columns(); ++j ) {
         auto pos = matrix_.find( i, idx(j) );
         if( pos != end ) {
            ++nonzeros;
         }
      }
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,false,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   size_t counter( 0UL );

   const ConstIterator last( end(j) );
   for( ConstIterator element=begin(j); element!=last; ++element ) {
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
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      reset( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value. Note that the
// capacity of the column remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::reset( size_t j )
{
   const size_t index( idx(j) );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( index+1UL )
                           :( index ) )
                        :( 0UL ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( index )
                           :( index+1UL ) )
                        :( rows() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      matrix_.erase( i, index );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column selection.
//
// \param nonzeros The new minimum capacity of the sparse column selection.
// \return void
//
// This function increases the capacity of the sparse column selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// columns are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::reserve( size_t nonzeros )
{
   MAYBE_UNUSED( nonzeros );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific column of the column selection.
//
// \param j The column index of the new element \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified column.
// \return void
//
// This function increases the capacity of column \a j of the sparse column selection to at
// least \a nonzeros elements, but not beyond the current number of columns, respectively. The
// current values of the sparse column selection and all other individual column capacities are
// preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,false,CCAs...>::reserve( size_t j, size_t nonzeros )
{
   MAYBE_UNUSED( j, nonzeros );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all column-specific reserve() calls.
// The function removes all excessive capacity from all columns. Note that this function does not
// remove the overall capacity but only reduces the capacity per column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,false,CCAs...>::trim()
{
   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific column of the sparse matrix.
//
// \param j The index of the column to be trimmed (\f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a column-specific reserve() call. It removes
// all excessive capacity from the specified column. The excessive capacity is assigned to the
// subsequent column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,false,CCAs...>::trim( size_t j )
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

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
/*!\brief Setting an element of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse column selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_, i, idx(j), matrix_.set( i, idx(j), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid column access index.
//
// This function inserts a new element into the sparse column selection. However, duplicate
// elements are not allowed. In case the sparse column selection already contains an element
// with row index \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_, i, idx(j), matrix_.insert( i, idx(j), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified column of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column selection with elements.
// It appends a new element to the end of the specified column without any additional memory
// allocation. Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified column of the sparse column selection
//  - the current number of non-zero elements in the column selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse column selection:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 42, 54 );
   auto B = columns( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in column 0 with column index 1
   B.finalize( 0 );        // Finalizing column 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in column 1 with column index 1
   B.finalize( 1 );        // Finalizing column 1
   B.finalize( 2 );        // Finalizing the empty column 2 to prepare column 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in column 3 with column index 0
   B.finalize( 3 );        // Finalizing column 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( i, idx(j), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a column.
//
// \param j The index of the column to be finalized \f$[0..N-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a column selection with
// elements. After completion of column \a j via the append() function, this function can be called
// to finalize column \a j and prepare the next column for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::finalize( size_t j )
{
   MAYBE_UNUSED( j );

   return;
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
/*!\brief Erasing an element from the sparse column selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,false,CCAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( i, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::erase( size_t j, Iterator pos )
{
   const size_t row( pos.row_ );

   if( row == rows() )
      return pos;

   matrix_.erase( row, pos.pos_ );
   return Iterator( matrix_, row+1UL, idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from column \a j of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::erase( size_t j, Iterator first, Iterator last )
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   for( ; first!=last; ++first ) {
      matrix_.erase( first.row_, first.pos_ );
   }

   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column selection. The elements are
// selected by the given unary predicate \a predicate, which is expected to accept a single
// argument of the type of the elements and to be pure. The following example demonstrates
// how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,false,false,false,CCAs...>::erase( Pred predicate )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      for( Iterator element=begin(j); element!=end(j); ++element ) {
         if( predicate( element->value() ) )
            matrix_.erase( element.row_, element.pos_ );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column selection.
//
// \param j The column index of the elements to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in column \a j of the sparse
// column selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( 2UL, C.begin(2UL), C.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,false,false,false,CCAs...>::erase( size_t j, Iterator first, Iterator last, Pred predicate )
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

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
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::find( size_t i, size_t j )
{
   const size_t index( idx(j) );
   const Iterator_t<MT> pos( matrix_.find( i, index ) );

   if( pos != matrix_.end( i ) )
      return Iterator( matrix_, i, index, pos );
   else
      return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::find( size_t i, size_t j ) const
{
   const size_t index( idx(j) );
   const ConstIterator_t<MT> pos( matrix_.find( i, index ) );

   if( pos != matrix_.end( i ) )
      return ConstIterator( matrix_, i, index, pos );
   else
      return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::lowerBound( size_t i, size_t j )
{
   const size_t index( idx(j) );

   for( ; i<rows(); ++i )
   {
      const Iterator_t<MT> pos( matrix_.find( i, index ) );

      if( pos != matrix_.end( i ) )
         return Iterator( matrix_, i, index, pos );
   }

   return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::lowerBound( size_t i, size_t j ) const
{
   const size_t index( idx(j) );

   for( ; i<rows(); ++i )
   {
      const ConstIterator_t<MT> pos( matrix_.find( i, index ) );

      if( pos != matrix_.end( i ) )
         return ConstIterator( matrix_, i, index, pos );
   }

   return end( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::Iterator
   Columns<MT,false,false,false,CCAs...>::upperBound( size_t i, size_t j )
{
   return lowerBound( i+1UL, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,false,CCAs...>::ConstIterator
   Columns<MT,false,false,false,CCAs...>::upperBound( size_t i, size_t j ) const
{
   return lowerBound( i+1UL, j );
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
/*!\brief In-place transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the column selection.
//
// \return Reference to the transposed column selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense column selection in-place. Note that this function can only
// be used for quadratic column selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( !tryAssign( matrix_, column( tmp, j, unchecked ), 0UL, idx(j) ) ) {
            BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the sparse column selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the sparse column selection.
//
// This function scales the column selection by applying the given scalar value \a scalar to each
// element of the column selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to scale
// a column selection on a lower or upper unitriangular matrix. The attempt to scale such a column
// selection results in a compile time error!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the scalar value
inline Columns<MT,false,false,false,CCAs...>&
   Columns<MT,false,false,false,CCAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i ) {
      const auto end( matrix_.end( i ) );
      for( size_t j=0UL; j<columns(); ++j ) {
         auto pos = matrix_.find( i, idx(j) );
         if( pos != end ) {
            pos->value() *= scalar;
         }
      }
   }

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
/*!\brief Returns whether the column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,false,false,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,false,false,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection can be used in SMP assignments.
//
// \return \a true in case the column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the column selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,false,false,CCAs...>::canSMPAssign() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side dense matrix
        , bool SO >           // Storage order of the right-hand side dense matrix
inline void Columns<MT,false,false,false,CCAs...>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         RT value( (*rhs)(i,j) );
         if( !isDefault<strict>( value ) )
            matrix_.set( i, idx(j), std::move( value ) );
         else matrix_.erase( i, idx(j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,false,false,CCAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         RT value( element->value() );
         if( !isDefault<strict>( value ) )
            matrix_.set( i, idx( element->index() ), std::move( value ) );
         else matrix_.erase( i, idx( element->index() ) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Columns<MT,false,false,false,CCAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<columns(); ++j ) {
      const size_t index( idx(j) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         RT value( element->value() );
         if( !isDefault<strict>( value ) )
            matrix_.set( element->index(), index, std::move( value ) );
         else matrix_.erase( element->index(), index );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a matrix.
//
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,false,false,false,CCAs...>::addAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const AddType tmp( serial( *this + (*rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a matrix.
//
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,false,false,false,CCAs...>::subAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const SubType tmp( serial( *this - (*rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a matrix.
//
// \param rhs The right-hand side matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Columns<MT,false,false,false,CCAs...>::schurAssign( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE ( SchurType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const SchurType tmp( serial( *this % (*rhs) ) );
   reset();
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
/*!\brief Specialization of Columns for column selections on symmetric row-major sparse matrices.
// \ingroup columns
//
// This specialization of Columns adapts the class template to the requirements of symmetric
// column-major sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
class Columns<MT,false,false,true,CCAs...>
   : public View< SparseMatrix< Columns<MT,false,false,true,CCAs...>, true > >
   , private ColumnsData<CCAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ColumnsData<CCAs...>;                 //!< The type of the ColumnsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Columns instance.
   using This = Columns<MT,false,false,true,CCAs...>;

   //! Base type of this Columns instance.
   using BaseType = View< SparseMatrix<This,true> >;

   using ViewedType    = MT;                           //!< The type viewed by this Columns instance.
   using ResultType    = ColumnsTrait_t<MT,N>;         //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the column elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Columns&;               //!< Data type for composite expression templates.

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
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RCAs >
   explicit inline Columns( MT& matrix, RCAs... args );

   Columns( const Columns& ) = default;
   Columns( Columns&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Columns() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   Columns& operator=( const Columns& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t rows() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t j, size_t nonzeros );
   inline void   trim();
   inline void   trim( size_t j );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const ElementType& value );
   inline Iterator insert  ( size_t i, size_t j, const ElementType& value );
   inline void     append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void     finalize( size_t j );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t j, Iterator pos );
   inline Iterator erase( size_t j, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t j, Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t i, size_t j );
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline Iterator      lowerBound( size_t i, size_t j );
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline Iterator      upperBound( size_t i, size_t j );
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the columns.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COLUMNS_TYPE     ( MT );
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
/*!\brief Constructor for column selections on symmetric row-major sparse matrices.
//
// \param matrix The matrix containing the columns.
// \param args The runtime column arguments.
// \exception std::invalid_argument Invalid column access index.
//
// By default, the provided column arguments are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than the number of columns of the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename... RCAs >  // Runtime column arguments
inline Columns<MT,false,false,true,CCAs...>::Columns( MT& matrix, RCAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the columns
{
   if( isChecked( args... ) ) {
      for( size_t j=0UL; j<columns(); ++j ) {
         if( matrix_.columns() <= idx(j) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
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
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Reference
   Columns<MT,false,false,true,CCAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(j),i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstReference
   Columns<MT,false,false,true,CCAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(j),i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Reference
   Columns<MT,false,false,true,CCAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected column elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstReference
   Columns<MT,false,false,true,CCAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.begin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
//
// This function returns a column iterator to the first non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cbegin( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.end( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cend( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
//
// This function returns an column iterator just past the last non-zero element of column \a j.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.cend( idx(j) );
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
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline MT& Columns<MT,false,false,true,CCAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the columns.
//
// \return The matrix containing the columns.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline const MT& Columns<MT,false,false,true,CCAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows of the column selection.
//
// \return The number of rows of the column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,true,CCAs...>::rows() const noexcept
{
   return matrix_.rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse column selection.
//
// \return The capacity of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,true,CCAs...>::capacity() const noexcept
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
//
// This function returns the current capacity of the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,true,CCAs...>::capacity( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.capacity( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse column selection.
//
// \return The number of non-zero elements in the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,true,CCAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j )
      nonzeros += nonZeros( j );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
//
// This function returns the current number of non-zero elements in the specified column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline size_t Columns<MT,false,false,true,CCAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.nonZeros( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::reset()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.reset( idx(j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
//
// This function resets the values in the specified column to their default value. Note that the
// capacity of the column remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::reset( size_t j )
{
   matrix_.reset( idx(j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse column selection.
//
// \param nonzeros The new minimum capacity of the sparse column selection.
// \return void
//
// This function increases the capacity of the sparse column selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// columns are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::reserve( size_t nonzeros )
{
   const size_t current( capacity() );

   if( nonzeros > current ) {
      matrix_.reserve( matrix_.capacity() + nonzeros - current );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific column of the column selection.
//
// \param j The column index of the new element \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified column.
// \return void
//
// This function increases the capacity of column \a j of the sparse column selection to at least
// \a nonzeros elements, but not beyond the current number of columns, respectively. The current
// values of the sparse column selection and all other individual column capacities are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,true,CCAs...>::reserve( size_t j, size_t nonzeros )
{
   matrix_.reserve( idx(j), nonzeros );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all column-specific reserve() calls.
// The function removes all excessive capacity from all columns. Note that this function does not
// remove the overall capacity but only reduces the capacity per column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,true,CCAs...>::trim()
{
   for( size_t j=0UL; j<columns(); ++j ) {
      trim( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific column of the sparse matrix.
//
// \param j The index of the column to be trimmed (\f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a column-specific reserve() call. It
// removes all excessive capacity from the specified column. The excessive capacity is assigned
// to the subsequent column.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
void Columns<MT,false,false,true,CCAs...>::trim( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.trim( idx(j) );
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
/*!\brief Setting an element of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse column selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return matrix_.set( idx(j), i, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid column access index.
//
// This function inserts a new element into the sparse column selection. However, duplicate
// elements are not allowed. In case the sparse column selection already contains an element
// with row index \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return matrix_.insert( idx(j), i, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified column of the sparse column selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse column selection with elements.
// It appends a new element to the end of the specified column without any additional memory
// allocation. Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified column of the sparse column selection
//  - the current number of non-zero elements in the column selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse column selection:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::rowMajor> > A( 42 );
   auto B = columns( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in column 0 with column index 1
   B.finalize( 0 );        // Finalizing column 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in column 1 with column index 1
   B.finalize( 1 );        // Finalizing column 1
   B.finalize( 2 );        // Finalizing the empty column 2 to prepare column 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in column 3 with column index 0
   B.finalize( 3 );        // Finalizing column 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( idx(j), i, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a column.
//
// \param j The index of the column to be finalized \f$[0..N-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a column selection with
// elements. After completion of column \a j via the append() function, this function can be
// called to finalize column \a j and prepare the next column for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::finalize( size_t j )
{
   MAYBE_UNUSED( j );

   return;
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
/*!\brief Erasing an element from the sparse column selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline void Columns<MT,false,false,true,CCAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::erase( size_t j, Iterator pos )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.erase( idx(j), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse column selection.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from column \a j of the sparse column selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::erase( size_t j, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.erase( idx(j), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse column selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse column selection. The elements are
// selected by the given unary predicate \a predicate, which is expected to accept a single
// argument of the type of the elements and to be pure. The following example demonstrates
// how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::rowMajor> > A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,false,false,true,CCAs...>::erase( Pred predicate )
{
   for( size_t j=0UL; j<columns(); ++j ) {
      matrix_.erase( idx(j), begin(j), end(j), predicate );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse column selection.
//
// \param j The column index of the elements to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in column \a j of the sparse
// column selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::rowMajor> > A;
   // ... Resizing and initialization

   auto C = columns( A, { 4UL, 3UL, 5UL, 7UL } );
   C.erase( 2UL, C.begin(2UL), C.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Pred >     // Type of the unary predicate
inline void Columns<MT,false,false,true,CCAs...>::erase( size_t j, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( idx(j), first, last, predicate );
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
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::find( size_t i, size_t j )
{
   return matrix_.find( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse column
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::find( size_t i, size_t j ) const
{
   return matrix_.find( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::lowerBound( size_t i, size_t j )
{
   return matrix_.lowerBound( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// column index. In combination with the upperBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::lowerBound( size_t i, size_t j ) const
{
   return matrix_.lowerBound( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::Iterator
   Columns<MT,false,false,true,CCAs...>::upperBound( size_t i, size_t j )
{
   return matrix_.upperBound( idx(j), i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// column index. In combination with the lowerBound() function this function can be used to create
// a pair of iterators specifying a range of indices. Note that the returned iterator is subject
// to invalidation due to inserting operations via the function call operator, the set() function
// or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline typename Columns<MT,false,false,true,CCAs...>::ConstIterator
   Columns<MT,false,false,true,CCAs...>::upperBound( size_t i, size_t j ) const
{
   return matrix_.upperBound( idx(j), i );
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
/*!\brief Returns whether the column selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse column selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,false,true,CCAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this column selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse column selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
template< typename Other >    // Data type of the foreign expression
inline bool Columns<MT,false,false,true,CCAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the column selection can be used in SMP assignments.
//
// \return \a true in case the column selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the column selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , typename... CCAs >  // Compile time column arguments
inline bool Columns<MT,false,false,true,CCAs...>::canSMPAssign() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
