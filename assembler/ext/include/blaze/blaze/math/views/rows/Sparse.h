//=================================================================================================
/*!
//  \file blaze/math/views/rows/Sparse.h
//  \brief Rows specialization for sparse matrices
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

#ifndef _BLAZE_MATH_VIEWS_ROWS_SPARSE_H_
#define _BLAZE_MATH_VIEWS_ROWS_SPARSE_H_


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
#include <blaze/math/constraints/Rows.h>
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
#include <blaze/math/traits/RowsTrait.h>
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
#include <blaze/math/views/rows/BaseTemplate.h>
#include <blaze/math/views/rows/RowsData.h>
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
//  CLASS TEMPLATE SPECIALIZATION FOR ROW-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on row-major sparse matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of row-major
// sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,true,false,SF,CRAs...>
   : public View< SparseMatrix< Rows<MT,true,false,SF,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Rows instance.
   using This = Rows<MT,true,false,SF,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< SparseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

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
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Rows& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Rows& operator=( const Rows& rhs );

   template< typename MT2, bool SO > inline Rows& operator= ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator+=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator-=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator%=( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t i, size_t nonzeros );
   inline void   trim();
   inline void   trim( size_t i );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const ElementType& value );
   inline Iterator insert  ( size_t i, size_t j, const ElementType& value );
   inline void     append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void     finalize( size_t i );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t i, Iterator pos );
   inline Iterator erase( size_t i, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t i, Iterator first, Iterator last, Pred predicate );
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
   inline Rows& transpose();
   inline Rows& ctranspose();

   template< typename Other > inline Rows& scale( const Other& scalar );
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
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE        ( MT );
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
/*!\brief Constructor for row selections on row-major sparse matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,true,false,SF,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Reference
   Rows<MT,true,false,SF,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstReference
   Rows<MT,true,false,SF,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Reference
   Rows<MT,true,false,SF,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstReference
   Rows<MT,true,false,SF,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.begin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.end( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cend( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cend( idx(i) );
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
// \exception std::invalid_argument Invalid assignment to row selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the row
// selection by means of an initializer list. The row selection elements are assigned the values
// from the given initializer list. Missing values are initialized as default. Note that if the
// size of the top-level initializer list does not match the number of rows of the row selection
// or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row selection" );
   }

   const InitializerMatrix<ElementType> tmp( list, columns() );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), i, 0UL ) ) {
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
/*!\brief Copy assignment operator for Rows.
//
// \param rhs Sparse row selection to be copied.
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse row selection is initialized as a copy of the given sparse row selection. In case
// the current sizes of the two row selections don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator=( const Rows& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( rhs, i, unchecked ), idx(i), 0UL ) ) {
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
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse row selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( right, i, unchecked ), idx(i), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the sparse row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator+=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the sparse row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator-=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
// \return Reference to the sparse row selection.
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::operator%=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,true,false,SF,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,true,false,SF,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row selection.
//
// \return The capacity of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::capacity() const noexcept
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::capacity( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.capacity( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse row selection.
//
// \return The number of non-zero elements in the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i )
      nonzeros += nonZeros( i );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.nonZeros( idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.reset( idx(i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value. Note that the
// capacity of the row remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::reset( size_t i )
{
   matrix_.reset( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row selection.
//
// \param nonzeros The new minimum capacity of the sparse row selection.
// \return void
//
// This function increases the capacity of the sparse row selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// rows are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::reserve( size_t nonzeros )
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
/*!\brief Setting the minimum capacity of a specific row of the row selection.
//
// \param i The row index of the new element \f$[0..M-1]\f$.
// \param nonzeros The new minimum capacity of the specified row.
// \return void
//
// This function increases the capacity of row \a i of the sparse row selection to at least
// \a nonzeros elements, but not beyond the current number of columns, respectively. The
// current values of the sparse row selection and all other individual row capacities are
// preserved.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,true,false,SF,CRAs...>::reserve( size_t i, size_t nonzeros )
{
   matrix_.reserve( idx(i), nonzeros );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all rows.
//
// \return void
//
// The trim() function can be used to reverse the effect of all row-specific reserve() calls.
// The function removes all excessive capacity from all rows. Note that this function does not
// remove the overall capacity but only reduces the capacity per row.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,true,false,SF,CRAs...>::trim()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      trim( i );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific row of the sparse matrix.
//
// \param i The index of the row to be trimmed (\f$[0..M-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a row-specific reserve() call. It removes
// all excessive capacity from the specified row. The excessive capacity is assigned to the
// subsequent row.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,true,false,SF,CRAs...>::trim( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   matrix_.trim( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new capacity for the specified sparse row.
//
// \param i The row index.
// \return The new sparse row capacity.
//
// This function calculates a new row capacity based on the current capacity of the specified
// sparse row. Note that the new capacity is restricted to the interval \f$[7..N]\f$.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,true,false,SF,CRAs...>::extendCapacity( size_t i ) const noexcept
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity( i )+1UL );
   nonzeros = max( nonzeros, 7UL );
   nonzeros = min( nonzeros, columns() );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity( i ), "Invalid capacity value" );

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
/*!\brief Setting an element of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse row selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return matrix_.set( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid row access index.
//
// This function inserts a new element into the sparse row selection. However, duplicate elements
// are not allowed. In case the sparse row selection already contains an element with row index
// \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return matrix_.insert( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified row of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row selection with elements. It
// appends a new element to the end of the specified row without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row of the sparse row selection
//  - the current number of non-zero elements in the row selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse row selection:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A( 42, 54 );
   auto B = rows( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );        // Finalizing row 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );        // Finalizing row 1
   B.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   B.finalize( 3 );        // Finalizing row 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a row.
//
// \param i The index of the row to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a row selection with
// elements. After completion of row \a i via the append() function, this function can be called
// to finalize row \a i and prepare the next row for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::finalize( size_t i )
{
   MAYBE_UNUSED( i );

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
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,true,false,SF,CRAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( idx(i), j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::erase( size_t i, Iterator pos )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.erase( idx(i), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from row \a i of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::erase( size_t i, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.erase( idx(i), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row selection. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,true,false,SF,CRAs...>::erase( Pred predicate )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.erase( idx(i), begin(i), end(i), predicate );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row selection.
//
// \param i The row index of the elements to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in row \a i of the sparse
// row selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( 2UL, R.begin(2UL), R.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,true,false,SF,CRAs...>::erase( size_t i, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   matrix_.erase( idx(i), first, last, predicate );
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
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::find( size_t i, size_t j )
{
   return matrix_.find( idx(i), j );
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
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::find( size_t i, size_t j ) const
{
   return matrix_.find( idx(i), j );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::lowerBound( size_t i, size_t j )
{
   return matrix_.lowerBound( idx(i), j );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::lowerBound( size_t i, size_t j ) const
{
   return matrix_.lowerBound( idx(i), j );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::Iterator
   Rows<MT,true,false,SF,CRAs...>::upperBound( size_t i, size_t j )
{
   return matrix_.upperBound( idx(i), j );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,true,false,SF,CRAs...>::ConstIterator
   Rows<MT,true,false,SF,CRAs...>::upperBound( size_t i, size_t j ) const
{
   return matrix_.upperBound( idx(i), j );
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
/*!\brief In-place transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief In-place conjugate transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief Scaling of the sparse row selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the sparse row selection.
//
// This function scales the row selection by applying the given scalar value \a scalar to each
// element of the row selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a row selection on a lower or upper unitriangular matrix. The attempt to scale such a
// row selection results in a compile time error!
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the scalar value
inline Rows<MT,true,false,SF,CRAs...>&
   Rows<MT,true,false,SF,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t i=0UL; i<rows(); ++i ) {
      const Iterator last( end(i) );
      for( Iterator element=begin(i); element!=last; ++element )
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
/*!\brief Returns whether the row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,true,false,SF,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,true,false,SF,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection can be used in SMP assignments.
//
// \return \a true in case the row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the row selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,true,false,SF,CRAs...>::canSMPAssign() const noexcept
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side dense matrix
        , bool SO >           // Storage order of the right-hand side dense matrix
inline void Rows<MT,true,false,SF,CRAs...>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      size_t remaining( matrix_.capacity( index ) );

      for( size_t j=0UL; j<columns(); ++j )
      {
         if( remaining == 0UL ) {
            matrix_.reserve( index, extendCapacity( i ) );
            remaining = matrix_.capacity( index ) - matrix_.nonZeros( index );
         }

         matrix_.append( index, j, (*rhs)(i,j), true );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,false,SF,CRAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t index( idx(i) );
      size_t remaining( matrix_.capacity( index ) );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
      {
         if( remaining == 0UL ) {
            matrix_.reserve( index, extendCapacity( i ) );
            remaining = matrix_.capacity( index ) - matrix_.nonZeros( index );
         }

         matrix_.append( index, element->index(), element->value(), true );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,true,false,SF,CRAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   // Counting the number of elements per row
   std::vector<size_t> rowLengths( rows(), 0UL );
   for( size_t j=0UL; j<columns(); ++j ) {
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         ++rowLengths[element->index()];
   }

   // Resizing the sparse matrix
   for( size_t i=0UL; i<rows(); ++i ) {
      reserve( i, rowLengths[i] );
   }

   // Appending the elements to the rows of the sparse row selection
   for( size_t j=0UL; j<columns(); ++j ) {
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         append( element->index(), j, element->value(), true );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,true,false,SF,CRAs...>::addAssign( const Matrix<MT2,SO>& rhs )
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,true,false,SF,CRAs...>::subAssign( const Matrix<MT2,SO>& rhs )
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,true,false,SF,CRAs...>::schurAssign( const Matrix<MT2,SO>& rhs )
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
//  CLASS TEMPLATE SPECIALIZATION FOR GENERAL COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on general column-major sparse matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of general
// column-major sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,false,false,false,CRAs...>
   : public View< SparseMatrix< Rows<MT,false,false,false,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the sparse matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Rows instance.
   using This = Rows<MT,false,false,false,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< SparseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

   //! Reference to a constant row value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant row value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;
   //**********************************************************************************************

   //**RowsElement class definition****************************************************************
   /*!\brief Access proxy for a specific element of the sparse row selection.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowsElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the RowsElement class.
      //
      // \param pos Iterator to the current position of the sparse row element.
      // \param column The column index.
      */
      inline RowsElement( IteratorType pos, size_t column )
         : pos_   ( pos    )  // Iterator to the current position of the sparse row element
         , column_( column )  // Index of the according column
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse row element.
      //
      // \param value The new value of the sparse row element.
      // \return Reference to the sparse row element.
      */
      template< typename T > inline RowsElement& operator=( const T& v ) {
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
      template< typename T > inline RowsElement& operator+=( const T& v ) {
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
      template< typename T > inline RowsElement& operator-=( const T& v ) {
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
      template< typename T > inline RowsElement& operator*=( const T& v ) {
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
      template< typename T > inline RowsElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse vector element at the current iterator position.
      //
      // \return Reference to the sparse vector element at the current iterator position.
      */
      inline const RowsElement* operator->() const {
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
      IteratorType pos_;     //!< Iterator to the current position of the sparse row element.
      size_t       column_;  //!< Index of the according column.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**RowsIterator class definition***************************************************************
   /*!\brief Iterator over the elements of a selected sparse row.
   */
   template< typename MatrixType      // Type of the sparse matrix
           , typename IteratorType >  // Type of the sparse matrix iterator
   class RowsIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;             //!< The iterator category.
      using ValueType        = RowsElement<MatrixType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                             //!< Pointer return type.
      using ReferenceType    = ValueType;                             //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                             //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the RowsIterator class.
      */
      inline RowsIterator()
         : matrix_( nullptr )  // The sparse matrix containing the selected row
         , row_   ( 0UL )      // The current row index
         , column_( 0UL )      // The current column index
         , pos_   ()           // Iterator to the current sparse element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the RowsIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      */
      inline RowsIterator( MatrixType& matrix, size_t row, size_t column )
         : matrix_( &matrix )  // The sparse matrix containing the selected row
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
      /*!\brief Constructor for the RowsIterator class.
      //
      // \param matrix The matrix containing the row.
      // \param row The row index.
      // \param column The column index.
      // \param pos Initial position of the iterator.
      */
      inline RowsIterator( MatrixType& matrix, size_t row, size_t column, IteratorType pos )
         : matrix_( &matrix )  // The sparse matrix containing the selected row
         , row_   ( row     )  // The current row index
         , column_( column  )  // The current column index
         , pos_   ( pos     )  // Iterator to the current sparse element
      {
         BLAZE_INTERNAL_ASSERT( matrix.find( row, column ) == pos, "Invalid initial iterator position" );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different RowsIterator instances.
      //
      // \param it The row iterator to be copied.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline RowsIterator( const RowsIterator<MatrixType2,IteratorType2>& it )
         : matrix_( it.matrix_ )  // The sparse matrix containing the selected row
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
      inline RowsIterator& operator++() {
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
      inline const RowsIterator operator++( int ) {
         const RowsIterator tmp( *this );
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
      /*!\brief Equality comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator==( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return column_ == rhs.column_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two RowsIterator objects.
      //
      // \param rhs The right-hand side row iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename MatrixType2, typename IteratorType2 >
      inline bool operator!=( const RowsIterator<MatrixType2,IteratorType2>& rhs ) const noexcept {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two row iterators.
      //
      // \param rhs The right-hand side row iterator.
      // \return The number of elements between the two row iterators.
      */
      inline DifferenceType operator-( const RowsIterator& rhs ) const {
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
      MatrixType*  matrix_;  //!< The sparse matrix containing the selected row.
      size_t       row_;     //!< The current row index.
      size_t       column_;  //!< The current column index.
      IteratorType pos_;     //!< Iterator to the current sparse element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename MatrixType2, typename IteratorType2 > friend class RowsIterator;
      template< typename MT2, bool SO2, bool DF2, bool SF2, typename... CRAs2 > friend class Rows;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = RowsIterator< const MT, ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, RowsIterator< MT, Iterator_t<MT> > >;
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
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Rows& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Rows& operator=( const Rows& rhs );

   template< typename MT2, bool SO > inline Rows& operator= ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator+=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator-=( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline Rows& operator%=( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t i, size_t nonzeros );
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
   inline void     finalize( size_t i );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t i, Iterator pos );
   inline Iterator erase( size_t i, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t i, Iterator first, Iterator last, Pred predicate );
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
   inline Rows& transpose();
   inline Rows& ctranspose();

   template< typename Other > inline Rows& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename MT2, bool SO > inline void assign     ( const DenseMatrix<MT2,SO>& rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,true>& rhs );
   template< typename MT2 >          inline void assign     ( const SparseMatrix<MT2,false>& rhs );
   template< typename MT2, bool SO > inline void addAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void subAssign  ( const Matrix<MT2,SO>& rhs );
   template< typename MT2, bool SO > inline void schurAssign( const Matrix<MT2,SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE            ( MT );
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
/*!\brief Constructor for row selections on general column-major sparse matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,false,false,false,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Reference
   Rows<MT,false,false,false,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstReference
   Rows<MT,false,false,false,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(idx(i),j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Reference
   Rows<MT,false,false,false,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstReference
   Rows<MT,false,false,false,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return Iterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return ConstIterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return ConstIterator( matrix_, idx(i), 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return Iterator( matrix_, idx(i), columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return ConstIterator( matrix_, idx(i), columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return ConstIterator( matrix_, idx(i), columns() );
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
// \exception std::invalid_argument Invalid assignment to row selection.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the row
// selection by means of an initializer list. The row selection elements are assigned the values
// from the given initializer list. Missing values are initialized as default. Note that if the
// size of the top-level initializer list does not match the number of rows of the row selection
// or the size of any nested list exceeds the number of columns, a \a std::invalid_argument
// exception is thrown. Also, if the underlying matrix \a MT is restricted and the assignment
// would violate an invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to row selection" );
   }

   const InitializerMatrix<ElementType> tmp( list, columns() );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), i, 0UL ) ) {
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
/*!\brief Copy assignment operator for Rows.
//
// \param rhs Sparse row selection to be copied.
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse row selection is initialized as a copy of the given sparse row selection. In case
// the current sizes of the two row selections don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator=( const Rows& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( rhs, i, unchecked ), idx(i), 0UL ) ) {
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
// \return Reference to the assigned row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The sparse row selection is initialized as a copy of the given matrix. In case the current
// sizes of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( right, i, unchecked ), idx(i), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be added to the row selection.
// \return Reference to the sparse row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator+=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
// \param rhs The right-hand side matrix to be subtracted from the row selection.
// \return Reference to the sparse row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator-=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
// \return Reference to the sparse row selection.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::operator%=( const Matrix<MT2,SO>& rhs )
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
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,false,false,false,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,false,false,false,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,false,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row selection.
//
// \return The capacity of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,false,CRAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,false,CRAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse row selection.
//
// \return The number of non-zero elements in the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,false,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<columns(); ++j ) {
      const auto end( matrix_.end( j ) );
      for( size_t i=0UL; i<rows(); ++i ) {
         auto pos = matrix_.find( idx(i), j );
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
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,false,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   size_t counter( 0UL );

   const ConstIterator last( end(i) );
   for( ConstIterator element=begin(i); element!=last; ++element ) {
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      reset( i );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value. Note that the
// capacity of the row remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::reset( size_t i )
{
   const size_t index( idx(i) );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( index+1UL )
                           :( index ) )
                        :( 0UL ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( index )
                           :( index+1UL ) )
                        :( columns() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      matrix_.erase( index, j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row selection.
//
// \param nonzeros The new minimum capacity of the sparse row selection.
// \return void
//
// This function increases the capacity of the sparse row selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// rows are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::reserve( size_t nonzeros )
{
   MAYBE_UNUSED( nonzeros );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific row of the row selection.
//
// \param i The row index of the new element \f$[0..M-1]\f$.
// \param nonzeros The new minimum capacity of the specified row.
// \return void
//
// This function increases the capacity of row \a i of the sparse row selection to at least
// \a nonzeros elements, but not beyond the current number of columns, respectively. The
// current values of the sparse row selection and all other individual row capacities are
// preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,false,CRAs...>::reserve( size_t i, size_t nonzeros )
{
   MAYBE_UNUSED( i, nonzeros );

   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all rows.
//
// \return void
//
// The trim() function can be used to reverse the effect of all row-specific reserve() calls.
// The function removes all excessive capacity from all rows. Note that this function does not
// remove the overall capacity but only reduces the capacity per row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,false,CRAs...>::trim()
{
   return;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific row of the sparse matrix.
//
// \param i The index of the row to be trimmed (\f$[0..M-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a row-specific reserve() call. It removes
// all excessive capacity from the specified row. The excessive capacity is assigned to the
// subsequent row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,false,CRAs...>::trim( size_t i )
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

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
/*!\brief Setting an element of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse row selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_, idx(i), j, matrix_.set( idx(i), j, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid row access index.
//
// This function inserts a new element into the sparse row selection. However, duplicate elements
// are not allowed. In case the sparse row selection already contains an element with row index
// \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return Iterator( matrix_, idx(i), j, matrix_.insert( idx(i), j, value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified row of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row selection with elements. It
// appends a new element to the end of the specified row without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row of the sparse row selection
//  - the current number of non-zero elements in the row selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse row selection:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A( 42, 54 );
   auto B = rows( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );        // Finalizing row 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );        // Finalizing row 1
   B.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   B.finalize( 3 );        // Finalizing row 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a row.
//
// \param i The index of the row to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a row selection with
// elements. After completion of row \a i via the append() function, this function can be called
// to finalize row \a i and prepare the next row for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::finalize( size_t i )
{
   MAYBE_UNUSED( i );

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
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,false,CRAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( idx(i), j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::erase( size_t i, Iterator pos )
{
   const size_t column( pos.column_ );

   if( column == columns() )
      return pos;

   matrix_.erase( column, pos.pos_ );
   return Iterator( matrix_, idx(i), column+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from row \a i of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::erase( size_t i, Iterator first, Iterator last )
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   for( ; first!=last; ++first ) {
      matrix_.erase( first.column_, first.pos_ );
   }

   return last;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row selection. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,false,false,false,CRAs...>::erase( Pred predicate )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      for( Iterator element=begin(i); element!=end(i); ++element ) {
         if( predicate( element->value() ) )
            matrix_.erase( element.column_, element.pos_ );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row selection.
//
// \param i The row index of the elements to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in row \a i of the sparse
// row selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( 2UL, R.begin(2UL), R.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,false,false,false,CRAs...>::erase( size_t i, Iterator first, Iterator last, Pred predicate )
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

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
/*!\brief Searches for a specific element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::find( size_t i, size_t j )
{
   const size_t index( idx(i) );
   const Iterator_t<MT> pos( matrix_.find( index, j ) );

   if( pos != matrix_.end( j ) )
      return Iterator( matrix_, index, j, pos );
   else
      return end( i );
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
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::find( size_t i, size_t j ) const
{
   const size_t index( idx(i) );
   const ConstIterator_t<MT> pos( matrix_.find( index, j ) );

   if( pos != matrix_.end( j ) )
      return ConstIterator( matrix_, index, j, pos );
   else
      return end( i );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::lowerBound( size_t i, size_t j )
{
   const size_t index( idx(i) );

   for( ; j<columns(); ++j )
   {
      const Iterator_t<MT> pos( matrix_.find( index, j ) );

      if( pos != matrix_.end( j ) )
         return Iterator( matrix_, index, j, pos );
   }

   return end( i );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::lowerBound( size_t i, size_t j ) const
{
   const size_t index( idx(i) );

   for( ; j<columns(); ++j )
   {
      const ConstIterator_t<MT> pos( matrix_.find( index, j ) );

      if( pos != matrix_.end( j ) )
         return ConstIterator( matrix_, index, j, pos );
   }

   return end( i );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::Iterator
   Rows<MT,false,false,false,CRAs...>::upperBound( size_t i, size_t j )
{
   return lowerBound( i, j+1UL );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,false,CRAs...>::ConstIterator
   Rows<MT,false,false,false,CRAs...>::upperBound( size_t i, size_t j ) const
{
   return lowerBound( i, j+1UL );
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
/*!\brief In-place transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the sparse row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::transpose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( trans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief In-place conjugate transpose of the row selection.
//
// \return Reference to the transposed row selection.
// \exception std::logic_error Invalid transpose of a non-quadratic matrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense row selection in-place. Note that this function can only
// be used for quadratic row selections, i.e. if the number of rows is equal to the number of
// columns. Also, the function fails if the invariants of an underlying, restricted matrix are
// violated. In all cases, a \a std::logic_error is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::ctranspose()
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic matrix" );
   }

   const ResultType tmp( ctrans( *this ) );

   if( IsRestricted_v<MT> ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( !tryAssign( matrix_, row( tmp, i, unchecked ), idx(i), 0UL ) ) {
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
/*!\brief Scaling of the sparse row selection by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the sparse row selection.
//
// This function scales the row selection by applying the given scalar value \a scalar to each
// element of the row selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used to
// scale a row selection on a lower or upper unitriangular matrix. The attempt to scale such a
// row selection results in a compile time error!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the scalar value
inline Rows<MT,false,false,false,CRAs...>&
   Rows<MT,false,false,false,CRAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   for( size_t j=0UL; j<columns(); ++j ) {
      const auto end( matrix_.end( j ) );
      for( size_t i=0UL; i<rows(); ++i ) {
         auto pos = matrix_.find( idx(i), j );
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
/*!\brief Returns whether the row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,false,false,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,false,false,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection can be used in SMP assignments.
//
// \return \a true in case the row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the row selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,false,false,CRAs...>::canSMPAssign() const noexcept
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side dense matrix
        , bool SO >           // Storage order of the right-hand side dense matrix
inline void Rows<MT,false,false,false,CRAs...>::assign( const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         RT value( (*rhs)(i,j) );
         if( !isDefault<strict>( value ) )
            matrix_.set( idx(i), j, std::move( value ) );
         else matrix_.erase( idx(i), j );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,false,false,CRAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t i=0UL; i<rows(); ++i ) {
      const size_t index( idx(i) );
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         RT value( element->value() );
         if( !isDefault<strict>( value ) )
            matrix_.set( index, element->index(), std::move( value ) );
         else matrix_.erase( index, element->index() );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2 >      // Type of the right-hand side sparse matrix
inline void Rows<MT,false,false,false,CRAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   using RT = If_t< IsComputation_v<MT2>, ElementType_t<MT>, const ElementType_t<MT2>& >;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         RT value( element->value() );
         if( !isDefault<strict>( value ) )
            matrix_.set( idx( element->index() ), j, std::move( value ) );
         else matrix_.erase( idx( element->index() ), j );
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,false,false,false,CRAs...>::addAssign( const Matrix<MT2,SO>& rhs )
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,false,false,false,CRAs...>::subAssign( const Matrix<MT2,SO>& rhs )
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
        , typename... CRAs >  // Compile time row arguments
template< typename MT2        // Type of the right-hand side matrix
        , bool SO >           // Storage order of the right-hand side matrix
inline void Rows<MT,false,false,false,CRAs...>::schurAssign( const Matrix<MT2,SO>& rhs )
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
//  CLASS TEMPLATE SPECIALIZATION FOR SYMMETRIC COLUMN-MAJOR SPARSE MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Rows for row selections on symmetric column-major sparse matrices.
// \ingroup rows
//
// This specialization of Rows adapts the class template to the requirements of symmetric
// column-major sparse matrices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
class Rows<MT,false,false,true,CRAs...>
   : public View< SparseMatrix< Rows<MT,false,false,true,CRAs...>, false > >
   , private RowsData<CRAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = RowsData<CRAs...>;                    //!< The type of the RowsData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Rows instance.
   using This = Rows<MT,false,false,true,CRAs...>;

   //! Base type of this Rows instance.
   using BaseType = View< SparseMatrix<This,false> >;

   using ViewedType    = MT;                           //!< The type viewed by this Rows instance.
   using ResultType    = RowsTrait_t<MT,N>;            //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;   //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;            //!< Type of the row elements.
   using ReturnType    = ReturnType_t<MT>;             //!< Return type for expression template evaluations
   using CompositeType = const Rows&;                  //!< Data type for composite expression templates.

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
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RRAs >
   explicit inline Rows( MT& matrix, RRAs... args );

   Rows( const Rows& ) = default;
   Rows( Rows&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Rows() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   Rows& operator=( const Rows& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idx;
   using DataType::idces;
   using DataType::rows;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t i, size_t nonzeros );
   inline void   trim();
   inline void   trim( size_t i );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const ElementType& value );
   inline Iterator insert  ( size_t i, size_t j, const ElementType& value );
   inline void     append  ( size_t i, size_t j, const ElementType& value, bool check=false );
   inline void     finalize( size_t i );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t i, size_t j );
   inline Iterator erase( size_t i, Iterator pos );
   inline Iterator erase( size_t i, Iterator first, Iterator last );

   template< typename Pred >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( size_t i, Iterator first, Iterator last, Pred predicate );
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
   Operand matrix_;  //!< The matrix containing the rows.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_SYMMETRIC_MATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ROWS_TYPE           ( MT );
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
/*!\brief Constructor for row selections on symmetric column-major sparse matrices.
//
// \param matrix The matrix containing the rows.
// \param args The runtime row arguments.
// \exception std::invalid_argument Invalid row access index.
//
// By default, the provided row arguments are checked at runtime. In case any row is not properly
// specified (i.e. if any specified index is greater than the number of rows of the given matrix)
// a \a std::invalid_argument exception is thrown. The checks can be skipped by providing the
// optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename... RRAs >  // Runtime row arguments
inline Rows<MT,false,false,true,CRAs...>::Rows( MT& matrix, RRAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the rows
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<rows(); ++i ) {
         if( matrix_.rows() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
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
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Reference
   Rows<MT,false,false,true,CRAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(j,idx(i));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the selected row elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstReference
   Rows<MT,false,false,true,CRAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(j,idx(i));
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Reference
   Rows<MT,false,false,true,CRAs...>::at( size_t i, size_t j )
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
/*!\brief Checked access to the selected row elements.
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstReference
   Rows<MT,false,false,true,CRAs...>::at( size_t i, size_t j ) const
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
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.begin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator to the first non-zero element of row \a i.
//
// This function returns a row iterator to the first non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cbegin( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.end( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cend( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last non-zero element of row \a i.
//
// This function returns an row iterator just past the last non-zero element of row \a i.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.cend( idx(i) );
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
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline MT& Rows<MT,false,false,true,CRAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the rows.
//
// \return The matrix containing the rows.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline const MT& Rows<MT,false,false,true,CRAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns of the row selection.
//
// \return The number of columns of the row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,true,CRAs...>::columns() const noexcept
{
   return matrix_.columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse row selection.
//
// \return The capacity of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,true,CRAs...>::capacity() const noexcept
{
   return nonZeros() + matrix_.capacity() - matrix_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,true,CRAs...>::capacity( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.capacity( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the sparse row selection.
//
// \return The number of non-zero elements in the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,true,CRAs...>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<rows(); ++i )
      nonzeros += nonZeros( i );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline size_t Rows<MT,false,false,true,CRAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.nonZeros( idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::reset()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.reset( idx(i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row to the default initial values.
//
// \param i The index of the row.
// \return void
//
// This function resets the values in the specified row to their default value. Note that the
// capacity of the row remains unchanged.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::reset( size_t i )
{
   matrix_.reset( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse row selection.
//
// \param nonzeros The new minimum capacity of the sparse row selection.
// \return void
//
// This function increases the capacity of the sparse row selection to at least \a nonzeros
// elements. The current values of the elements and the individual capacities of the selected
// rows are preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::reserve( size_t nonzeros )
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
/*!\brief Setting the minimum capacity of a specific row of the row selection.
//
// \param i The row index of the new element \f$[0..M-1]\f$.
// \param nonzeros The new minimum capacity of the specified row.
// \return void
//
// This function increases the capacity of row \a i of the sparse row selection to at least
// \a nonzeros elements, but not beyond the current number of columns, respectively. The
// current values of the sparse row selection and all other individual row capacities are
// preserved.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,true,CRAs...>::reserve( size_t i, size_t nonzeros )
{
   matrix_.reserve( idx(i), nonzeros );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all rows.
//
// \return void
//
// The trim() function can be used to reverse the effect of all row-specific reserve() calls.
// The function removes all excessive capacity from all rows. Note that this function does not
// remove the overall capacity but only reduces the capacity per row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,true,CRAs...>::trim()
{
   for( size_t i=0UL; i<rows(); ++i ) {
      trim( i );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific row of the sparse matrix.
//
// \param i The index of the row to be trimmed (\f$[0..M-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a row-specific reserve() call. It removes
// all excessive capacity from the specified row. The excessive capacity is assigned to the
// subsequent row.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
void Rows<MT,false,false,true,CRAs...>::trim( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   matrix_.trim( idx(i) );
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
/*!\brief Setting an element of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the sparse row selection. In case the sparse
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::set( size_t i, size_t j, const ElementType& value )
{
   return matrix_.set( j, idx(i), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid row access index.
//
// This function inserts a new element into the sparse row selection. However, duplicate elements
// are not allowed. In case the sparse row selection already contains an element with row index
// \a i and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::insert( size_t i, size_t j, const ElementType& value )
{
   return matrix_.insert( j, idx(i), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified row of the sparse row selection.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse row selection with elements. It
// appends a new element to the end of the specified row without any additional memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row of the sparse row selection
//  - the current number of non-zero elements in the row selection must be smaller than the
//    capacity of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a sparse row selection:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A( 42 );
   auto B = rows( A, { 10, 16, 4, 3 } );

   B.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   B.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );        // Finalizing row 0
   B.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );        // Finalizing row 1
   B.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   B.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   B.finalize( 3 );        // Finalizing row 3
   \endcode

// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::append( size_t i, size_t j, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      matrix_.insert( j, idx(i), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a row.
//
// \param i The index of the row to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a row selection with
// elements. After completion of row \a i via the append() function, this function can be called
// to finalize row \a i and prepare the next row for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::finalize( size_t i )
{
   MAYBE_UNUSED( i );

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
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline void Rows<MT,false,false,true,CRAs...>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.erase( j, idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::erase( size_t i, Iterator pos )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.erase( idx(i), pos );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse row selection.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from row \a i of the sparse row selection.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::erase( size_t i, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return matrix_.erase( idx(i), first, last );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse row selection.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse row selection. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,false,false,true,CRAs...>::erase( Pred predicate )
{
   for( size_t i=0UL; i<rows(); ++i ) {
      matrix_.erase( idx(i), begin(i), end(i), predicate );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse row selection.
//
// \param i The row index of the elements to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements in row \a i of the sparse
// row selection. The elements are selected by the given unary predicate \a predicate, which
// is expected to accept a single argument of the type of the elements and to be pure. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::SymmetricMatrix< blaze::CompressedMatrix<double,blaze::columnMajor> > A;
   // ... Resizing and initialization

   auto R = rows( A, { 4UL, 3UL, 5UL, 7UL } );
   R.erase( 2UL, R.begin(2UL), R.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Pred >     // Type of the unary predicate
inline void Rows<MT,false,false,true,CRAs...>::erase( size_t i, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   matrix_.erase( idx(i), first, last, predicate );
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
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::find( size_t i, size_t j )
{
   return matrix_.find( j, idx(i) );
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
// This function can be used to check whether a specific element is contained in the sparse row
// selection. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of row \a i (the end() iterator) is returned.
// Note that the returned iterator is subject to invalidation due to inserting operations via
// the function call operator, the set() function or the insert() function!
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::find( size_t i, size_t j ) const
{
   return matrix_.find( j, idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::lowerBound( size_t i, size_t j )
{
   return matrix_.lowerBound( j, idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::lowerBound( size_t i, size_t j ) const
{
   return matrix_.lowerBound( j, idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::Iterator
   Rows<MT,false,false,true,CRAs...>::upperBound( size_t i, size_t j )
{
   return matrix_.upperBound( j, idx(i) );
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
        , typename... CRAs >  // Compile time row arguments
inline typename Rows<MT,false,false,true,CRAs...>::ConstIterator
   Rows<MT,false,false,true,CRAs...>::upperBound( size_t i, size_t j ) const
{
   return matrix_.upperBound( j, idx(i) );
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
/*!\brief Returns whether the row selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address can alias with the sparse row selection. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,false,true,CRAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this row selection, \a false if not.
//
// This function returns whether the given address is aliased with the sparse row selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
template< typename Other >    // Data type of the foreign expression
inline bool Rows<MT,false,false,true,CRAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the row selection can be used in SMP assignments.
//
// \return \a true in case the row selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the row selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename MT         // Type of the sparse matrix
        , typename... CRAs >  // Compile time row arguments
inline bool Rows<MT,false,false,true,CRAs...>::canSMPAssign() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
