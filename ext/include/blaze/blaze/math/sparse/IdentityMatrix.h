//=================================================================================================
/*!
//  \file blaze/math/sparse/IdentityMatrix.h
//  \brief Implementation of an identity matrix
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

#ifndef _BLAZE_MATH_SPARSE_IDENTITYMATRIX_H_
#define _BLAZE_MATH_SPARSE_IDENTITYMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclStrLowTrait.h>
#include <blaze/math/traits/DeclStrUppTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUniLowTrait.h>
#include <blaze/math/traits/DeclUniUppTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/YieldsIdentity.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsSame.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup identity_matrix IdentityMatrix
// \ingroup sparse_matrix
*/
/*!\brief Efficient implementation of an \f$ N \times N \f$ identity matrix.
// \ingroup identity_matrix
//
// The IdentityMatrix class template is the representation of an immutable, arbitrary sized
// identity matrix with \f$ N \cdot N \f$ elements of arbitrary type. The type of the elements,
// the storage order, and the group tag of the matrix can be specified via the three template
// parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class IdentityMatrix;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the matrix elements. IdentityMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::defaultStorageOrder.
//  - Tag : optional type parameter to tag the matrix. The default type is \a blaze::Group0.
//          See \ref grouping_tagging for details.
//
// It is not possible to insert, erase or modify the elements of an identity matrix. It is only
// possible to read from the elements:

   \code
   using blaze::rowMajor;

   // Creating a row-major 4x4 identity matrix with 4 rows and 4 columns
   IdentityMatrix<double,rowMajor> A( 4 );

   // The function call operator provides access to all possible elements of the identity matrix,
   // including the zero elements.
   A(1,2) = 2.0;       // Compilation error: It is not possible to write to an identity matrix
   double d = A(2,1);  // Access to the element (2,1)

   // In order to traverse all non-zero elements currently stored in the matrix, the begin()
   // and end() functions can be used. In the example, all non-zero elements of the 2nd row
   // of A are traversed.
   for( IdentityMatrix<double,rowMajor>::Iterator i=A.begin(1); i!=A.end(1); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of IdentityMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, ...) can be performed on all possible combinations of row-major and column-major
// dense and sparse matrices with fitting element types. The following example gives an impression
// of the use of IdentityMatrix:

   \code
   using blaze::IdentityMatrix;
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   IdentityMatrix<double,rowMajor> A( 3 );  // Row-major 3x3 identity matrix

   DynamicMatrix<double,columnMajor> B( 3, 3 );  // Column-major 3x3 dynamic dense matrix
   CompressedMatrix<double,rowMajor> C( 3, 3 );  // Row-major 3x3 compressed sparse matrix
   CompressedMatrix<double,rowMajor> D( 3, 5 );  // Row-major 3x5 compressed sparse matrix
   // ... Initialization of B, C, and D

   DynamicMatrix<double,rowMajor>       E( A );  // Creation of a new row-major matrix as a copy of A
   CompressedMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;    // Addition of an identity matrix and a dense matrix
   E = C - A;    // Subtraction of a sparse matrix and an identity matrix
   F = A * D;    // Matrix multiplication between two matrices of different element types

   E = 2.0 * A;  // Scaling of an identity matrix
   F = A * 2.0;  // Scaling of an identity matrix
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
class IdentityMatrix
   : public Expression< SparseMatrix< IdentityMatrix<Type,SO,Tag>, SO > >
{
 public:
   //**Type definitions****************************************************************************
   using This     = IdentityMatrix<Type,SO,Tag>;          //!< Type of this IdentityMatrix instance.
   using BaseType = Expression< SparseMatrix<This,SO> >;  //!< Base type of this IdentityMatrix instance.

   //!< Result type for expression template evaluations.
   using ResultType = This;

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = IdentityMatrix<Type,!SO,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = IdentityMatrix<Type,!SO,Tag>;

   using ElementType    = Type;         //!< Type of the identity matrix elements.
   using TagType        = Tag;          //!< Tag type of this IdentityMatrix instance.
   using ReturnType     = const Type;   //!< Return type for expression template evaluations.
   using CompositeType  = const This&;  //!< Data type for composite expression templates.
   using Reference      = const Type;   //!< Reference to a identity matrix element.
   using ConstReference = const Type;   //!< Reference to a constant identity matrix element.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain an IdentityMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = IdentityMatrix<NewType,SO,Tag>;  //!< The type of the other IdentityMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a IdentityMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = IdentityMatrix<Type,SO,Tag>;  //!< The type of the other IdentityMatrix.
   };
   //**********************************************************************************************

   //**ConstIterator class definition**************************************************************
   /*!\brief Iterator over the elements of the identity matrix.
   */
   class ConstIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! Element type of the identity matrix.
      using Element = ValueIndexPair<Type>;

      using IteratorCategory = std::forward_iterator_tag;  //!< The iterator category.
      using ValueType        = Element;                    //!< Type of the underlying pointers.
      using PointerType      = ValueType*;                 //!< Pointer return type.
      using ReferenceType    = ValueType&;                 //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                  //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying pointers.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Default constructor**********************************************************************
      /*!\brief Default constructor for the ConstIterator class.
      */
      constexpr ConstIterator() noexcept
         : index_()  // Index to the current identity matrix element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the ConstIterator class.
      //
      // \param index Index to the initial matrix element.
      */
      constexpr ConstIterator( size_t index ) noexcept
         : index_( index )  // Index to the current identity matrix element
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      constexpr ConstIterator& operator++() noexcept {
         ++index_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      constexpr ConstIterator operator++( int ) noexcept {
         ConstIterator tmp( *this );
         ++index_;
         return tmp;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return The current value of the sparse element.
      */
      constexpr const Element operator*() const noexcept {
         return Element( Type(1), index_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse matrix element at the current iterator position.
      //
      // \return Reference to the sparse matrix element at the current iterator position.
      */
      constexpr const ConstIterator* operator->() const noexcept {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse element.
      //
      // \return The current value of the sparse element.
      */
      constexpr Type value() const noexcept {
         return Type(1);
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      constexpr size_t index() const noexcept {
         return index_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      constexpr bool operator==( const ConstIterator& rhs ) const noexcept {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      constexpr bool operator!=( const ConstIterator& rhs ) const noexcept {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two ConstIterator objects.
      //
      // \param rhs The right-hand side ConstIterator object.
      // \return The number of elements between the two ConstIterator objects.
      */
      constexpr DifferenceType operator-( const ConstIterator& rhs ) const noexcept {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      size_t index_;  //!< Index to the current identity matrix element.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   using Iterator = ConstIterator;  //!< Iterator over non-constant elements.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
            constexpr IdentityMatrix() noexcept;
   explicit constexpr IdentityMatrix( size_t n ) noexcept;

   template< typename MT, bool SO2 >
   inline IdentityMatrix( const Matrix<MT,SO2>& m );

   IdentityMatrix( const IdentityMatrix& ) = default;
   IdentityMatrix( IdentityMatrix&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~IdentityMatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline    ConstReference at( size_t i, size_t j ) const;
   constexpr ConstIterator  begin ( size_t i ) const noexcept;
   constexpr ConstIterator  cbegin( size_t i ) const noexcept;
   constexpr ConstIterator  end   ( size_t i ) const noexcept;
   constexpr ConstIterator  cend  ( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename MT, bool SO2 >
   inline IdentityMatrix& operator=( const Matrix<MT,SO2>& rhs ) &;

   IdentityMatrix& operator=( const IdentityMatrix& ) & = default;
   IdentityMatrix& operator=( IdentityMatrix&& ) & = default;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   constexpr size_t rows() const noexcept;
   constexpr size_t columns() const noexcept;
   constexpr size_t capacity() const noexcept;
   constexpr size_t capacity( size_t i ) const noexcept;
   constexpr size_t nonZeros() const noexcept;
   constexpr size_t nonZeros( size_t i ) const noexcept;
   constexpr void   clear() noexcept;
   constexpr void   resize( size_t n ) noexcept;
   constexpr void   swap( IdentityMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline ConstIterator find      ( size_t i, size_t j ) const;
   inline ConstIterator lowerBound( size_t i, size_t j ) const;
   inline ConstIterator upperBound( size_t i, size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   constexpr IdentityMatrix& transpose() noexcept;
   constexpr IdentityMatrix& ctranspose() noexcept;
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
   size_t n_;  //!< The current number of rows and columns of the identity matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE       ( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for IdentityMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr IdentityMatrix<Type,SO,Tag>::IdentityMatrix() noexcept
   : n_( 0UL )  // The current number of rows and columns of the identity matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an identity matrix of size \f$ N \times N \f$.
//
// \param n The number of rows and columns of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr IdentityMatrix<Type,SO,Tag>::IdentityMatrix( size_t n ) noexcept
   : n_( n )  // The current number of rows and columns of the identity matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor for different identity matrices.
//
// \param m Identity matrix to be copied.
// \exception std::invalid_argument Invalid setup of identity matrix.
//
// The matrix is sized according to the given \f$ N \times N \f$ identity matrix and
// initialized as a copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
template< typename MT     // Type of the foreign identity matrix
        , bool SO2 >      // Storage order of the foreign identity matrix
inline IdentityMatrix<Type,SO,Tag>::IdentityMatrix( const Matrix<MT,SO2>& m )
   : n_( (*m).rows() )  // The current number of rows and columns of the identity matrix
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( !IsIdentity_v<MT> && !isIdentity( *m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of identity matrix" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the identity matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr typename IdentityMatrix<Type,SO,Tag>::ConstReference
   IdentityMatrix<Type,SO,Tag>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid identity matrix column access index" );

   if( i == j )
      return Type( 1 );
   else
      return Type( 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
inline typename IdentityMatrix<Type,SO,Tag>::ConstReference
   IdentityMatrix<Type,SO,Tag>::at( size_t i, size_t j ) const
{
   if( i >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i+1UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return ConstIterator( i+1UL );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different identity matrices.
//
// \param rhs Identity matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid assignment to identity matrix.
//
// The matrix is resized according to the given \f$ N \times N \f$ identity matrix and
// initialized as a copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
template< typename MT     // Type of the right-hand side identity matrix
        , bool SO2 >      // Storage order of the right-hand side identity matrix
inline IdentityMatrix<Type,SO,Tag>&
   IdentityMatrix<Type,SO,Tag>::operator=( const Matrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( !IsIdentity_v<MT> && !isIdentity( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment of identity matrix" );
   }

   n_ = (*rhs).rows();

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the identity matrix.
//
// \return The number of rows of the identity matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::rows() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the identity matrix.
//
// \return The number of columns of the identity matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the identity matrix.
//
// \return The capacity of the identity matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::capacity() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row/column.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return 1UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the identity matrix
//
// \return The number of non-zero elements in the identity matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::nonZeros() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr size_t IdentityMatrix<Type,SO,Tag>::nonZeros( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < n_, "Invalid identity matrix row/column access index" );

   return 1UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the identity matrix.
//
// \return void
//
// After the clear() function, the size of the identity matrix is 0.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void IdentityMatrix<Type,SO,Tag>::clear() noexcept
{
   n_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the identity matrix.
//
// \param n The new number of rows and columns of the identity matrix.
// \return void
//
// This function resizes the matrix using the given size to \f$ n \times n \f$. Note that this
// function may invalidate all existing views (submatrices, rows, columns, ...) on the matrix if
// it is used to shrink the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
void constexpr IdentityMatrix<Type,SO,Tag>::resize( size_t n ) noexcept
{
   n_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two sparse matrices.
//
// \param m The identity matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void IdentityMatrix<Type,SO,Tag>::swap( IdentityMatrix& m ) noexcept
{
   const size_t tmp( n_ );
   n_ = m.n_;
   m.n_ = tmp;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific matrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// matrix. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an row/column iterator to the element.
// Otherwise an iterator just past the last non-zero element of row \a i or column \a j (the
// end() iterator) is returned.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
inline typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::find( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( i == j )
      return begin( i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less then the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
inline typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::lowerBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( ( !SO && j <= i ) || ( SO && i <= j ) )
      return begin( SO ? j : i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater then the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater then the given row
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
inline typename IdentityMatrix<Type,SO,Tag>::ConstIterator
   IdentityMatrix<Type,SO,Tag>::upperBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( SO  || i < rows()   , "Invalid identity matrix row access index"    );
   BLAZE_USER_ASSERT( !SO || j < columns(), "Invalid identity matrix column access index" );

   if( ( !SO && j < i ) || ( SO && i < j ) )
      return begin( SO ? j : i );
   else
      return end( SO ? j : i );
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief In-place transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr IdentityMatrix<Type,SO,Tag>& IdentityMatrix<Type,SO,Tag>::transpose() noexcept
{
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr IdentityMatrix<Type,SO,Tag>& IdentityMatrix<Type,SO,Tag>::ctranspose() noexcept
{
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the matrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address can alias with the matrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool IdentityMatrix<Type,SO,Tag>::canAlias( const Other* alias ) const noexcept
{
   MAYBE_UNUSED( alias );

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address is aliased with the matrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool IdentityMatrix<Type,SO,Tag>::isAliased( const Other* alias ) const noexcept
{
   MAYBE_UNUSED( alias );

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix can be used in SMP assignments.
//
// \return \a true in case the matrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the matrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the matrix).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
inline bool IdentityMatrix<Type,SO,Tag>::canSMPAssign() const noexcept
{
   return false;
}
//*************************************************************************************************








//=================================================================================================
//
//  IDENTITYMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name IdentityMatrix operators */
//@{
template< typename Type, bool SO, typename Tag >
constexpr void reset( IdentityMatrix<Type,SO,Tag>& m ) noexcept;

template< typename Type, bool SO, typename Tag >
constexpr void reset( IdentityMatrix<Type,SO,Tag>& m, size_t i ) noexcept;

template< typename Type, bool SO, typename Tag >
constexpr void clear( IdentityMatrix<Type,SO,Tag>& m ) noexcept;

template< RelaxationFlag RF, typename Type, bool SO, typename Tag >
constexpr bool isDefault( const IdentityMatrix<Type,SO,Tag>& m ) noexcept;

template< typename Type, bool SO, typename Tag >
constexpr bool isIntact( const IdentityMatrix<Type,SO,Tag>& m ) noexcept;

template< typename Type, bool SO, typename Tag >
constexpr void swap( IdentityMatrix<Type,SO,Tag>& a, IdentityMatrix<Type,SO,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void reset( IdentityMatrix<Type,SO,Tag>& m ) noexcept
{
   MAYBE_UNUSED( m );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given identity matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void reset( IdentityMatrix<Type,SO,Tag>& m, size_t i ) noexcept
{
   MAYBE_UNUSED( m, i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given identity matrix.
// \ingroup identity_matrix
//
// \param m The matrix to be cleared.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void clear( IdentityMatrix<Type,SO,Tag>& m ) noexcept
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given identity matrix is in default state.
// \ingroup identity_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the identity matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   blaze::IdentityMatrix<int> I;
   // ... Resizing and initialization
   if( isDefault( I ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( I ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the matrix
        , bool SO            // Storage order
        , typename Tag >     // Tag type
constexpr bool isDefault( const IdentityMatrix<Type,SO,Tag>& m ) noexcept
{
   return ( m.rows() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given identity matrix are intact.
// \ingroup identity_matrix
//
// \param m The identity matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the identity matrix are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::IdentityMatrix<int> I;
   // ... Resizing and initialization
   if( isIntact( I ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr bool isIntact( const IdentityMatrix<Type,SO,Tag>& m ) noexcept
{
   MAYBE_UNUSED( m );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two identity matrices.
// \ingroup identity_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr void swap( IdentityMatrix<Type,SO,Tag>& a, IdentityMatrix<Type,SO,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Declares the given matrix expression \a m as identity matrix.
// \ingroup identity_matrix
//
// \param m The input matrix.
// \return The redeclared matrix.
// \exception std::invalid_argument Invalid identity matrix specification.
//
// The \a declid function declares the given dense or sparse matrix expression \a m as identity
// matrix. In case the given matrix is not a square matrix, a \a std::invalid_argument exception
// is thrown.\n
// The following example demonstrates the use of the \a declid function:

   \code
   blaze::CompressedMatrix<double> A, B;
   // ... Resizing and initialization
   B = declid( A );
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
inline IdentityMatrix<ElementType_t<MT>,SO>
   declid( const Matrix<MT,SO>& m )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid identity matrix specification" );
   }

   return IdentityMatrix<ElementType_t<MT>,SO>( (*m).rows() );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, typename Tag >
struct IsSquare< IdentityMatrix<MT,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, typename Tag >
struct IsSymmetric< IdentityMatrix<MT,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, typename Tag >
struct IsHermitian< IdentityMatrix<MT,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, typename Tag >
struct IsUniLower< IdentityMatrix<MT,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, typename Tag >
struct IsUniUpper< IdentityMatrix<MT,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct IsResizable< IdentityMatrix<T,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SchurTraitEval1< T1, T2
                      , EnableIf_t< IsMatrix_v<T1> &&
                                    IsMatrix_v<T2> &&
                                    ( ( IsIdentity_v<T1> && IsIdentity_v<T2> ) ||
                                      ( IsIdentity_v<T1> && IsUniTriangular_v<T2> ) ||
                                      ( IsUniTriangular_v<T1> && IsIdentity_v<T2> ) ||
                                      ( IsUniLower_v<T1> && IsUniUpper_v<T2> ) ||
                                      ( IsUniUpper_v<T1> && IsUniLower_v<T2> ) ) &&
                                    !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( ( IsDenseMatrix_v<T1> && IsDenseMatrix_v<T2> ) ||
                                ( IsSparseMatrix_v<T1> && IsSparseMatrix_v<T2> )
                                ? SO1 && SO2
                                : ( IsSparseMatrix_v<T1>
                                    ? SO1
                                    : SO2 ) );

   using Type = IdentityMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                              , SO
                              , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   !IsIdentity_v<T1> && !IsZero_v<T1> && IsIdentity_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T1>, MultTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   IsIdentity_v<T1> && !IsIdentity_v<T2> && !IsZero_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T2>, MultTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   IsIdentity_v<T1> && IsIdentity_v<T2> > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = IdentityMatrix< AddTrait_t<MultType,MultType>
                              , StorageOrder_v<T1>
                              , AddTrait_t<MultTag,MultTag> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  KRONTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct KronTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   IsIdentity_v<T1> && IsIdentity_v<T2> > >
{
   using Type = IdentityMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                              , StorageOrder_v<T2>
                              , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
struct UnaryMapTraitEval1< T, OP
                         , EnableIf_t< IsMatrix_v<T> &&
                                       YieldsIdentity_v<OP,T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = IdentityMatrix< EvaluateTrait_t<ElementType>
                              , StorageOrder_v<T>
                              , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< IsMatrix_v<T1> &&
                                        IsMatrix_v<T2> &&
                                        YieldsIdentity_v<OP,T1,T2> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   using Type = IdentityMatrix< EvaluateTrait_t<ElementType>
                              , ( StorageOrder_v<T1> && StorageOrder_v<T2> )
                              , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSYMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclSymTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLHERMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclHermTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclLowTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNILOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclUniLowTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclStrLowTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = INVALID_TYPE;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclUppTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNIUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclUniUppTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclStrUppTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = INVALID_TYPE;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLDIAGTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct DeclDiagTrait< IdentityMatrix<T,SO,Tag> >
{
   using Type = IdentityMatrix<T,SO,Tag>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename Tag, typename T2 >
struct HighType< IdentityMatrix<T1,SO,Tag>, IdentityMatrix<T2,SO,Tag> >
{
   using Type = IdentityMatrix< typename HighType<T1,T2>::Type, SO, Tag >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool SO, typename Tag, typename T2 >
struct LowType< IdentityMatrix<T1,SO,Tag>, IdentityMatrix<T2,SO,Tag> >
{
   using Type = IdentityMatrix< typename LowType<T1,T2>::Type, SO, Tag >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
