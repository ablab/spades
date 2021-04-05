//=================================================================================================
/*!
//  \file blaze/math/sparse/CompressedMatrix.h
//  \brief Implementation of a compressed MxN matrix
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

#ifndef _BLAZE_MATH_SPARSE_COMPRESSEDMATRIX_H_
#define _BLAZE_MATH_SPARSE_COMPRESSEDMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <utility>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/sparse/MatrixAccessProxy.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/ColumnsTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/ExpandTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/RowsTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/system/Thresholds.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/algorithms/Transfer.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameSize.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup compressed_matrix CompressedMatrix
// \ingroup sparse_matrix
*/
/*!\brief Efficient implementation of a \f$ M \times N \f$ compressed matrix.
// \ingroup compressed_matrix
//
// The CompressedMatrix class template is the representation of an arbitrary sized sparse
// matrix with \f$ M \cdot N \f$ dynamically allocated elements of arbitrary type. The type
// of the elements, the storage order, and the group tag of the matrix can be specified via
// the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class CompressedMatrix;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the matrix elements. CompressedMatrix can be used with
//          any non-cv-qualified, non-reference, non-pointer element type.
//  - SO  : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//          The default value is blaze::defaultStorageOrder.
//  - Tag : optional type parameter to tag the matrix. The default type is \a blaze::Group0.
//          See \ref grouping_tagging for details.
//
// Inserting/accessing elements in a compressed matrix can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Creating a row-major 4x3 compressed matrix with 4 rows and 3 columns
   CompressedMatrix<double,rowMajor> A( 4, 3 );

   // The function call operator provides access to all possible elements of the compressed matrix,
   // including the zero elements. In case the function call operator is used to access an element
   // that is currently not stored in the sparse matrix, the element is inserted into the matrix.
   A(1,2) = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the matrix it is inserted into the matrix, if it is already contained
   // in the matrix its value is modified.
   A.set( 2, 0, -1.2 );

   // An alternative for inserting elements into the matrix is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the matrix.
   A.insert( 2, 1, 3.7 );

   // A very efficient way to add new elements to a sparse matrix is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the specified row and that the matrix's capacity
   // is large enough to hold the new element.
   A.reserve( 3, 2 );       // Reserving space for 2 non-zero elements in row 3
   A.append( 3, 1, -2.1 );  // Appending the value -2.1 at column index 1 in row 3
   A.append( 3, 2,  1.4 );  // Appending the value 1.4 at column index 2 in row 3

   // The most efficient way to fill a (newly created) sparse matrix with elements, however, is
   // a combination of reserve(), append(), and the finalize() function.
   CompressedMatrix<double,rowMajor> B( 4, 3 );
   B.reserve( 3 );       // Reserving enough space for 3 non-zero elements
   B.append( 0, 1, 1 );  // Appending the value 1 in row 0 with column index 1
   B.finalize( 0 );      // Finalizing row 0
   B.append( 1, 1, 2 );  // Appending the value 2 in row 1 with column index 1
   B.finalize( 1 );      // Finalizing row 1
   B.append( 2, 0, 3 );  // Appending the value 3 in row 2 with column index 0
   B.finalize( 2 );      // Finalizing row 2

   // In order to traverse all non-zero elements currently stored in the matrix, the begin()
   // and end() functions can be used. In the example, all non-zero elements of the 2nd row
   // of A are traversed.
   for( CompressedMatrix<double,rowMajor>::Iterator i=A.begin(1); i!=A.end(1); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of CompressedMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of row-major and
// column-major dense and sparse matrices with fitting element types. The following example gives
// an impression of the use of CompressedMatrix:

   \code
   using blaze::CompressedMatrix;
   using blaze::DynamicMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   CompressedMatrix<double,rowMajor> A( 2, 3 );  // Empty row-major compressed double precision 2x3 matrix
   A(0,0) = 1.0; A(0,2) = 3.0; A(1,1) = 5.0;     // Element initialization

   CompressedMatrix<float,columnMajor> B( 2, 3 );  // Empty column-major compressed single precision 2x3 matrix
   B(0,1) = 3.0F; B(1,0) = 2.0F; B(1,2) = 6.0F;    // Element initialization

   DynamicMatrixMatrix<float> C( 2, 3, 4.0F );  // Directly, homogeneously initialized single precision dense 2x3 matrix
   CompressedMatrix<float>    D( 3, 2 );        // Empty row-major compressed single precision 3x2 matrix

   CompressedMatrix<double,rowMajor>    E( A );  // Creation of a new row-major matrix as a copy of A
   CompressedMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;     // Matrix addition and assignment to a row-major matrix
   E = A - C;     // Matrix subtraction and assignment to a column-major matrix
   F = A * D;     // Matrix multiplication between two matrices of different element types

   A *= 2.0;      // In-place scaling of matrix A
   E  = 2.0 * B;  // Scaling of matrix B
   F  = D * 2.0;  // Scaling of matrix D

   E += A - B;    // Addition assignment
   E -= A + C;    // Subtraction assignment
   F *= A * D;    // Multiplication assignment
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
class CompressedMatrix
   : public SparseMatrix< CompressedMatrix<Type,SO,Tag>, SO >
{
 private:
   //**Type definitions****************************************************************************
   using ElementBase  = ValueIndexPair<Type>;  //!< Base class for the compressed matrix element.
   using IteratorBase = ElementBase*;          //!< Iterator over non-constant base elements.
   //**********************************************************************************************

   //**Private class Element***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Value-index-pair for the CompressedMatrix class.
   //
   // This struct grants access to the data members of the base class and adapts the copy and
   // move semantics of the value-index-pair.
   */
   struct Element
      : public ElementBase
   {
      //**Constructors*****************************************************************************
      Element() = default;
      Element( const Element& rhs ) = default;
      Element( Element&& rhs ) = default;
      //*******************************************************************************************

      //**Assignment operators*********************************************************************
      inline Element& operator=( const Element& rhs )
      {
         this->value_ = rhs.value_;
         return *this;
      }

      inline Element& operator=( Element&& rhs )
      {
         this->value_ = std::move( rhs.value_ );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& rhs )
         -> EnableIf_t< IsSparseElement_v<Other>, Element& >
      {
         this->value_ = rhs.value();
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& rhs )
         -> EnableIf_t< IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( rhs.value() );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& v )
         -> EnableIf_t< !IsSparseElement_v<Other>, Element& >
      {
         this->value_ = v;
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& v )
         -> EnableIf_t< !IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( v );
         return *this;
      }
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      friend class CompressedMatrix;
      //*******************************************************************************************
   };
   /*! \endcond */
   //**********************************************************************************************

   //**Private class Uninitialized*****************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary helper class for the construction of compressed matrices.
   */
   struct Uninitialized {};
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This       = CompressedMatrix<Type,SO,Tag>;  //!< Type of this CompressedMatrix instance.
   using BaseType   = SparseMatrix<This,SO>;          //!< Base type of this CompressedMatrix instance.
   using ResultType = This;                           //!< Result type for expression template evaluations.

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = CompressedMatrix<Type,!SO,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = CompressedMatrix<Type,!SO,Tag>;

   using ElementType    = Type;                     //!< Type of the compressed matrix elements.
   using TagType        = Tag;                      //!< Tag type of this CompressedMatrix instance.
   using ReturnType     = const Type&;              //!< Return type for expression template evaluations.
   using CompositeType  = const This&;              //!< Data type for composite expression templates.
   using Reference      = MatrixAccessProxy<This>;  //!< Reference to a compressed matrix value.
   using ConstReference = const Type&;              //!< Reference to a constant compressed matrix value.
   using Iterator       = Element*;                 //!< Iterator over non-constant elements.
   using ConstIterator  = const Element*;           //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CompressedMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = CompressedMatrix<NewType,SO,Tag>;  //!< The type of the other CompressedMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CompressedMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = CompressedMatrix<Type,SO,Tag>;  //!< The type of the other CompressedMatrix.
   };
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
   inline CompressedMatrix();
   inline CompressedMatrix( size_t m, size_t n );
   inline CompressedMatrix( size_t m, size_t n, size_t nonzeros );
          CompressedMatrix( size_t m, size_t n, const std::vector<size_t>& nonzeros );
   inline CompressedMatrix( initializer_list< initializer_list<Type> > list );

   inline CompressedMatrix( const CompressedMatrix& sm );
   inline CompressedMatrix( CompressedMatrix&& sm ) noexcept;

   template< typename MT, bool SO2 > inline CompressedMatrix( const DenseMatrix<MT,SO2>&  dm );
   template< typename MT, bool SO2 > inline CompressedMatrix( const SparseMatrix<MT,SO2>& sm );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CompressedMatrix();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j ) noexcept;
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t i ) noexcept;
   inline ConstIterator  begin ( size_t i ) const noexcept;
   inline ConstIterator  cbegin( size_t i ) const noexcept;
   inline Iterator       end   ( size_t i ) noexcept;
   inline ConstIterator  end   ( size_t i ) const noexcept;
   inline ConstIterator  cend  ( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline CompressedMatrix& operator=( initializer_list< initializer_list<Type> > list ) &;
   inline CompressedMatrix& operator=( const CompressedMatrix& rhs ) &;
   inline CompressedMatrix& operator=( CompressedMatrix&& rhs ) & noexcept;

   template< typename MT, bool SO2 > inline CompressedMatrix& operator= ( const DenseMatrix<MT,SO2>&  rhs ) &;
   template< typename MT, bool SO2 > inline CompressedMatrix& operator= ( const SparseMatrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline CompressedMatrix& operator+=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline CompressedMatrix& operator-=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline CompressedMatrix& operator%=( const DenseMatrix<MT,SO2>&  rhs ) &;
   template< typename MT, bool SO2 > inline CompressedMatrix& operator%=( const SparseMatrix<MT,SO2>& rhs ) &;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows() const noexcept;
   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   inline void   clear();
          void   resize ( size_t m, size_t n, bool preserve=true );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t i, size_t nonzeros );
   inline void   trim   ();
   inline void   trim   ( size_t i );
   inline void   shrinkToFit();
   inline void   swap( CompressedMatrix& sm ) noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const Type& value );
   inline Iterator insert  ( size_t i, size_t j, const Type& value );
   inline void     append  ( size_t i, size_t j, const Type& value, bool check=false );
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
   inline CompressedMatrix& transpose();
   inline CompressedMatrix& ctranspose();

   template< typename Other > inline CompressedMatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename MT, bool SO2 > inline void assign     ( const DenseMatrix<MT,SO2>&  rhs );
   template< typename MT >           inline void assign     ( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT >           inline void assign     ( const SparseMatrix<MT,!SO>& rhs );
   template< typename MT, bool SO2 > inline void addAssign  ( const DenseMatrix<MT,SO2>&  rhs );
   template< typename MT, bool SO2 > inline void addAssign  ( const SparseMatrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline void subAssign  ( const DenseMatrix<MT,SO2>&  rhs );
   template< typename MT, bool SO2 > inline void subAssign  ( const SparseMatrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline void schurAssign( const DenseMatrix<MT,SO2>&  rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline CompressedMatrix( size_t m, size_t n, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const noexcept;
          void   reserveElements( size_t nonzeros );

   inline Iterator     castDown( IteratorBase it ) const noexcept;
   inline IteratorBase castUp  ( Iterator     it ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   Iterator insert( Iterator pos, size_t i, size_t j, const Type& value );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;         //!< The current number of rows of the compressed matrix.
   size_t n_;         //!< The current number of columns of the compressed matrix.
   size_t capacity_;  //!< The current capacity of the pointer array.
   Iterator* begin_;  //!< Pointers to the first non-zero element of each row.
   Iterator* end_;    //!< Pointers one past the last non-zero element of each row.

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE       ( ElementBase, Element );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
const Type CompressedMatrix<Type,SO,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for CompressedMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix()
   : m_       ( 0UL )      // The current number of rows of the compressed matrix
   , n_       ( 0UL )      // The current number of columns of the compressed matrix
   , capacity_( 0UL )      // The current capacity of the pointer array
   , begin_   ( nullptr )  // Pointers to the first non-zero element of each row
   , end_     ( nullptr )  // Pointers one past the last non-zero element of each row
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
//
// The matrix is initialized to the zero matrix and has no free capacity.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( size_t m, size_t n )
   : CompressedMatrix( m, n, Uninitialized() )
{
   for( size_t i=1UL; i<2UL*m_+2UL; ++i )
      begin_[i] = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nonzeros The number of expected non-zero elements.
//
// The matrix is initialized to the zero matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( size_t m, size_t n, size_t nonzeros )
   : CompressedMatrix( m, n, Uninitialized() )
{
   begin_[0UL] = allocate<Element>( nonzeros );
   for( size_t i=1UL; i<(2UL*m_+1UL); ++i )
      begin_[i] = begin_[0UL];
   end_[m_] = begin_[0UL]+nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nonzeros The expected number of non-zero elements in each row/column.
//
// The matrix is initialized to the zero matrix and will have the specified capacity in each
// row/column. Note that in case of a row-major matrix the given vector must have at least
// \a m elements, in case of a column-major matrix at least \a n elements.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
CompressedMatrix<Type,SO,Tag>::CompressedMatrix( size_t m, size_t n, const std::vector<size_t>& nonzeros )
   : CompressedMatrix( m, n, Uninitialized() )
{
   BLAZE_USER_ASSERT( nonzeros.size() == m, "Size of capacity vector and number of rows don't match" );

   size_t newCapacity( 0UL );
   for( auto it=nonzeros.begin(); it!=nonzeros.end(); ++it )
      newCapacity += *it;

   begin_[0UL] = end_[0UL] = allocate<Element>( newCapacity );
   for( size_t i=0UL; i<m_; ++i ) {
      begin_[i+1UL] = end_[i+1UL] = begin_[i] + nonzeros[i];
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all matrix elements.
//
// \param list The initializer list.
//
// This constructor provides the option to explicitly initialize the elements of the matrix by
// means of an initializer list:

   \code
   using blaze::rowMajor;

   blaze::CompressedMatrix<int,rowMajor> A{ { 1, 2, 3 },
                                            { 4, 5 },
                                            { 7, 8, 9 } };
   \endcode

// The matrix is sized according to the size of the initializer list and all its elements are
// initialized by the values of the given initializer list. Missing values are considered to
// be default values.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( initializer_list< initializer_list<Type> > list )
   : CompressedMatrix( list.size(), determineColumns( list ), blaze::nonZeros( list ) )
{
   size_t i( 0UL );

   for( const auto& rowList : list )
   {
      size_t j( 0UL );

      for( const Type& element : rowList ) {
         if( !isDefault<strict>( element ) )
            append( i, j, element );
         ++j;
      }

      finalize( i );
      ++i;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for CompressedMatrix.
//
// \param sm Sparse matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( const CompressedMatrix& sm )
   : CompressedMatrix( sm.m_, sm.n_, Uninitialized() )
{
   const size_t nonzeros( sm.nonZeros() );

   begin_[0UL] = allocate<Element>( nonzeros );
   for( size_t i=0UL; i<m_; ++i ) {
      end_[i] = castDown( std::copy( sm.begin(i), sm.end(i), castUp( begin_[i] ) ) );
      begin_[i+1UL] = end_[i];
   }
   end_[m_] = begin_[0UL]+nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for CompressedMatrix.
//
// \param sm The compressed matrix to be moved into this instance.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( CompressedMatrix&& sm ) noexcept
   : m_       ( sm.m_ )         // The current number of rows of the compressed matrix
   , n_       ( sm.n_ )         // The current number of columns of the compressed matrix
   , capacity_( sm.capacity_ )  // The current capacity of the pointer array
   , begin_   ( sm.begin_ )     // Pointers to the first non-zero element of each row
   , end_     ( sm.end_ )       // Pointers one past the last non-zero element of each row
{
   sm.m_        = 0UL;
   sm.n_        = 0UL;
   sm.capacity_ = 0UL;
   sm.begin_    = nullptr;
   sm.end_      = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from dense matrices.
//
// \param dm Dense matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign dense matrix
        , bool SO2 >      // Storage order of the foreign dense matrix
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( const DenseMatrix<MT,SO2>& dm )
   : CompressedMatrix( (*dm).rows(), (*dm).columns() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   assign( *this, *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different sparse matrices.
//
// \param sm Sparse matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign compressed matrix
        , bool SO2 >      // Storage order of the foreign compressed matrix
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( const SparseMatrix<MT,SO2>& sm )
   : CompressedMatrix( (*sm).rows(), (*sm).columns(), (*sm).nonZeros() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   assign( *this, *sm );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an uninitialized matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::CompressedMatrix( size_t m, size_t n, Uninitialized )
   : m_       ( m )                     // The current number of rows of the compressed matrix
   , n_       ( n )                     // The current number of columns of the compressed matrix
   , capacity_( m )                     // The current capacity of the pointer array
   , begin_( new Iterator[2UL*m+2UL] )  // Pointers to the first non-zero element of each row
   , end_  ( begin_+(m+1UL) )           // Pointers one past the last non-zero element of each row
{
   begin_[0] = nullptr;
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for CompressedMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>::~CompressedMatrix()
{
   if( begin_ != nullptr ) {
      deallocate( begin_[0UL] );
      delete[] begin_;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the compressed matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function returns a reference to the accessed value at position (\a i,\a j). In case the
// compressed matrix does not yet store an element at position (\a i,\a j) , a new element is
// inserted into the compressed matrix. Note that this function only performs an index check in
// case BLAZE_USER_ASSERT() is active. In contrast, the at() function is guaranteed to perform a
// check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Reference
   CompressedMatrix<Type,SO,Tag>::operator()( size_t i, size_t j ) noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return Reference( *this, i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief 2D-access to the compressed matrix elements.
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstReference
   CompressedMatrix<Type,SO,Tag>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const ConstIterator pos( lowerBound( i, j ) );

   if( pos == end_[i] || pos->index_ != j )
      return zero_;
   else
      return pos->value_;
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
// This function returns a reference to the accessed value at position (\a i,\a j). In case the
// compressed matrix does not yet store an element at position (\a i,\a j) , a new element is
// inserted into the compressed matrix. In contrast to the function call operator this function
// always performs a check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Reference
   CompressedMatrix<Type,SO,Tag>::at( size_t i, size_t j )
{
   if( i >= m_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
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
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstReference
   CompressedMatrix<Type,SO,Tag>::at( size_t i, size_t j ) const
{
   if( i >= m_ ) {
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::begin( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return begin_[i];
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return begin_[i];
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return begin_[i];
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::end( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return end_[i];
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return end_[i];
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid compressed matrix row access index" );
   return end_[i];
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief List assignment to all matrix elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::rowMajor;

   blaze::CompressedMatrix<int,rowMajor> A;
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix is resized according to the given initializer list and all its elements are
// assigned the values from the given initializer list. Missing values are considered to
// be default values.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator=( initializer_list< initializer_list<Type> > list ) &
{
   using blaze::nonZeros;

   resize( list.size(), determineColumns( list ), false );
   reserve( nonZeros( list ) );

   size_t i( 0UL );

   for( const auto& rowList : list )
   {
      size_t j( 0UL );

      for( const Type& element : rowList ) {
         if( !isDefault<strict>( element ) )
            append( i, j, element );
         ++j;
      }

      finalize( i );
      ++i;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for CompressedMatrix.
//
// \param rhs Sparse matrix to be copied.
// \return Reference to the assigned compressed matrix.
//
// The compressed matrix is resized according to the given compressed matrix and initialized
// as a copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator=( const CompressedMatrix& rhs ) &
{
   using std::swap;

   if( &rhs == this ) return *this;

   const size_t nonzeros( rhs.nonZeros() );

   if( rhs.m_ > capacity_ || nonzeros > capacity() )
   {
      Iterator* newBegin( new Iterator[2UL*rhs.m_+2UL] );
      Iterator* newEnd  ( newBegin+(rhs.m_+1UL) );

      newBegin[0UL] = allocate<Element>( nonzeros );
      for( size_t i=0UL; i<rhs.m_; ++i ) {
         newEnd[i] = castDown( std::copy( rhs.begin_[i], rhs.end_[i], castUp( newBegin[i] ) ) );
         newBegin[i+1UL] = newEnd[i];
      }
      newEnd[rhs.m_] = newBegin[0UL]+nonzeros;

      swap( begin_, newBegin );
      end_ = newEnd;
      capacity_ = rhs.m_;

      if( newBegin != nullptr ) {
         deallocate( newBegin[0UL] );
         delete[] newBegin;
      }
   }
   else {
      for( size_t i=0UL; i<rhs.m_; ++i ) {
         end_[i] = castDown( std::copy( rhs.begin_[i], rhs.end_[i], castUp( begin_[i] ) ) );
         begin_[i+1UL] = end_[i];
      }
   }

   m_ = rhs.m_;
   n_ = rhs.n_;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for CompressedMatrix.
//
// \param rhs The compressed matrix to be moved into this instance.
// \return Reference to the assigned compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator=( CompressedMatrix&& rhs ) & noexcept
{
   if( begin_ != nullptr ) {
      deallocate( begin_[0UL] );
      delete[] begin_;
   }

   m_        = rhs.m_;
   n_        = rhs.n_;
   capacity_ = rhs.capacity_;
   begin_    = rhs.begin_;
   end_      = rhs.end_;

   rhs.m_        = 0UL;
   rhs.n_        = 0UL;
   rhs.capacity_ = 0UL;
   rhs.begin_    = nullptr;
   rhs.end_      = nullptr;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for dense matrices.
//
// \param rhs Dense matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator=( const DenseMatrix<MT,SO2>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).canAlias( this ) ) {
      CompressedMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      assign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different sparse matrices.
//
// \param rhs Sparse matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO2 >      // Storage order of the right-hand side compressed matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator=( const SparseMatrix<MT,SO2>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).canAlias( this ) ||
       (*rhs).rows()     > capacity_ ||
       (*rhs).nonZeros() > capacity() ) {
      CompressedMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      reset();

      if( !IsZero_v<MT> ) {
         assign( *this, *rhs );
      }
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator+=( const Matrix<MT,SO2>& rhs ) &
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      addAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator-=( const Matrix<MT,SO2>& rhs ) &
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      subAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Schur product assignment operator for the multiplication of a dense matrix
//        (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator%=( const DenseMatrix<MT,SO2>& rhs ) &
{
   using blaze::schurAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      CompressedMatrix tmp( *this % (*rhs) );
      swap( tmp );
   }
   else {
      CompositeType_t<MT> tmp( *rhs );
      schurAssign( *this, tmp );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Schur product assignment operator for the multiplication of a sparse matrix
//        (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side sparse matrix
        , bool SO2 >      // Storage order of the right-hand side sparse matrix
inline CompressedMatrix<Type,SO,Tag>&
   CompressedMatrix<Type,SO,Tag>::operator%=( const SparseMatrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      CompressedMatrix tmp( *this % (*rhs) );
      swap( tmp );
   }
   else {
      reset();
   }

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the compressed matrix.
//
// \return The number of rows of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::rows() const noexcept
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the compressed matrix.
//
// \return The number of columns of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the compressed matrix.
//
// \return The capacity of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::capacity() const noexcept
{
   if( begin_ != nullptr )
      return end_[m_] - begin_[0UL];
   else return 0UL;
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
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::capacity( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return begin_[i+1UL] - begin_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the compressed matrix
//
// \return The number of non-zero elements in the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<m_; ++i )
      nonzeros += nonZeros( i );

   return nonzeros;
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
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return end_[i] - begin_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::reset()
{
   for( size_t i=0UL; i<m_; ++i )
      end_[i] = begin_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// This function resets the values in the specified row/column to their default value. In case
// the storage order is set to \a rowMajor the function resets the values in row \a i, in case
// the storage order is set to \a columnMajor the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::reset( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   end_[i] = begin_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the compressed matrix.
//
// \return void
//
// After the clear() function, the size of the compressed matrix is 0.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::clear()
{
   if( end_ != nullptr )
      end_[0UL] = end_[m_];
   m_ = 0UL;
   n_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the compressed matrix.
//
// \param m The new number of rows of the compressed matrix.
// \param n The new number of columns of the compressed matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Note that this function may invalidate all existing views (submatrices, rows, columns,
// ...) on the matrix if it is used to shrink the matrix. Additionally, the resize operation
// potentially changes all matrix elements. In order to preserve the old matrix values, the
// \a preserve flag can be set to \a true.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
void CompressedMatrix<Type,SO,Tag>::resize( size_t m, size_t n, bool preserve )
{
   using std::swap;

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( begin_ == nullptr || size_t( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );

   if( m == m_ && n == n_ ) return;

   if( begin_ == nullptr )
   {
      begin_ = new Iterator[2UL*m+2UL];
      end_   = begin_+m+1UL;

      for( size_t i=0UL; i<2UL*m+2UL; ++i ) {
         begin_[i] = nullptr;
      }

      capacity_ = m;
   }
   else if( m > capacity_ )
   {
      Iterator* newBegin( new Iterator[2UL*m+2UL] );
      Iterator* newEnd  ( newBegin+m+1UL );

      newBegin[0UL] = begin_[0UL];

      if( preserve ) {
         for( size_t i=0UL; i<m_; ++i ) {
            newEnd  [i]     = end_  [i];
            newBegin[i+1UL] = begin_[i+1UL];
         }
         for( size_t i=m_; i<m; ++i ) {
            newBegin[i+1UL] = newEnd[i] = begin_[m_];
         }
      }
      else {
         for( size_t i=0UL; i<m; ++i ) {
            newBegin[i+1UL] = newEnd[i] = begin_[0UL];
         }
      }

      newEnd[m] = end_[m_];

      swap( newBegin, begin_ );
      delete[] newBegin;
      end_ = newEnd;
      capacity_ = m;
   }
   else if( m > m_ )
   {
      end_[m] = end_[m_];

      if( !preserve ) {
         for( size_t i=0UL; i<m_; ++i )
            end_[i] = begin_[i];
      }

      for( size_t i=m_; i<m; ++i ) {
         begin_[i+1UL] = end_[i] = begin_[m_];
      }
   }
   else
   {
      if( preserve ) {
         for( size_t i=0UL; i<m; ++i )
            end_[i] = lowerBound( i, n );
      }
      else {
         for( size_t i=0UL; i<m; ++i )
            end_[i] = begin_[i];
      }

      end_[m] = end_[m_];
   }

   m_ = m;
   n_ = n;

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( size_t( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the compressed matrix.
//
// \param nonzeros The new minimum capacity of the compressed matrix.
// \return void
//
// This function increases the capacity of the compressed matrix to at least \a nonzeros elements.
// The current values of the matrix elements and the individual capacities of the matrix rows
// are preserved.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::reserve( size_t nonzeros )
{
   if( nonzeros > capacity() )
      reserveElements( nonzeros );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of a specific row/column of the compressed matrix.
//
// \param i The row/column index \f$[0..M-1]\f$ or \f$[0..N-1]\f$.
// \param nonzeros The new minimum capacity of the specified row/column.
// \return void
//
// This function increases the capacity of row/column \a i of the compressed matrix to at least
// \a nonzeros elements. The current values of the compressed matrix and all other individual
// row/column capacities are preserved. In case the storage order is set to \a rowMajor, the
// function reserves capacity for row \a i and the index has to be in the range \f$[0..M-1]\f$.
// In case the storage order is set to \a columnMajor, the function reserves capacity for column
// \a i and the index has to be in the range \f$[0..N-1]\f$.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
void CompressedMatrix<Type,SO,Tag>::reserve( size_t i, size_t nonzeros )
{
   using std::swap;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( static_cast<size_t>( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );

   const size_t current( capacity(i) );

   if( current >= nonzeros ) return;

   const ptrdiff_t additional( nonzeros - current );

   if( end_[m_] - begin_[m_] < additional )
   {
      const size_t newCapacity( begin_[m_] - begin_[0UL] + additional );
      BLAZE_INTERNAL_ASSERT( newCapacity > capacity(), "Invalid capacity value" );

      Iterator* newBegin( new Iterator[2UL*m_+2UL] );
      Iterator* newEnd  ( newBegin+m_+1UL );

      newBegin[0UL] = allocate<Element>( newCapacity );
      newEnd  [m_ ] = newBegin[0UL]+newCapacity;

      for( size_t k=0UL; k<i; ++k ) {
         newEnd  [k    ] = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
         newBegin[k+1UL] = newBegin[k] + capacity(k);
      }
      newEnd  [i    ] = castDown( transfer( begin_[i], end_[i], castUp( newBegin[i] ) ) );
      newBegin[i+1UL] = newBegin[i] + nonzeros;
      for( size_t k=i+1UL; k<m_; ++k ) {
         newEnd  [k    ] = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
         newBegin[k+1UL] = newBegin[k] + capacity(k);
      }

      BLAZE_INTERNAL_ASSERT( newBegin[m_] == newEnd[m_], "Invalid pointer calculations" );

      swap( newBegin, begin_ );
      deallocate( newBegin[0UL] );
      delete[] newBegin;
      end_ = newEnd;
      capacity_ = m_;
   }
   else
   {
      begin_[m_] += additional;
      for( size_t j=m_-1UL; j>i; --j ) {
         begin_[j]  = castDown( std::move_backward( begin_[j], end_[j], castUp( end_[j]+additional ) ) );
         end_  [j] += additional;
      }
   }

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( static_cast<size_t>( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity from all rows/columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all row/column-specific reserve()
// calls. The function removes all excessive capacity from all rows (in case of a rowMajor
// matrix) or columns (in case of a columnMajor matrix). Note that this function does not
// remove the overall capacity but only reduces the capacity per row/column.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::trim()
{
   for( size_t i=0UL; i<m_; ++i )
      trim( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Removing all excessive capacity of a specific row/column of the compressed matrix.
//
// \param i The index of the row/column to be trimmed (\f$[0..M-1]\f$ or \f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a row/column-specific reserve() call.
// It removes all excessive capacity from the specified row (in case of a rowMajor matrix)
// or column (in case of a columnMajor matrix). The excessive capacity is assigned to the
// subsequent row/column.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::trim( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   if( i < ( m_ - 1UL ) )
      end_[i+1] = castDown( std::move( begin_[i+1], end_[i+1], castUp( end_[i] ) ) );
   begin_[i+1] = end_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the matrix by removing unused capacity. Please note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this matrix are invalidated.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::shrinkToFit()
{
   if( nonZeros() < capacity() ) {
      CompressedMatrix( *this ).swap( *this );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two sparse matrices.
//
// \param sm The compressed matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::swap( CompressedMatrix& sm ) noexcept
{
   using std::swap;

   swap( m_, sm.m_ );
   swap( n_, sm.n_ );
   swap( capacity_, sm.capacity_ );
   swap( begin_, sm.begin_ );
   swap( end_  , sm.end_   );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating a new matrix capacity.
//
// \return The new compressed matrix capacity.
//
// This function calculates a new matrix capacity based on the current capacity of the sparse
// matrix. Note that the new capacity is restricted to the interval \f$[7..M \cdot N]\f$.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,SO,Tag>::extendCapacity() const noexcept
{
   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = blaze::max( nonzeros, 7UL   );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reserving the specified number of compressed matrix elements.
//
// \param nonzeros The number of matrix elements to be reserved.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
void CompressedMatrix<Type,SO,Tag>::reserveElements( size_t nonzeros )
{
   using std::swap;

   Iterator* newBegin = new Iterator[2UL*capacity_+2UL];
   Iterator* newEnd   = newBegin+capacity_+1UL;

   newBegin[0UL] = allocate<Element>( nonzeros );

   for( size_t k=0UL; k<m_; ++k ) {
      BLAZE_INTERNAL_ASSERT( begin_[k] <= end_[k], "Invalid row pointers" );
      newEnd  [k]     = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
      newBegin[k+1UL] = newBegin[k] + ( begin_[k+1UL] - begin_[k] );
   }

   newEnd[m_] = newBegin[0UL]+nonzeros;

   swap( newBegin, begin_ );
   end_ = newEnd;

   if( newBegin != nullptr ) {
      deallocate( newBegin[0UL] );
      delete[] newBegin;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs a down-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs a down-cast of the given iterator to base elements to an iterator to
// derived elements.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::castDown( IteratorBase it ) const noexcept
{
   return static_cast<Iterator>( it );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs an up-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs an up-cast of the given iterator to derived elements to an iterator
// to base elements.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::IteratorBase
   CompressedMatrix<Type,SO,Tag>::castUp( Iterator it ) const noexcept
{
   return static_cast<IteratorBase>( it );
}
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting an element of the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the compressed matrix. In case the compressed
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::set( size_t i, size_t j, const Type& value )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( lowerBound( i, j ) );

   if( pos != end_[i] && pos->index_ == j ) {
       pos->value() = value;
       return pos;
   }
   else return insert( pos, i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid compressed matrix access index.
//
// This function inserts a new element into the compressed matrix. However, duplicate elements
// are not allowed. In case the compressed matrix already contains an element with row index \a i
// and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::insert( size_t i, size_t j, const Type& value )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( lowerBound( i, j ) );

   if( pos != end_[i] && pos->index_ == j ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Bad access index" );
   }

   return insert( pos, i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the compressed matrix.
//
// \param pos The position of the new element.
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid compressed matrix access index.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::insert( Iterator pos, size_t i, size_t j, const Type& value )
{
   using std::swap;

   if( begin_[i+1UL] - end_[i] != 0 ) {
      std::move_backward( pos, end_[i], castUp( end_[i]+1UL ) );
      pos->value_ = value;
      pos->index_ = j;
      ++end_[i];

      return pos;
   }
   else if( end_[m_] - begin_[m_] != 0 ) {
      std::move_backward( pos, end_[m_-1UL], castUp( end_[m_-1UL]+1UL ) );

      pos->value_ = value;
      pos->index_ = j;

      for( size_t k=i+1UL; k<m_+1UL; ++k ) {
         ++begin_[k];
         ++end_[k-1UL];
      }

      return pos;
   }
   else {
      size_t newCapacity( extendCapacity() );

      Iterator* newBegin = new Iterator[2UL*capacity_+2UL];
      Iterator* newEnd   = newBegin+capacity_+1UL;

      newBegin[0UL] = allocate<Element>( newCapacity );

      for( size_t k=0UL; k<i; ++k ) {
         const size_t nonzeros( end_[k] - begin_[k] );
         const size_t total( begin_[k+1UL] - begin_[k] );
         newEnd  [k]     = newBegin[k] + nonzeros;
         newBegin[k+1UL] = newBegin[k] + total;
      }
      newEnd  [i]     = newBegin[i] + ( end_[i] - begin_[i] ) + 1;
      newBegin[i+1UL] = newBegin[i] + ( begin_[i+1] - begin_[i] ) + 1;
      for( size_t k=i+1UL; k<m_; ++k ) {
         const size_t nonzeros( end_[k] - begin_[k] );
         const size_t total( begin_[k+1UL] - begin_[k] );
         newEnd  [k]     = newBegin[k] + nonzeros;
         newBegin[k+1UL] = newBegin[k] + total;
      }

      newEnd[m_] = newEnd[capacity_] = newBegin[0UL]+newCapacity;

      Iterator tmp = castDown( std::move( begin_[0UL], pos, castUp( newBegin[0UL] ) ) );
      tmp->value_ = value;
      tmp->index_ = j;
      std::move( pos, end_[m_-1UL], castUp( tmp+1UL ) );

      swap( newBegin, begin_ );
      end_ = newEnd;
      deallocate( newBegin[0UL] );
      delete[] newBegin;

      return tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Appending an element to the specified row/column of the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a compressed matrix with elements. It
// appends a new element to the end of the specified row/column without any additional memory
// allocation. Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified row/column of the compressed matrix
//  - the current number of non-zero elements in the matrix must be smaller than the capacity
//    of the matrix
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a (newly created) compressed matrix:

   \code
   using blaze::rowMajor;

   // Setup of the compressed row-major matrix
   //
   //       ( 0 1 0 )
   //   A = ( 0 2 0 )
   //       ( 0 0 0 )
   //       ( 3 0 0 )
   //
   blaze::CompressedMatrix<double,rowMajor> A( 4, 3 );

   A.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   A.append( 0, 1, 1.0 );  // Appending the value 1 in row 0 with column index 1
   A.finalize( 0 );        // Finalizing row 0
   A.append( 1, 1, 2.0 );  // Appending the value 2 in row 1 with column index 1
   A.finalize( 1 );        // Finalizing row 1
   A.finalize( 2 );        // Finalizing the empty row 2 to prepare row 3
   A.append( 3, 0, 3.0 );  // Appending the value 3 in row 3 with column index 0
   A.finalize( 3 );        // Finalizing row 3
   \endcode

// \note The \c finalize() function has to be explicitly called for each row/column, even
// for empty ones!
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::append( size_t i, size_t j, const Type& value, bool check )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_USER_ASSERT( end_[i] < end_[m_], "Not enough reserved capacity left" );
   BLAZE_USER_ASSERT( begin_[i] == end_[i] || j > ( end_[i]-1UL )->index_, "Index is not strictly increasing" );

   end_[i]->value_ = value;

   if( !check || !isDefault<strict>( end_[i]->value_ ) ) {
      end_[i]->index_ = j;
      ++end_[i];
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Finalizing the element insertion of a row/column.
//
// \param i The index of the row/column to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill a matrix with elements.
// After completion of row/column \a i via the append() function, this function can be called to
// finalize row/column \a i and prepare the next row/column for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::finalize( size_t i )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid row access index" );

   begin_[i+1UL] = end_[i];
   if( i != m_-1UL )
      end_[i+1UL] = end_[i];
}
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Erasing an element from the compressed matrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,SO,Tag>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( find( i, j ) );
   if( pos != end_[i] )
      end_[i] = castDown( std::move( pos+1, end_[i], castUp( pos ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the compressed matrix.
//
// \param i The row/column index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the compressed matrix. In case the storage order is set
// to \a rowMajor the function erases an element from row \a i, in case the storage flag is set
// to \a columnMajor the function erases an element from column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::erase( size_t i, Iterator pos )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_USER_ASSERT( pos >= begin_[i] && pos <= end_[i], "Invalid compressed matrix iterator" );

   if( pos != end_[i] )
      end_[i] = castDown( std::move( pos+1, end_[i], castUp( pos ) ) );

   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the compressed matrix.
//
// \param i The row/column index of the elements to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the compressed matrix. In case the storage order
// is set to \a rowMajor the function erases a range of elements from row \a i, in case the storage
// flag is set to \a columnMajor the function erases a range of elements from column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::erase( size_t i, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index" );
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range"   );
   BLAZE_USER_ASSERT( first >= begin_[i] && first <= end_[i], "Invalid compressed matrix iterator" );
   BLAZE_USER_ASSERT( last  >= begin_[i] && last  <= end_[i], "Invalid compressed matrix iterator" );

   if( first != last )
      end_[i] = castDown( std::move( last, end_[i], castUp( first ) ) );

   return first;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing specific elements from the compressed matrix.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the compressed matrix. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   A.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the matrix
        , bool SO          // Storage order
        , typename Tag >   // Type tag
template< typename Pred >  // Type of the unary predicate
inline void CompressedMatrix<Type,SO,Tag>::erase( Pred predicate )
{
   for( size_t i=0UL; i<m_; ++i ) {
      end_[i] = castDown( std::remove_if( castUp( begin_[i] ), castUp( end_[i] ),
                                          [predicate=predicate]( const ElementBase& element) {
                                             return predicate( element.value() );
                                          } ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing specific elements from a range of the compressed matrix.
//
// \param i The row/column index of the elements to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements of the compressed matrix. The
// elements are selected by the given unary predicate \a predicate, which is expected to accept
// a single argument of the type of the elements and to be pure. In case the storage order is
// set to \a rowMajor the function erases a range of elements from row \a i, in case the storage
// flag is set to \a columnMajor the function erases a range of elements from column \a i. The
// following example demonstrates how to remove all elements that are smaller than a certain
// threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   A.erase( 2UL, A.begin(2UL), A.end(2UL), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the matrix
        , bool SO          // Storage order
        , typename Tag >   // Type tag
template< typename Pred >  // Type of the unary predicate
inline void CompressedMatrix<Type,SO,Tag>::erase( size_t i, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index" );
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range"   );
   BLAZE_USER_ASSERT( first >= begin_[i] && first <= end_[i], "Invalid compressed matrix iterator" );
   BLAZE_USER_ASSERT( last  >= begin_[i] && last  <= end_[i], "Invalid compressed matrix iterator" );

   const auto pos = std::remove_if( castUp( first ), castUp( last ),
                                    [predicate=predicate]( const ElementBase& element ) {
                                       return predicate( element.value() );
                                    } );

   end_[i] = castDown( std::move( last, end_[i], pos ) );
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
// end() iterator) is returned. Note that the returned compressed matrix iterator is subject
// to invalidation due to inserting operations via the function call operator, the set()
// function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::find( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).find( i, j ) );
}
//*************************************************************************************************


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
// end() iterator) is returned. Note that the returned compressed matrix iterator is subject
// to invalidation due to inserting operations via the function call operator, the set()
// function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::find( size_t i, size_t j ) const
{
   const ConstIterator pos( lowerBound( i, j ) );
   if( pos != end_[i] && pos->index_ == j )
      return pos;
   else return end_[i];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less than the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less than the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less than the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::lowerBound( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).lowerBound( i, j ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less than the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index not less than the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index not less than the given row
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::lowerBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return std::lower_bound( begin_[i], end_[i], j,
                            []( const Element& element, size_t index )
                            {
                               return element.index() < index;
                            } );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater than the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater than the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater than the given row
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator,
// the set() function or or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::Iterator
   CompressedMatrix<Type,SO,Tag>::upperBound( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).upperBound( i, j ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater than the given index, end() iterator otherwise.
//
// In case of a row-major matrix, this function returns a row iterator to the first element with
// an index greater than the given column index. In case of a column-major matrix, the function
// returns a column iterator to the first element with an index greater than the given row
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed matrix
// iterator is subject to invalidation due to inserting operations via the function call operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,SO,Tag>::ConstIterator
   CompressedMatrix<Type,SO,Tag>::upperBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return std::upper_bound( begin_[i], end_[i], j,
                            []( size_t index, const Element& element )
                            {
                               return index < element.index();
                            } );
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
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>& CompressedMatrix<Type,SO,Tag>::transpose()
{
   CompressedMatrix tmp( trans( *this ) );
   swap( tmp );
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
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,SO,Tag>& CompressedMatrix<Type,SO,Tag>::ctranspose()
{
   CompressedMatrix tmp( ctrans( *this ) );
   swap( tmp );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the compressed matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the compressed matrix.
//
// This function scales the matrix by applying the given scalar value \a scalar to each element
// of the matrix. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::CompressedMatrix<int> A;
   // ... Resizing and initialization
   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline CompressedMatrix<Type,SO,Tag>& CompressedMatrix<Type,SO,Tag>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<m_; ++i )
      for( auto element=begin_[i]; element!=end_[i]; ++element )
         element->value_ *= scalar;

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
inline bool CompressedMatrix<Type,SO,Tag>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
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
inline bool CompressedMatrix<Type,SO,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
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
        , typename Tag >  // Type tag
inline bool CompressedMatrix<Type,SO,Tag>::canSMPAssign() const noexcept
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,SO,Tag>::assign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   if( m_ == 0UL || n_ == 0UL )
      return;

   size_t nonzeros( 0UL );

   for( size_t i=1UL; i<=m_; ++i )
      begin_[i] = end_[i] = end_[m_];

   for( size_t i=0UL; i<m_; ++i )
   {
      begin_[i] = end_[i] = begin_[0UL]+nonzeros;

      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? i+1UL : i )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( IsStrictlyLower_v<MT> ? i : i+1UL )
                           :( n_ ) );

      for( size_t j=jbegin; j<jend; ++j )
      {
         if( nonzeros == capacity() ) {
            reserveElements( extendCapacity() );
            for( size_t k=i+1UL; k<=m_; ++k )
               begin_[k] = end_[k] = end_[m_];
         }

         end_[i]->value_ = (*rhs)(i,j);

         if( !isDefault<strict>( end_[i]->value_ ) ) {
            end_[i]->index_ = j;
            ++end_[i];
            ++nonzeros;
         }
      }
   }

   begin_[m_] = begin_[0UL]+nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a row-major compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side compressed matrix
inline void CompressedMatrix<Type,SO,Tag>::assign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );
   BLAZE_INTERNAL_ASSERT( capacity() >= (*rhs).nonZeros(), "Invalid capacity detected" );

   if( m_ == 0UL || begin_[0] == nullptr )
      return;

   for( size_t i=0UL; i<m_; ++i ) {
      end_[i] = castDown( std::copy( (*rhs).begin(i), (*rhs).end(i), castUp( begin_[i] ) ) );
      begin_[i+1UL] = end_[i];
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side compressed matrix
inline void CompressedMatrix<Type,SO,Tag>::assign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );
   BLAZE_INTERNAL_ASSERT( capacity() >= (*rhs).nonZeros(), "Invalid capacity detected" );

   // Counting the number of elements per row
   std::vector<size_t> rowLengths( m_, 0UL );
   for( size_t j=0UL; j<n_; ++j ) {
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         ++rowLengths[element->index()];
   }

   // Resizing the compressed matrix
   for( size_t i=0UL; i<m_; ++i ) {
      begin_[i+1UL] = end_[i+1UL] = begin_[i] + rowLengths[i];
   }

   // Appending the elements to the rows of the compressed matrix
   for( size_t j=0UL; j<n_; ++j ) {
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         append( element->index(), j, element->value() );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,SO,Tag>::addAssign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this + (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO2 >      // Storage order of the right-hand side compressed matrix
inline void CompressedMatrix<Type,SO,Tag>::addAssign( const SparseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this + (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,SO,Tag>::subAssign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this - (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO2 >      // Storage order of the right-hand compressed matrix
inline void CompressedMatrix<Type,SO,Tag>::subAssign( const SparseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this - (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the Schur product assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO2 >      // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,SO,Tag>::schurAssign( const DenseMatrix<MT,SO2>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   for( size_t i=0UL; i<m_; ++i ) {
      const Iterator last( end(i) );
      for( auto element=begin(i); element!=last; ++element )
         element->value_ *= (*rhs)(i,element->index_);
   }
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of CompressedMatrix for column-major matrices.
// \ingroup compressed_matrix
//
// This specialization of CompressedMatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
class CompressedMatrix<Type,true,Tag>
   : public SparseMatrix< CompressedMatrix<Type,true,Tag>, true >
{
 private:
   //**Type definitions****************************************************************************
   using ElementBase  = ValueIndexPair<Type>;  //!< Base class for the compressed matrix element.
   using IteratorBase = ElementBase*;          //!< Iterator over non-constant base elements.
   //**********************************************************************************************

   //**Private class Element***********************************************************************
   /*!\brief Value-index-pair for the CompressedMatrix class.
   //
   // This struct grants access to the data members of the base class and adapts the copy and
   // move semantics of the value-index-pair.
   */
   struct Element
      : public ElementBase
   {
      //**Constructors*****************************************************************************
      Element() = default;
      Element( const Element& rhs ) = default;
      Element( Element&& rhs ) = default;
      //*******************************************************************************************

      //**Assignment operators*********************************************************************
      inline Element& operator=( const Element& rhs )
      {
         this->value_ = rhs.value_;
         return *this;
      }

      inline Element& operator=( Element&& rhs )
      {
         this->value_ = std::move( rhs.value_ );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& rhs )
         -> EnableIf_t< IsSparseElement_v<Other>, Element& >
      {
         this->value_ = rhs.value();
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& rhs )
         -> EnableIf_t< IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( rhs.value() );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& v )
         -> EnableIf_t< !IsSparseElement_v<Other>, Element& >
      {
         this->value_ = v;
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& v )
         -> EnableIf_t< !IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( v );
         return *this;
      }
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      friend class CompressedMatrix;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Private class Uninitialized*****************************************************************
   /*!\brief Auxiliary helper class for the construction of compressed matrices.
   */
   struct Uninitialized {};
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This       = CompressedMatrix<Type,true,Tag>;  //!< Type of this CompressedMatrix instance.
   using BaseType   = SparseMatrix<This,true>;          //!< Base type of this CompressedMatrix instance.
   using ResultType = This;                             //!< Result type for expression template evaluations.

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = CompressedMatrix<Type,false,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = CompressedMatrix<Type,false,Tag>;

   using ElementType    = Type;                     //!< Type of the compressed matrix elements.
   using TagType        = Tag;                      //!< Tag type of this CompressedMatrix instance.
   using ReturnType     = const Type&;              //!< Return type for expression template evaluations.
   using CompositeType  = const This&;              //!< Data type for composite expression templates.
   using Reference      = MatrixAccessProxy<This>;  //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;              //!< Reference to a constant matrix value.
   using Iterator       = Element*;                 //!< Iterator over non-constant elements.
   using ConstIterator  = const Element*;           //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CompressedMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = CompressedMatrix<NewType,true,Tag>;  //!< The type of the other CompressedMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CompressedMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = CompressedMatrix<Type,true,Tag>;  //!< The type of the other CompressedMatrix.
   };
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
   inline CompressedMatrix();
   inline CompressedMatrix( size_t m, size_t n );
   inline CompressedMatrix( size_t m, size_t n, size_t nonzeros );
          CompressedMatrix( size_t m, size_t n, const std::vector<size_t>& nonzeros );
   inline CompressedMatrix( initializer_list< initializer_list<Type> > list );

   inline CompressedMatrix( const CompressedMatrix& sm );
   inline CompressedMatrix( CompressedMatrix&& sm ) noexcept;

   template< typename MT, bool SO > inline CompressedMatrix( const DenseMatrix<MT,SO>&  dm );
   template< typename MT, bool SO > inline CompressedMatrix( const SparseMatrix<MT,SO>& sm );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CompressedMatrix();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j ) noexcept;
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Iterator       begin ( size_t j ) noexcept;
   inline ConstIterator  begin ( size_t j ) const noexcept;
   inline ConstIterator  cbegin( size_t j ) const noexcept;
   inline Iterator       end   ( size_t j ) noexcept;
   inline ConstIterator  end   ( size_t j ) const noexcept;
   inline ConstIterator  cend  ( size_t j ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline CompressedMatrix& operator=( initializer_list< initializer_list<Type> > list ) &;
   inline CompressedMatrix& operator=( const CompressedMatrix& rhs ) &;
   inline CompressedMatrix& operator=( CompressedMatrix&& rhs ) & noexcept;

   template< typename MT, bool SO > inline CompressedMatrix& operator= ( const DenseMatrix<MT,SO>&  rhs ) &;
   template< typename MT, bool SO > inline CompressedMatrix& operator= ( const SparseMatrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline CompressedMatrix& operator+=( const Matrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline CompressedMatrix& operator-=( const Matrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline CompressedMatrix& operator%=( const DenseMatrix<MT,SO>&  rhs ) &;
   template< typename MT, bool SO > inline CompressedMatrix& operator%=( const SparseMatrix<MT,SO>& rhs ) &;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows() const noexcept;
   inline size_t columns() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   inline void   clear();
          void   resize ( size_t m, size_t n, bool preserve=true );
   inline void   reserve( size_t nonzeros );
          void   reserve( size_t j, size_t nonzeros );
   inline void   trim   ();
   inline void   trim   ( size_t j );
   inline void   shrinkToFit();
   inline void   swap( CompressedMatrix& sm ) noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set     ( size_t i, size_t j, const Type& value );
   inline Iterator insert  ( size_t i, size_t j, const Type& value );
   inline void     append  ( size_t i, size_t j, const Type& value, bool check=false );
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
   inline CompressedMatrix& transpose();
   inline CompressedMatrix& ctranspose();

   template< typename Other > inline CompressedMatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename MT, bool SO > inline void assign     ( const DenseMatrix<MT,SO>&     rhs );
   template< typename MT >          inline void assign     ( const SparseMatrix<MT,true>&  rhs );
   template< typename MT >          inline void assign     ( const SparseMatrix<MT,false>& rhs );
   template< typename MT, bool SO > inline void addAssign  ( const DenseMatrix<MT,SO>&     rhs );
   template< typename MT, bool SO > inline void addAssign  ( const SparseMatrix<MT,SO>&    rhs );
   template< typename MT, bool SO > inline void subAssign  ( const DenseMatrix<MT,SO>&     rhs );
   template< typename MT, bool SO > inline void subAssign  ( const SparseMatrix<MT,SO>&    rhs );
   template< typename MT, bool SO > inline void schurAssign( const DenseMatrix<MT,SO>&     rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline CompressedMatrix( size_t m, size_t n, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t extendCapacity() const noexcept;
          void   reserveElements( size_t nonzeros );

   inline Iterator     castDown( IteratorBase it ) const noexcept;
   inline IteratorBase castUp  ( Iterator     it ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   Iterator insert( Iterator pos, size_t i, size_t j, const Type& value );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;         //!< The current number of rows of the compressed matrix.
   size_t n_;         //!< The current number of columns of the compressed matrix.
   size_t capacity_;  //!< The current capacity of the pointer array.
   Iterator* begin_;  //!< Pointers to the first non-zero element of each column.
   Iterator* end_;    //!< Pointers one past the last non-zero element of each column.

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE       ( ElementBase, Element );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
const Type CompressedMatrix<Type,true,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The default constructor for CompressedMatrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix()
   : m_       ( 0UL )      // The current number of rows of the compressed matrix
   , n_       ( 0UL )      // The current number of columns of the compressed matrix
   , capacity_( 0UL )      // The current capacity of the pointer array
   , begin_   ( nullptr )  // Pointers to the first non-zero element of each column
   , end_     ( nullptr )  // Pointers one past the last non-zero element of each column
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
//
// The matrix is initialized to the zero matrix and has no free capacity.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( size_t m, size_t n )
   : CompressedMatrix( m, n, Uninitialized() )
{
   for( size_t j=1UL; j<2UL*n_+2UL; ++j )
      begin_[j] = nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nonzeros The number of expected non-zero elements.
//
// The matrix is initialized to the zero matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( size_t m, size_t n, size_t nonzeros )
   : CompressedMatrix( m, n, Uninitialized() )
{
   begin_[0UL] = allocate<Element>( nonzeros );
   for( size_t j=1UL; j<(2UL*n_+1UL); ++j )
      begin_[j] = begin_[0UL];
   end_[n_] = begin_[0UL]+nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nonzeros The expected number of non-zero elements in each column.
//
// The matrix is initialized to the zero matrix and will have the specified capacity in each
// column. Note that the given vector must have at least \a n elements.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
CompressedMatrix<Type,true,Tag>::CompressedMatrix( size_t m, size_t n, const std::vector<size_t>& nonzeros )
   : CompressedMatrix( m, n, Uninitialized() )
{
   BLAZE_USER_ASSERT( nonzeros.size() == n, "Size of capacity vector and number of columns don't match" );

   size_t newCapacity( 0UL );
   for( auto it=nonzeros.begin(); it!=nonzeros.end(); ++it )
      newCapacity += *it;

   begin_[0UL] = end_[0UL] = allocate<Element>( newCapacity );
   for( size_t j=0UL; j<n_; ++j ) {
      begin_[j+1UL] = end_[j+1UL] = begin_[j] + nonzeros[j];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List initialization of all matrix elements.
//
// \param list The initializer list.
//
// This constructor provides the option to explicitly initialize the elements of the matrix by
// means of an initializer list:

   \code
   using blaze::columnMajor;

   blaze::CompressedMatrix<int,columnMajor> A{ { 1, 2, 3 },
                                               { 4, 5 },
                                               { 7, 8, 9 } };
   \endcode

// The matrix is sized according to the size of the initializer list and all its elements are
// initialized by the values of the given initializer list. Missing values are considered to
// be default values.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( initializer_list< initializer_list<Type> > list )
   : CompressedMatrix( list.size(), determineColumns( list ), blaze::nonZeros( list ) )
{
   for( size_t j=0UL; j<n_; ++j )
   {
      size_t i( 0UL );

      for( const auto& rowList : list )
      {
         if( rowList.size() <= j ) {
            ++i;
            continue;
         }

         auto pos( rowList.begin() );
         std::advance( pos, j );
         if( !isDefault<strict>( *pos ) )
            append( i, j, *pos );
         ++i;
      }

      finalize( j );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for CompressedMatrix.
//
// \param sm Sparse matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( const CompressedMatrix& sm )
   : CompressedMatrix( sm.m_, sm.n_, Uninitialized() )
{
   const size_t nonzeros( sm.nonZeros() );

   begin_[0UL] = allocate<Element>( nonzeros );
   for( size_t j=0UL; j<n_; ++j ) {
      end_[j] = castDown( std::copy( sm.begin(j), sm.end(j), castUp( begin_[j] ) ) );
      begin_[j+1UL] = end_[j];
   }
   end_[n_] = begin_[0UL]+nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The move constructor for CompressedMatrix.
//
// \param sm The compressed matrix to be moved into this instance.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( CompressedMatrix&& sm ) noexcept
   : m_       ( sm.m_ )         // The current number of rows of the compressed matrix
   , n_       ( sm.n_ )         // The current number of columns of the compressed matrix
   , capacity_( sm.capacity_ )  // The current capacity of the pointer array
   , begin_   ( sm.begin_ )     // Pointers to the first non-zero element of each column
   , end_     ( sm.end_ )       // Pointers one past the last non-zero element of each column
{
   sm.m_        = 0UL;
   sm.n_        = 0UL;
   sm.capacity_ = 0UL;
   sm.begin_    = nullptr;
   sm.end_      = nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from dense matrices.
//
// \param dm Dense matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign dense matrix
        , bool SO >       // Storage order of the foreign dense matrix
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( const DenseMatrix<MT,SO>& dm )
   : CompressedMatrix( (*dm).rows(), (*dm).columns() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   assign( *this, *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from different sparse matrices.
//
// \param sm Sparse matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign compressed matrix
        , bool SO >       // Storage order of the foreign compressed matrix
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( const SparseMatrix<MT,SO>& sm )
   : CompressedMatrix( (*sm).rows(), (*sm).columns(), (*sm).nonZeros() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   assign( *this, *sm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for an uninitialized matrix of size \f$ M \times N \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::CompressedMatrix( size_t m, size_t n, Uninitialized )
   : m_       ( m )                     // The current number of rows of the compressed matrix
   , n_       ( n )                     // The current number of columns of the compressed matrix
   , capacity_( n )                     // The current capacity of the pointer array
   , begin_( new Iterator[2UL*n+2UL] )  // Pointers to the first non-zero element of each column
   , end_  ( begin_+(n+1UL) )           // Pointers one past the last non-zero element of each column
{
   begin_[0UL] = nullptr;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The destructor for CompressedMatrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>::~CompressedMatrix()
{
   if( begin_ != nullptr ) {
      deallocate( begin_[0UL] );
      delete[] begin_;
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
/*!\brief 2D-access to the compressed matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function returns a reference to the accessed value at position (\a i,\a j). In case the
// compressed matrix does not yet store an element at position (\a i,\a j) , a new element is
// inserted into the compressed matrix. Note that this function only performs an index check in
// case BLAZE_USER_ASSERT() is active. In contrast, the at() function is guaranteed to perform a
// check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Reference
   CompressedMatrix<Type,true,Tag>::operator()( size_t i, size_t j ) noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return Reference( *this, i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the compressed matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstReference
   CompressedMatrix<Type,true,Tag>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const ConstIterator pos( lowerBound( i, j ) );

   if( pos == end_[j] || pos->index_ != i )
      return zero_;
   else
      return pos->value_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// This function returns a reference to the accessed value at position (\a i,\a j). In case the
// compressed matrix does not yet store an element at position (\a i,\a j) , a new element is
// inserted into the compressed matrix. In contrast to the subscript operator this function
// always performs a check of the given access indices.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Reference
   CompressedMatrix<Type,true,Tag>::at( size_t i, size_t j )
{
   if( i >= m_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= n_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstReference
   CompressedMatrix<Type,true,Tag>::at( size_t i, size_t j ) const
{
   if( i >= m_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= n_ ) {
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
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::begin( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::begin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::cbegin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::end( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return end_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::end( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return end_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::cend( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid compressed matrix column access index" );
   return end_[j];
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
/*!\brief List assignment to all matrix elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::columnMajor;

   blaze::CompressedMatrix<int,columnMajor> A;
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix is resized according to the given initializer list and all its elements are
// assigned the values from the given initializer list. Missing values are considered to
// be default values.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator=( initializer_list< initializer_list<Type> > list ) &
{
   using blaze::nonZeros;

   resize( list.size(), determineColumns( list ), false );
   reserve( nonZeros( list ) );

   for( size_t j=0UL; j<n_; ++j )
   {
      size_t i( 0UL );

      for( const auto& rowList : list )
      {
         if( rowList.size() <= j ) {
            ++i;
            continue;
         }

         auto pos( rowList.begin() );
         std::advance( pos, j );
         if( !isDefault<strict>( *pos ) )
            append( i, j, *pos );
         ++i;
      }

      finalize( j );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for CompressedMatrix.
//
// \param rhs Sparse matrix to be copied.
// \return Reference to the assigned compressed matrix.
//
// The compressed matrix is resized according to the given compressed matrix and initialized
// as a copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator=( const CompressedMatrix& rhs ) &
{
   using std::swap;

   if( &rhs == this ) return *this;

   const size_t nonzeros( rhs.nonZeros() );

   if( rhs.n_ > capacity_ || nonzeros > capacity() )
   {
      Iterator* newBegin( new Iterator[2UL*rhs.n_+2UL] );
      Iterator* newEnd  ( newBegin+(rhs.n_+1UL) );

      newBegin[0UL] = allocate<Element>( nonzeros );
      for( size_t j=0UL; j<rhs.n_; ++j ) {
         newEnd[j] = castDown( std::copy( rhs.begin_[j], rhs.end_[j], castUp( newBegin[j] ) ) );
         newBegin[j+1UL] = newEnd[j];
      }
      newEnd[rhs.n_] = newBegin[0UL]+nonzeros;

      swap( begin_, newBegin );
      end_ = newEnd;
      capacity_ = rhs.n_;

      if( newBegin != nullptr ) {
         deallocate( newBegin[0UL] );
         delete[] newBegin;
      }
   }
   else {
      for( size_t j=0UL; j<rhs.n_; ++j ) {
         end_[j] = castDown( std::copy( rhs.begin_[j], rhs.end_[j], castUp( begin_[j] ) ) );
         begin_[j+1UL] = end_[j];
      }
   }

   m_ = rhs.m_;
   n_ = rhs.n_;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Move assignment operator for CompressedMatrix.
//
// \param rhs The compressed matrix to be moved into this instance.
// \return Reference to the assigned compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator=( CompressedMatrix&& rhs ) & noexcept
{
   if( begin_ != nullptr ) {
      deallocate( begin_[0UL] );
      delete[] begin_;
   }

   m_        = rhs.m_;
   n_        = rhs.n_;
   capacity_ = rhs.capacity_;
   begin_    = rhs.begin_;
   end_      = rhs.end_;

   rhs.m_        = 0UL;
   rhs.n_        = 0UL;
   rhs.capacity_ = 0UL;
   rhs.begin_    = nullptr;
   rhs.end_      = nullptr;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for dense matrices.
//
// \param rhs Dense matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator=( const DenseMatrix<MT,SO>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).canAlias( this ) ) {
      CompressedMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      assign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different sparse matrices.
//
// \param rhs Sparse matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO >       // Storage order of the right-hand side compressed matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator=( const SparseMatrix<MT,SO>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).canAlias( this ) ||
       (*rhs).columns()  > capacity_ ||
       (*rhs).nonZeros() > capacity() ) {
      CompressedMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      reset();

      if( !IsZero_v<MT> ) {
         assign( *this, *rhs );
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO >       // Storage order of the right-hand side matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator+=( const Matrix<MT,SO>& rhs ) &
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      addAssign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the matrix.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO >       // Storage order of the right-hand side matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator-=( const Matrix<MT,SO>& rhs ) &
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      subAssign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a dense matrix
//        (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator%=( const DenseMatrix<MT,SO>& rhs ) &
{
   using blaze::schurAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      CompressedMatrix tmp( *this % (*rhs) );
      swap( tmp );
   }
   else {
      CompositeType_t<MT> tmp( *rhs );
      schurAssign( *this, tmp );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a sparse matrix
//        (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side sparse matrix
        , bool SO >       // Storage order of the right-hand side sparse matrix
inline CompressedMatrix<Type,true,Tag>&
   CompressedMatrix<Type,true,Tag>::operator%=( const SparseMatrix<MT,SO>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsZero_v<MT> ) {
      CompressedMatrix tmp( *this % (*rhs) );
      swap( tmp );
   }
   else {
      reset();
   }

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
/*!\brief Returns the current number of rows of the compressed matrix.
//
// \return The number of rows of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::rows() const noexcept
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current number of columns of the compressed matrix.
//
// \return The number of columns of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::columns() const noexcept
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the compressed matrix.
//
// \return The capacity of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::capacity() const noexcept
{
   if( begin_ != nullptr )
      return end_[n_] - begin_[0UL];
   else return 0UL;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::capacity( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return begin_[j+1UL] - begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the compressed matrix
//
// \return The number of non-zero elements in the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<n_; ++j )
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
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return end_[j] - begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::reset()
{
   for( size_t j=0UL; j<n_; ++j )
      end_[j] = begin_[j];
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
// This function reset the values in the specified column to their default value. Note that
// the capacity of the column remains unchanged.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::reset( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   end_[j] = begin_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the compressed matrix.
//
// \return void
//
// After the clear() function, the size of the compressed matrix is 0.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::clear()
{
   if( end_ != nullptr )
      end_[0UL] = end_[n_];
   m_ = 0UL;
   n_ = 0UL;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Changing the size of the compressed matrix.
//
// \param m The new number of rows of the compressed matrix.
// \param n The new number of columns of the compressed matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Note that this function may invalidate all existing views (submatrices, rows, columns,
// ...) on the matrix if it is used to shrink the matrix. Additionally, the resize operation
// potentially changes all matrix elements. In order to preserve the old matrix values, the
// \a preserve flag can be set to \a true.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
void CompressedMatrix<Type,true,Tag>::resize( size_t m, size_t n, bool preserve )
{
   using std::swap;

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( begin_ == nullptr || size_t( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );

   if( m == m_ && n == n_ ) return;

   if( begin_ == nullptr )
   {
      begin_ = new Iterator[2UL*n+2UL];
      end_   = begin_+n+1UL;

      for( size_t j=0UL; j<2UL*n+2UL; ++j ) {
         begin_[j] = nullptr;
      }

      capacity_ = n;
   }
   else if( n > capacity_ )
   {
      Iterator* newBegin( new Iterator[2UL*n+2UL] );
      Iterator* newEnd  ( newBegin+n+1UL );

      newBegin[0UL] = begin_[0UL];

      if( preserve ) {
         for( size_t j=0UL; j<n_; ++j ) {
            newEnd  [j]     = end_  [j];
            newBegin[j+1UL] = begin_[j+1UL];
         }
         for( size_t j=n_; j<n; ++j ) {
            newBegin[j+1UL] = newEnd[j] = begin_[n_];
         }
      }
      else {
         for( size_t j=0UL; j<n; ++j ) {
            newBegin[j+1UL] = newEnd[j] = begin_[0UL];
         }
      }

      newEnd[n] = end_[n_];

      swap( newBegin, begin_ );
      delete[] newBegin;
      end_ = newEnd;
      capacity_ = n;
   }
   else if( n > n_ )
   {
      end_[n] = end_[n_];

      if( !preserve ) {
         for( size_t j=0UL; j<n_; ++j )
            end_[j] = begin_[j];
      }

      for( size_t j=n_; j<n; ++j ) {
         begin_[j+1UL] = end_[j] = begin_[n_];
      }
   }
   else
   {
      if( preserve ) {
         for( size_t j=0UL; j<n; ++j )
            end_[j] = lowerBound( m, j );
      }
      else {
         for( size_t j=0UL; j<n; ++j )
            end_[j] = begin_[j];
      }

      end_[n] = end_[n_];
   }

   m_ = m;
   n_ = n;

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( size_t( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the compressed matrix.
//
// \param nonzeros The new minimum capacity of the compressed matrix.
// \return void
//
// This function increases the capacity of the compressed matrix to at least \a nonzeros elements.
// The current values of the matrix elements and the individual capacities of the matrix rows
// are preserved.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::reserve( size_t nonzeros )
{
   if( nonzeros > capacity() )
      reserveElements( nonzeros );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of a specific column of the compressed matrix.
//
// \param j The column index. The index has to be in the range \f$[0..M-1]\f$.
// \param nonzeros The new minimum capacity of the specified column.
// \return void
//
// This function increases the capacity of column \a j of the compressed matrix to at least
// \a nonzeros elements. The current values of the compressed matrix and all other individual
// column capacities are preserved.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
void CompressedMatrix<Type,true,Tag>::reserve( size_t j, size_t nonzeros )
{
   using std::swap;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( static_cast<size_t>( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );

   const size_t current( capacity(j) );

   if( current >= nonzeros ) return;

   const ptrdiff_t additional( nonzeros - current );

   if( end_[n_] - begin_[n_] < additional )
   {
      const size_t newCapacity( begin_[n_] - begin_[0UL] + additional );
      BLAZE_INTERNAL_ASSERT( newCapacity > capacity(), "Invalid capacity value" );

      Iterator* newBegin( new Iterator[2UL*n_+2UL] );
      Iterator* newEnd  ( newBegin+n_+1UL );

      newBegin[0UL] = allocate<Element>( newCapacity );
      newEnd  [n_ ] = newBegin[0UL]+newCapacity;

      for( size_t k=0UL; k<j; ++k ) {
         newEnd  [k    ] = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
         newBegin[k+1UL] = newBegin[k] + capacity(k);
      }
      newEnd  [j    ] = castDown( transfer( begin_[j], end_[j], castUp( newBegin[j] ) ) );
      newBegin[j+1UL] = newBegin[j] + nonzeros;
      for( size_t k=j+1UL; k<n_; ++k ) {
         newEnd  [k    ] = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
         newBegin[k+1UL] = newBegin[k] + capacity(k);
      }

      BLAZE_INTERNAL_ASSERT( newBegin[n_] == newEnd[n_], "Invalid pointer calculations" );

      swap( newBegin, begin_ );
      deallocate( newBegin[0UL] );
      delete[] newBegin;
      end_ = newEnd;
      capacity_ = n_;
   }
   else
   {
      begin_[n_] += additional;
      for( size_t k=n_-1UL; k>j; --k ) {
         begin_[k]  = castDown( std::move_backward( begin_[k], end_[k], castUp( end_[k]+additional ) ) );
         end_  [k] += additional;
      }
   }

   BLAZE_INTERNAL_ASSERT( end_ >= begin_, "Invalid internal storage detected" );
   BLAZE_INTERNAL_ASSERT( static_cast<size_t>( end_ - begin_ ) == capacity_ + 1UL, "Invalid storage setting detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity from all columns.
//
// \return void
//
// The trim() function can be used to reverse the effect of all column-specific reserve() calls
// It removes all excessive capacity from all columns. Note that this function does not remove
// the overall capacity but only reduces the capacity per column.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
void CompressedMatrix<Type,true,Tag>::trim()
{
   for( size_t j=0UL; j<n_; ++j )
      trim( j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removing all excessive capacity of a specific column of the compressed matrix.
//
// \param j The index of the column to be trimmed (\f$[0..N-1]\f$).
// \return void
//
// This function can be used to reverse the effect of a column-specific reserve() call. It
// removes all excessive capacity from the specified column. The excessive capacity is assigned
// to the subsequent column.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
void CompressedMatrix<Type,true,Tag>::trim( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   if( j < ( n_ - 1UL ) )
      end_[j+1] = castDown( std::move( begin_[j+1], end_[j+1], castUp( end_[j] ) ) );
   begin_[j+1] = end_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the matrix by removing unused capacity. Please note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this matrix are invalidated.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::shrinkToFit()
{
   if( nonZeros() < capacity() ) {
      CompressedMatrix( *this ).swap( *this );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Swapping the contents of two sparse matrices.
//
// \param sm The compressed matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::swap( CompressedMatrix& sm ) noexcept
{
   using std::swap;

   swap( m_, sm.m_ );
   swap( n_, sm.n_ );
   swap( capacity_, sm.capacity_ );
   swap( begin_, sm.begin_ );
   swap( end_  , sm.end_   );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating a new matrix capacity.
//
// \return The new compressed matrix capacity.
//
// This function calculates a new matrix capacity based on the current capacity of the sparse
// matrix. Note that the new capacity is restricted to the interval \f$[7..M \cdot N]\f$.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t CompressedMatrix<Type,true,Tag>::extendCapacity() const noexcept
{
   size_t nonzeros( 2UL*capacity()+1UL );
   nonzeros = blaze::max( nonzeros, 7UL );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity(), "Invalid capacity value" );

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reserving the specified number of compressed matrix elements.
//
// \param nonzeros The number of matrix elements to be reserved.
// \return void
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
void CompressedMatrix<Type,true,Tag>::reserveElements( size_t nonzeros )
{
   using std::swap;

   Iterator* newBegin = new Iterator[2UL*capacity_+2UL];
   Iterator* newEnd   = newBegin+capacity_+1UL;

   newBegin[0UL] = allocate<Element>( nonzeros );

   for( size_t k=0UL; k<n_; ++k ) {
      BLAZE_INTERNAL_ASSERT( begin_[k] <= end_[k], "Invalid column pointers" );
      newEnd  [k]     = castDown( transfer( begin_[k], end_[k], castUp( newBegin[k] ) ) );
      newBegin[k+1UL] = newBegin[k] + ( begin_[k+1UL] - begin_[k] );
   }

   newEnd[n_] = newBegin[0UL]+nonzeros;

   swap( newBegin, begin_ );
   end_ = newEnd;

   if( newBegin != nullptr ) {
      deallocate( newBegin[0UL] );
      delete[] newBegin;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs a down-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs a down-cast of the given iterator to base elements to an iterator to
// derived elements.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::castDown( IteratorBase it ) const noexcept
{
   return static_cast<Iterator>( it );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs an up-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs an up-cast of the given iterator to derived elements to an iterator
// to base elements.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::IteratorBase
   CompressedMatrix<Type,true,Tag>::castUp( Iterator it ) const noexcept
{
   return static_cast<IteratorBase>( it );
}
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Iterator to the set element.
//
// This function sets the value of an element of the compressed matrix. In case the compressed
// matrix already contains an element with row index \a i and column index \a j its value is
// modified, else a new element with the given \a value is inserted.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::set( size_t i, size_t j, const Type& value )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( lowerBound( i, j ) );

   if( pos != end_[j] && pos->index_ == i ) {
      pos->value() = value;
      return pos;
   }
   else return insert( pos, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid compressed matrix access index.
//
// This function inserts a new element into the compressed matrix. However, duplicate elements
// are not allowed. In case the compressed matrix already contains an element with row index \a i
// and column index \a j, a \a std::invalid_argument exception is thrown.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::insert( size_t i, size_t j, const Type& value )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( lowerBound( i, j ) );

   if( pos != end_[j] && pos->index_ == i ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Bad access index" );
   }

   return insert( pos, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the compressed matrix.
//
// \param pos The position of the new element.
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Iterator to the newly inserted element.
// \exception std::invalid_argument Invalid compressed matrix access index.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::insert( Iterator pos, size_t i, size_t j, const Type& value )
{
   using std::swap;

   if( begin_[j+1UL] - end_[j] != 0 ) {
      std::move_backward( pos, end_[j], castUp( end_[j]+1UL ) );
      pos->value_ = value;
      pos->index_ = i;
      ++end_[j];

      return pos;
   }
   else if( end_[n_] - begin_[n_] != 0 ) {
      std::move_backward( pos, end_[n_-1UL], castUp( end_[n_-1]+1UL ) );

      pos->value_ = value;
      pos->index_ = i;

      for( size_t k=j+1UL; k<n_+1UL; ++k ) {
         ++begin_[k];
         ++end_[k-1UL];
      }

      return pos;
   }
   else {
      size_t newCapacity( extendCapacity() );

      Iterator* newBegin = new Iterator[2UL*capacity_+2UL];
      Iterator* newEnd   = newBegin+capacity_+1UL;

      newBegin[0UL] = allocate<Element>( newCapacity );

      for( size_t k=0UL; k<j; ++k ) {
         const size_t nonzeros( end_[k] - begin_[k] );
         const size_t total( begin_[k+1UL] - begin_[k] );
         newEnd  [k]     = newBegin[k] + nonzeros;
         newBegin[k+1UL] = newBegin[k] + total;
      }
      newEnd  [j]     = newBegin[j] + ( end_[j] - begin_[j] ) + 1;
      newBegin[j+1UL] = newBegin[j] + ( begin_[j+1UL] - begin_[j] ) + 1;
      for( size_t k=j+1UL; k<n_; ++k ) {
         const size_t nonzeros( end_[k] - begin_[k] );
         const size_t total( begin_[k+1UL] - begin_[k] );
         newEnd  [k]     = newBegin[k] + nonzeros;
         newBegin[k+1UL] = newBegin[k] + total;
      }

      newEnd[n_] = newEnd[capacity_] = newBegin[0UL]+newCapacity;

      Iterator tmp = castDown( std::move( begin_[0UL], pos, castUp( newBegin[0UL] ) ) );
      tmp->value_ = value;
      tmp->index_ = i;
      std::move( pos, end_[n_-1UL], castUp( tmp+1UL ) );

      swap( newBegin, begin_ );
      end_ = newEnd;
      deallocate( newBegin[0UL] );
      delete[] newBegin;

      return tmp;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the specified column of the compressed matrix.
//
// \param i The row index of the new element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a compressed matrix with elements. It
// appends a new element to the end of the specified column without any additional memory
// allocation. Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the specified column of the compressed matrix
//  - the current number of non-zero elements in the matrix must be smaller than the capacity of
//    the matrix.
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// In combination with the reserve() and the finalize() function, append() provides the most
// efficient way to add new elements to a (new created) compressed matrix:

   \code
   using blaze::columnMajor;

   // Setup of the compressed column-major matrix
   //
   //       ( 0 0 0 3 )
   //   A = ( 1 2 0 0 )
   //       ( 0 0 0 0 )
   //
   blaze::CompressedMatrix<double,columnMajor> A( 3, 4 );

   A.reserve( 3 );         // Reserving enough capacity for 3 non-zero elements
   A.append( 1, 0, 1.0 );  // Appending the value 1 in column 0 with row index 1
   A.finalize( 0 );        // Finalizing column 0
   A.append( 1, 1, 2.0 );  // Appending the value 2 in column 1 with row index 1
   A.finalize( 1 );        // Finalizing column 1
   A.finalize( 2 );        // Finalizing the empty column 2 to prepare column 3
   A.append( 0, 3, 3.0 );  // Appending the value 3 in column 3 with row index 0
   A.finalize( 3 );        // Finalizing column 3
   \endcode

// \note The \c finalize() function has to be explicitly called for each column, even for
// empty ones!
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::append( size_t i, size_t j, const Type& value, bool check )
{
   BLAZE_USER_ASSERT( i < m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_USER_ASSERT( end_[j] < end_[n_], "Not enough reserved capacity left" );
   BLAZE_USER_ASSERT( begin_[j] == end_[j] || i > ( end_[j]-1UL )->index_, "Index is not strictly increasing" );

   end_[j]->value_ = value;

   if( !check || !isDefault<strict>( end_[j]->value_ ) ) {
      end_[j]->index_ = i;
      ++end_[j];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Finalizing the element insertion of a column.
//
// \param j The index of the column to be finalized \f$[0..M-1]\f$.
// \return void
//
// This function is part of the low-level interface to efficiently fill the matrix with elements.
// After completion of column \a j via the append() function, this function can be called to
// finalize column \a j and prepare the next column for insertion process via append().
//
// \note Although finalize() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::finalize( size_t j )
{
   BLAZE_USER_ASSERT( j < n_, "Invalid column access index" );

   begin_[j+1UL] = end_[j];
   if( j != n_-1UL )
      end_[j+1UL] = end_[j];
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
/*!\brief Erasing an element from the compressed matrix.
//
// \param i The row index of the element to be erased. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void CompressedMatrix<Type,true,Tag>::erase( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const Iterator pos( find( i, j ) );
   if( pos != end_[j] )
      end_[j] = castDown( std::move( pos+1, end_[j], castUp( pos ) ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the compressed matrix.
//
// \param j The column index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from column \a j of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::erase( size_t j, Iterator pos )
{
   BLAZE_USER_ASSERT( j < columns()   , "Invalid column access index"    );
   BLAZE_USER_ASSERT( pos >= begin_[j] && pos <= end_[j], "Invalid compressed matrix iterator" );

   if( pos != end_[j] )
      end_[j] = castDown( std::move( pos+1, end_[j], castUp( pos ) ) );

   return pos;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the compressed matrix.
//
// \param j The column index of the elements to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return void
//
// This function erases a range of elements from column \a j of the compressed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::erase( size_t j, Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range"   );
   BLAZE_USER_ASSERT( first >= begin_[j] && first <= end_[j], "Invalid compressed matrix iterator" );
   BLAZE_USER_ASSERT( last  >= begin_[j] && last  <= end_[j], "Invalid compressed matrix iterator" );

   if( first != last )
      end_[j] = castDown( std::move( last, end_[j], castUp( first ) ) );

   return first;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the compressed matrix.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the compressed matrix. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   A.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the matrix
        , typename Tag >   // Type tag
template< typename Pred >  // Type of the unary predicate
inline void CompressedMatrix<Type,true,Tag>::erase( Pred predicate )
{
   for( size_t j=0UL; j<n_; ++j ) {
      end_[j] = castDown( std::remove_if( castUp( begin_[j] ), castUp( end_[j] ),
                                          [predicate=predicate]( const ElementBase& element ) {
                                             return predicate( element.value() );
                                          } ) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the compressed matrix.
//
// \param j The column index of the elements to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements of the compressed matrix.
// The elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements and to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedMatrix<double,blaze::columnMajor> A;
   // ... Resizing and initialization

   A.erase( 2UL, A.begin(2UL), A.end(2UL), []( double vaue ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the matrix
        , typename Tag >   // Type tag
template< typename Pred >  // Type of the unary predicate
inline void CompressedMatrix<Type,true,Tag>::erase( size_t j, Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range"   );
   BLAZE_USER_ASSERT( first >= begin_[j] && first <= end_[j], "Invalid compressed matrix iterator" );
   BLAZE_USER_ASSERT( last  >= begin_[j] && last  <= end_[j], "Invalid compressed matrix iterator" );

   const auto pos = std::remove_if( castUp( first ), castUp( last ),
                                    [predicate=predicate]( const ElementBase& element ) {
                                       return predicate( element.value() );
                                    } );

   end_[j] = castDown( std::move( last, end_[j], pos ) );
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
/*!\brief Searches for a specific matrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// matrix. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned compressed matrix iterator is subject to invalidation due to inserting
// operations via the subscript operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::find( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).find( i, j ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific matrix element.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// matrix. It specifically searches for the element with row index \a i and column index \a j.
// In case the element is found, the function returns an iterator to the element. Otherwise an
// iterator just past the last non-zero element of column \a j (the end() iterator) is returned.
// Note that the returned compressed matrix iterator is subject to invalidation due to inserting
// operations via the subscript operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::find( size_t i, size_t j ) const
{
   const ConstIterator pos( lowerBound( i, j ) );
   if( pos != end_[j] && pos->index_ == i )
      return pos;
   else return end_[j];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less than the given index, end() iterator otherwise.
//
// The function returns a column iterator to the first element with an index not less than the
// given row index. In combination with the upperBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// matrix iterator is subject to invalidation due to inserting operations via the function call
// operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::lowerBound( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).lowerBound( i, j ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less than the given index, end() iterator otherwise.
//
// The function returns a column iterator to the first element with an index not less than the
// given row index. In combination with the upperBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// matrix iterator is subject to invalidation due to inserting operations via the function call
// operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::lowerBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return std::lower_bound( begin_[j], end_[j], i,
                            []( const Element& element, size_t index )
                            {
                               return element.index() < index;
                            } );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater than the given index, end() iterator otherwise.
//
// The function returns a column iterator to the first element with an index greater than the
// given row index. In combination with the lowerBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// matrix iterator is subject to invalidation due to inserting operations via the function call
// operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::Iterator
   CompressedMatrix<Type,true,Tag>::upperBound( size_t i, size_t j )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).upperBound( i, j ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater than the given index.
//
// \param i The row index of the search element. The index has to be in the range \f$[0..M-1]\f$.
// \param j The column index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater than the given index, end() iterator otherwise.
//
// The function returns a column iterator to the first element with an index greater than the
// given row index. In combination with the lowerBound() function this function can be used to
// create a pair of iterators specifying a range of indices. Note that the returned compressed
// matrix iterator is subject to invalidation due to inserting operations via the function call
// operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename CompressedMatrix<Type,true,Tag>::ConstIterator
   CompressedMatrix<Type,true,Tag>::upperBound( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return std::upper_bound( begin_[j], end_[j], i,
                            []( size_t index, const Element& element )
                            {
                               return index < element.index();
                            } );
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
/*!\brief In-place transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>& CompressedMatrix<Type,true,Tag>::transpose()
{
   CompressedMatrix tmp( trans( *this ) );
   swap( tmp );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the matrix.
//
// \return Reference to the transposed matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline CompressedMatrix<Type,true,Tag>& CompressedMatrix<Type,true,Tag>::ctranspose()
{
   CompressedMatrix tmp( ctrans( *this ) );
   swap( tmp );
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the compressed matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the compressed matrix.
//
// This function scales the matrix by applying the given scalar value \a scalar to each element
// of the matrix. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::CompressedMatrix<int> A;
   // ... Resizing and initialization
   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline CompressedMatrix<Type,true,Tag>& CompressedMatrix<Type,true,Tag>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<n_; ++j )
      for( auto element=begin_[j]; element!=end_[j]; ++element )
         element->value_ *= scalar;

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
/*!\brief Returns whether the matrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool CompressedMatrix<Type,true,Tag>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this matrix, \a false if not.
//
// This function returns whether the given address is aliased with the vector. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool CompressedMatrix<Type,true,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
        , typename Tag >  // Type tag
inline bool CompressedMatrix<Type,true,Tag>::canSMPAssign() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,true,Tag>::assign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   if( m_ == 0UL || n_ == 0UL )
      return;

   size_t nonzeros( 0UL );

   for( size_t j=1UL; j<=n_; ++j )
      begin_[j] = end_[j] = end_[n_];

   for( size_t j=0UL; j<n_; ++j )
   {
      begin_[j] = end_[j] = begin_[0UL]+nonzeros;

      const size_t ibegin( ( IsLower_v<MT> )
                           ?( IsStrictlyLower_v<MT> ? j+1UL : j )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( m_ ) );

      for( size_t i=ibegin; i<iend; ++i )
      {
         if( nonzeros == capacity() ) {
            reserveElements( extendCapacity() );
            for( size_t k=j+1UL; k<=n_; ++k )
               begin_[k] = end_[k] = end_[n_];
         }

         end_[j]->value_ = (*rhs)(i,j);

         if( !isDefault<strict>( end_[j]->value_ ) ) {
            end_[j]->index_ = i;
            ++end_[j];
            ++nonzeros;
         }
      }
   }

   begin_[n_] = begin_[0UL]+nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side compressed matrix
inline void CompressedMatrix<Type,true,Tag>::assign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );
   BLAZE_INTERNAL_ASSERT( capacity() >= (*rhs).nonZeros(), "Invalid capacity detected" );

   if( n_ == 0UL || begin_[0] == nullptr )
      return;

   for( size_t j=0UL; j<n_; ++j ) {
      end_[j] = castDown( std::copy( (*rhs).begin(j), (*rhs).end(j), castUp( begin_[j] ) ) );
      begin_[j+1UL] = end_[j];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side compressed matrix
inline void CompressedMatrix<Type,true,Tag>::assign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );
   BLAZE_INTERNAL_ASSERT( capacity() >= (*rhs).nonZeros(), "Invalid capacity detected" );

   // Counting the number of elements per column
   std::vector<size_t> columnLengths( n_, 0UL );
   for( size_t i=0UL; i<m_; ++i ) {
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         ++columnLengths[element->index()];
   }

   // Resizing the compressed matrix
   for( size_t j=0UL; j<n_; ++j ) {
      begin_[j+1UL] = end_[j+1UL] = begin_[j] + columnLengths[j];
   }

   // Appending the elements to the columns of the compressed matrix
   for( size_t i=0UL; i<m_; ++i ) {
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         append( i, element->index(), element->value() );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,true,Tag>::addAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this + (*rhs) ) );
   swap( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO >       // Storage order of the right-hand side compressed matrix
inline void CompressedMatrix<Type,true,Tag>::addAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this + (*rhs) ) );
   swap( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,true,Tag>::subAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this - (*rhs) ) );
   swap( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a compressed matrix.
//
// \param rhs The right-hand side compressed matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side compressed matrix
        , bool SO >       // Storage order of the right-hand side compressed matrix
inline void CompressedMatrix<Type,true,Tag>::subAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   CompressedMatrix tmp( serial( *this - (*rhs) ) );
   swap( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side dense matrix
        , bool SO >       // Storage order of the right-hand side dense matrix
inline void CompressedMatrix<Type,true,Tag>::schurAssign( const DenseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   for( size_t j=0UL; j<n_; ++j ) {
      const Iterator last( end(j) );
      for( auto element=begin(j); element!=last; ++element )
         element->value_ *= (*rhs)(element->index_,j);
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  COMPRESSEDMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name CompressedMatrix operators */
//@{
template< typename Type, bool SO, typename Tag >
void reset( CompressedMatrix<Type,SO,Tag>& m );

template< typename Type, bool SO, typename Tag >
void reset( CompressedMatrix<Type,SO,Tag>& m, size_t i );

template< typename Type, bool SO, typename Tag >
void clear( CompressedMatrix<Type,SO,Tag>& m );

template< RelaxationFlag RF, typename Type, bool SO, typename Tag >
bool isDefault( const CompressedMatrix<Type,SO,Tag>& m );

template< typename Type, bool SO, typename Tag >
bool isIntact( const CompressedMatrix<Type,SO,Tag>& m );

template< typename Type, bool SO, typename Tag >
void swap( CompressedMatrix<Type,SO,Tag>& a, CompressedMatrix<Type,SO,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given compressed matrix.
// \ingroup compressed_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void reset( CompressedMatrix<Type,SO,Tag>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given compressed matrix.
// \ingroup compressed_matrix
//
// \param m The matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given compressed matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void reset( CompressedMatrix<Type,SO,Tag>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given compressed matrix.
// \ingroup compressed_matrix
//
// \param m The matrix to be cleared.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void clear( CompressedMatrix<Type,SO,Tag>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given compressed matrix is in default state.
// \ingroup compressed_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the compressed matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   blaze::CompressedMatrix<int> A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( A ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the matrix
        , bool SO            // Storage order
        , typename Tag >     // Type tag
inline bool isDefault( const CompressedMatrix<Type,SO,Tag>& m )
{
   return ( m.rows() == 0UL && m.columns() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given compressed matrix are intact.
// \ingroup compressed_matrix
//
// \param m The compressed matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the compressed matrix are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::CompressedMatrix<int> A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline bool isIntact( const CompressedMatrix<Type,SO,Tag>& m )
{
   return ( m.nonZeros() <= m.capacity() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two compressed matrices.
// \ingroup compressed_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline void swap( CompressedMatrix<Type,SO,Tag>& a, CompressedMatrix<Type,SO,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct IsResizable< CompressedMatrix<T,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool SO, typename Tag >
struct IsShrinkable< CompressedMatrix<T,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct AddTraitEval2< T1, T2
                    , EnableIf_t< IsSparseMatrix_v<T1> && IsSparseMatrix_v<T2> > >
{
   using Type = CompressedMatrix< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , ( StorageOrder_v<T1> && StorageOrder_v<T2> )
                                , AddTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SubTraitEval2< T1, T2
                    , EnableIf_t< IsSparseMatrix_v<T1> && IsSparseMatrix_v<T2> > >
{
   using Type = CompressedMatrix< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , ( StorageOrder_v<T1> && StorageOrder_v<T2> )
                                , SubTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
struct SchurTraitEval2< T1, T2
                      , EnableIf_t< IsMatrix_v<T1> &&
                                    IsMatrix_v<T2> &&
                                    ( IsSparseMatrix_v<T1> || IsSparseMatrix_v<T2> ) > >
{
   static constexpr bool SO = ( IsSparseMatrix_v<T1> && IsSparseMatrix_v<T2>
                                ? ( StorageOrder_v<T1> && StorageOrder_v<T2> )
                                : ( IsSparseMatrix_v<T1>
                                    ? StorageOrder_v<T1>
                                    : StorageOrder_v<T2> ) );

   using Type = CompressedMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
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
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsSparseMatrix_v<T1> && IsScalar_v<T2> > >
{
   using Type = CompressedMatrix< MultTrait_t< ElementType_t<T1>, T2 >
                                , StorageOrder_v<T1>
                                , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsScalar_v<T1> && IsSparseMatrix_v<T2> > >
{
   using Type = CompressedMatrix< MultTrait_t< T1, ElementType_t<T2> >
                                , StorageOrder_v<T2>
                                , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< ( IsSparseVector_v<T1> ||
                                     IsSparseVector_v<T2> ) &&
                                   IsColumnVector_v<T1> &&
                                   IsRowVector_v<T2> > >
{
   using Type = CompressedMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , ( IsSparseVector_v<T2> ? rowMajor : columnMajor )
                                , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsSparseMatrix_v<T1> &&
                                   IsSparseMatrix_v<T2> &&
                                   !( IsIdentity_v<T1> && IsIdentity_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = CompressedMatrix< AddTrait_t<MultType,MultType>
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
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
struct KronTraitEval2< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsSparseMatrix_v<T1> || IsSparseMatrix_v<T2> ) > >
{
   using Type = CompressedMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , ( IsDenseMatrix_v<T2> ? StorageOrder_v<T1> : StorageOrder_v<T2> )
                                , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsSparseMatrix_v<T1> && IsScalar_v<T2> > >
{
   using Type = CompressedMatrix< DivTrait_t< ElementType_t<T1>, T2 >
                                , StorageOrder_v<T1>
                                , DivTrait_t< TagType_t<T1>, T2 > >;
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
struct UnaryMapTraitEval2< T, OP
                         , EnableIf_t< IsSparseMatrix_v<T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = CompressedMatrix< EvaluateTrait_t<ElementType>
                                , StorageOrder_v<T>
                                , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T  // Type to be expanded
        , size_t E >  // Compile time expansion
struct ExpandTraitEval2< T, E
                       , EnableIf_t< IsSparseVector_v<T> > >
{
   using Type = CompressedMatrix< ElementType_t<T>
                                , ( IsColumnVector_v<T> ? columnMajor : rowMajor )
                                , TagType_t<T> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REPEATTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, size_t R0, size_t R1 >
struct RepeatTraitEval2< T, R0, R1, inf
                       , EnableIf_t< IsSparseMatrix_v<T> > >
{
   using Type = CompressedMatrix< ElementType_t<T>
                                , StorageOrder_v<T>
                                , TagType_t<T> >;
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
struct HighType< CompressedMatrix<T1,SO,Tag>, CompressedMatrix<T2,SO,Tag> >
{
   using Type = CompressedMatrix< typename HighType<T1,T2>::Type, SO, Tag >;
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
struct LowType< CompressedMatrix<T1,SO,Tag>, CompressedMatrix<T2,SO,Tag> >
{
   using Type = CompressedMatrix< typename LowType<T1,T2>::Type, SO, Tag >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I, size_t J, size_t M, size_t N >
struct SubmatrixTraitEval2< MT, I, J, M, N
                          , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedMatrix< RemoveConst_t< ElementType_t<MT> >
                                , StorageOrder_v<MT>
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t M >
struct RowsTraitEval2< MT, M
                     , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedMatrix< RemoveConst_t< ElementType_t<MT> >
                                , false
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t N >
struct ColumnsTraitEval2< MT, N
                        , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedMatrix< RemoveConst_t< ElementType_t<MT> >
                                , true
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
