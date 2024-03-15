//=================================================================================================
/*!
//  \file blaze/math/dense/InitializerMatrix.h
//  \brief Header file for the implementation of a matrix representation of an initializer list
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

#ifndef _BLAZE_MATH_DENSE_INITIALIZERMATRIX_H_
#define _BLAZE_MATH_DENSE_INITIALIZERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/dense/InitializerIterator.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsInitializer.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup initializer_matrix InitializerMatrix
// \ingroup dense_matrix
*/
/*!\brief Dense matrix representation of an initializer list.
// \ingroup initializer_matrix
//
// The InitializerMatrix class template is a dense matrix representation of an (extended)
// initializer list of arbitrary type. The type of the elements and the group tag of the matrix
// can be specified via the two template parameters:

   \code
   namespace blaze {

   template< typename Type, typename Tag >
   class InitializerMatrix;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the matrix elements. InitializerMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// On construction, an InitializerMatrix is immediately bound to an initializer list:

   \code
   const blaze::initializer_list< initializer_list<int> > list = { { 2, 6, -1 },
                                                                   { 3, 5 } };

   blaze::InitializerMatrix<int> A( list );  // Representation of the initializer list as dense matrix
   \endcode

// It is possible to only represent an extended initializer list by explicitly specifying the
// number of columns:

   \code
   const initializer_list< initializer_list<int> > list = { { 2, 6, -1 },
                                                            { 3, 5 } };

   blaze::InitializerVector<int> B( list, 3UL );  // Representation of the original initializer list
   blaze::InitializerVector<int> C( list, 4UL );  // Representing the initializer list { { 2, 6, -1, 0 }, { 3, 5, 0, 0 } }
   \endcode

// Since an InitializerMatrix represents a specific initializer list, its lifetime is bound to
// the lifetime of the according initializer list. When the initializer list goes out of scope
// access to the initializer list via an InitializerMatrix results in undefined behavior:

   \code
   blaze::InitializerMatrix<int> D{ { 1, 2, 3 }, { 4, 5, 6 } };         // Undefined behavior!
   blaze::InitializerMatrix<int> E( { { 0, 3, 2 }, { -1, 1 } }, 3UL );  // Undefined behavior!
   \endcode

// Also, an InitializerMatrix can only be used on the right of an assignment as its elements
// are considered to be immutable. The following example gives an impression on the usage of an
// InitializerMatrix:

   \code
   const blaze::initializer_list< initializer_list<int> > list = { { 2, 6, -1 },
                                                                   { 3, 5 } };

   blaze::InitializerMatrix<int> F( list );  // Representation of the initializer list as dense matrix
   blaze::DynamicMatrix<int> G;

   G = F;  // Initialize vector G via vector F
   F = G;  // Compilation error! Cannot assign to an initializer matrix
   \endcode

// An initializer matrix can be used as operand in arithmetic operations. All operations (addition,
// subtraction, multiplication, scaling, ...) can be performed on all possible combinations of
// row-major and column-major dense and sparse matrices with fitting element types. The following
// example gives an impression of the use of InitializerMatrix:

   \code
   using blaze::initializer_list;
   using blaze::InitializerMatrix;
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   const blaze::initializer_list< initializer_list<double> > list = { { 1.0, 2.0, 3.0 },
                                                                      { 4.0, 5.0, 6.0 } };

   InitializerMatrix<double> A( list );

   DynamicMatrix<float,columnMajor> B( 2, 3 );  // Default constructed column-major single precision 2x3 matrix
   B(0,0) = 1.0; B(0,1) = 3.0; B(0,2) = 5.0;    // Initialization of the first row
   B(1,0) = 2.0; B(1,1) = 4.0; B(1,2) = 6.0;    // Initialization of the second row

   CompressedMatrix<float> C( 2, 3 );        // Empty row-major sparse single precision matrix
   DynamicMatrix<float>    D( 3, 2, 4.0F );  // Directly, homogeneously initialized single precision 3x2 matrix

   DynamicMatrix<double,rowMajor>    E( A );  // Creation of a new row-major matrix as a copy of A
   DynamicMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;    // Matrix addition and assignment to a row-major matrix
   F = A - C;    // Matrix subtraction and assignment to a column-major matrix
   F = A * D;    // Matrix multiplication between two matrices of different element types

   E = 2.0 * B;  // Scaling of matrix B
   F = D * 2.0;  // Scaling of matrix D

   E += A - B;   // Addition assignment
   E -= A + C;   // Subtraction assignment
   F *= A * D;   // Multiplication assignment
   \endcode
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
class InitializerMatrix
   : public DenseMatrix< InitializerMatrix<Type,Tag>, false >
{
 public:
   //**Type definitions****************************************************************************
   using This     = InitializerMatrix<Type,Tag>;  //!< Type of this InitializerMatrix instance.
   using BaseType = DenseMatrix<This,false>;      //!< Base type of this InitializerMatrix instance.

   //! Result type for expression template evaluations.
   using ResultType = DynamicMatrix<Type,false>;

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType  = DynamicMatrix<Type,true>;

   using TransposeType = DynamicMatrix<Type,true>;  //!< Transpose type for expression template evaluations.
   using ElementType   = Type;                      //!< Type of the matrix elements.
   using TagType       = Tag;                       //!< Tag type of this InitializerVector instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = const Type&;  //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = const Type*;  //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = InitializerIterator<Type>;  //!< Iterator over non-constant elements.
   using ConstIterator = InitializerIterator<Type>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a InitializerMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = InitializerMatrix<NewType,Tag>;  //!< The type of the other InitializerMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a InitializerMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = InitializerMatrix<Type,Tag>;  //!< The type of the other InitializerMatrix.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the matrix is involved
       in can be optimized via SIMD operations. In case the element type of the matrix is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = false;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline InitializerMatrix( initializer_list< initializer_list<Type> > list ) noexcept;
   inline InitializerMatrix( initializer_list< initializer_list<Type> > list, size_t n );

   InitializerMatrix( const InitializerMatrix& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~InitializerMatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline ConstReference at( size_t i, size_t j ) const;
   inline ConstPointer   data  () const noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline ConstIterator  begin ( size_t i ) const noexcept;
   inline ConstIterator  cbegin( size_t i ) const noexcept;
   inline ConstIterator  end   ( size_t i ) const noexcept;
   inline ConstIterator  cend  ( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   InitializerMatrix& operator=( const InitializerMatrix& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows() const noexcept;
   inline size_t columns() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   swap( InitializerMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using ListType = initializer_list< initializer_list<Type> >;  //!< Type of the represented initializer list.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;       //!< The current number of rows of the matrix.
   size_t n_;       //!< The current number of columns of the matrix.
   ListType list_;  //!< The initializer list represented by the matrix.
                    /*!< Access to the matrix elements is gained via the function call
                         operator. The memory layout of the elements is
                         \f[\left(\begin{array}{*{5}{c}}
                         0            & 1             & 2             & \cdots & N-1         \\
                         N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                         \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                         M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                         \end{array}\right)\f]. */

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
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
        , typename Tag >  // Type tag
const Type InitializerMatrix<Type,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for InitializerMatrix.
//
// \param list The initializer list represented by the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline InitializerMatrix<Type,Tag>::InitializerMatrix( initializer_list< initializer_list<Type> > list ) noexcept
   : m_   ( list.size() )               // The current number of rows of the matrix
   , n_   ( determineColumns( list ) )  // The current number of columns of the matrix
   , list_( list )                      // The initializer list represented by the matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for InitializerMatrix.
//
// \param list The initializer list represented by the matrix.
// \param n The number of columns of the matrix.
// \exception std::invalid_argument Invalid initializer list dimension.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline InitializerMatrix<Type,Tag>::InitializerMatrix( initializer_list< initializer_list<Type> > list, size_t n )
   : m_   ( list.size() )  // The current number of rows of the matrix
   , n_   ( n    )         // The current number of columns of the matrix
   , list_( list )         // The initializer list represented by the matrix
{
   if( n < determineColumns( list ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid initializer list dimension" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief 2D-access to the matrix elements.
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
inline typename InitializerMatrix<Type,Tag>::ConstReference
   InitializerMatrix<Type,Tag>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );

   const initializer_list<Type>& list( list_.begin()[i] );

   if( j < list.size() )
      return list.begin()[j];
   else
      return zero_;
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
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstReference
   InitializerMatrix<Type,Tag>::at( size_t i, size_t j ) const
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
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data. Whereas the number of
// elements within a row/column are given by the \c rows() and \c columns() member functions,
// respectively, the total number of elements including padding is given by the \c spacing()
// member function.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstPointer
   InitializerMatrix<Type,Tag>::data() const noexcept
{
   return list_.begin()->begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row \a i.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstPointer
   InitializerMatrix<Type,Tag>::data( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return list_.begin()[i].begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row \a i.
//
// \param i The row index.
// \return Iterator to the first element of row \a i.
//
// This function returns a row iterator to the first element of row \a i.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstIterator
   InitializerMatrix<Type,Tag>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( 0UL, list_.begin()[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row \a i.
//
// \param i The row index.
// \return Iterator to the first element of row \a i.
//
// This function returns a row iterator to the first element of row \a i
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstIterator
   InitializerMatrix<Type,Tag>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( 0UL, list_.begin()[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last element of row \a i.
//
// This function returns an row iterator just past the last element of row \a i.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstIterator
   InitializerMatrix<Type,Tag>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( n_, list_.begin()[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row \a i.
//
// \param i The row index.
// \return Iterator just past the last element of row \a i.
//
// This function returns an row iterator just past the last element of row \a i.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline typename InitializerMatrix<Type,Tag>::ConstIterator
   InitializerMatrix<Type,Tag>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( n_, list_.begin()[i] );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::rows() const noexcept
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows.
//
// \return The spacing between the beginning of two rows.
//
// This function returns the spacing between the beginning of two rows, i.e. the total number
// of elements of a row.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::spacing() const noexcept
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::capacity() const noexcept
{
   return m_ * n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current capacity of the specified row.
//
// \param i The index of the row.
// \return The current capacity of row \a i.
//
// This function returns the current capacity of the specified row.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the dense matrix.
//
// This function returns the number of non-zero elements in the matrix (i.e. the elements that
// compare unequal to their default value). Note that the number of non-zero elements is always
// less than or equal to the total number of elements in the matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( const auto& rowList : list_ ) {
      for( size_t i=0UL; i<rowList.size(); ++i ) {
         if( !isDefault<strict>( rowList.begin()[i] ) )
            ++nonzeros;
      }
   }

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row.
//
// \param i The index of the row.
// \return The number of non-zero elements of row \a i.
//
// This function returns the current number of non-zero elements in the specified row (i.e. the
// elements that compare unequal to their default value).
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline size_t InitializerMatrix<Type,Tag>::nonZeros( size_t i ) const
{
   using blaze::nonZeros;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return nonZeros( list_.begin()[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void InitializerMatrix<Type,Tag>::swap( InitializerMatrix& m ) noexcept
{
   using std::swap;

   swap( m_   , m.m_    );
   swap( n_   , m.n_    );
   swap( list_, m.list_ );
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
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool InitializerMatrix<Type,Tag>::canAlias( const Other* alias ) const noexcept
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
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool InitializerMatrix<Type,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************




//=================================================================================================
//
//  INITIALIZERMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name InitializerMatrix operators */
//@{
template< typename Type, typename Tag >
bool isIntact( const InitializerMatrix<Type,Tag>& m ) noexcept;

template< typename Type, typename Tag >
void swap( InitializerMatrix<Type,Tag>& a, InitializerMatrix<Type,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given initializer matrix are intact.
// \ingroup initializer_matrix
//
// \param m The initializer matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the initializer matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::InitializerMatrix<int> A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline bool isIntact( const InitializerMatrix<Type,Tag>& m ) noexcept
{
   MAYBE_UNUSED( m );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two initializer matrices.
// \ingroup initializer_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , typename Tag >  // Type tag
inline void swap( InitializerMatrix<Type,Tag>& a, InitializerMatrix<Type,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename Tag >
struct HasConstDataAccess< InitializerMatrix<T,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISINITIALIZER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename Tag >
struct IsInitializer< InitializerMatrix<T,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename Tag, typename T2 >
struct HighType< InitializerMatrix<T1,Tag>, InitializerMatrix<T2,Tag> >
{
   using Type = InitializerMatrix< typename HighType<T1,T2>::Type, Tag >;
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
template< typename T1, typename Tag, typename T2 >
struct LowType< InitializerMatrix<T1,Tag>, InitializerMatrix<T2,Tag> >
{
   using Type = InitializerMatrix< typename LowType<T1,T2>::Type, Tag >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
