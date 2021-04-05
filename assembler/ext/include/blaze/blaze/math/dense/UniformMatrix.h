//=================================================================================================
/*!
//  \file blaze/math/dense/UniformMatrix.h
//  \brief Header file for the implementation of a uniform matrix
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

#ifndef _BLAZE_MATH_DENSE_UNIFORMMATRIX_H_
#define _BLAZE_MATH_DENSE_UNIFORMMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/dense/UniformIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/Forward.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
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
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/typetraits/RemoveConst.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup uniform_matrix UniformMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient implementation of a uniform matrix.
// \ingroup uniform_matrix
//
// The UniformMatrix class template is the representation of an arbitrary sized uniform matrix
// with elements of arbitrary type. The type of the elements, the storage order, and the group
// tag of the matrix can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Tag >
   class UniformMatrix;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the matrix elements. UniformMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the matrix.
//          The default value is \c blaze::defaultStorageOrder.
//  - Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// Depending on the storage order, the matrix elements are either stored in a row-wise fashion
// or in a column-wise fashion. Given the 2x3 matrix

                          \f[\left(\begin{array}{*{3}{c}}
                          1 & 2 & 3 \\
                          4 & 5 & 6 \\
                          \end{array}\right)\f]\n

// in case of row-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 2 & 3 & 4 & 5 & 6. \\
                          \end{array}\right)\f]

// In case of column-major order the elements are stored in the order

                          \f[\left(\begin{array}{*{6}{c}}
                          1 & 4 & 2 & 5 & 3 & 6. \\
                          \end{array}\right)\f]

// The use of UniformMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of row-major and
// column-major dense and sparse matrices with fitting element types. The following example gives
// an impression of the use of UniformMatrix:

   \code
   using blaze::UniformMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   UniformMatrix<double,rowMajor> A( 2, 3 );  // Default initialized, row-major 2x3 uniform matrix
   A = 1.0;                                   // Assignment to all elements of the uniform matrix

   UniformMatrix<float,columnMajor> B( 2, 3, 2.0F );  // Directly, uniformly initialized 2x3 matrix
   CompressedMatrix<float> C( 2, 3 );                 // Empty row-major sparse single precision matrix

   UniformMatrix<double,rowMajor> D;  // Creation of a new row-major matrix as a copy of A

   D = A + B;     // Matrix addition and assignment to a row-major matrix
   D = A - C;     // Matrix subtraction and assignment to a column-major matrix
   D = A * B;     // Matrix multiplication between two matrices of different element types

   A *= 2.0;      // In-place scaling of matrix
   D  = 2.0 * A;  // Scaling of matrix A
   D  = A * 2.0;  // Scaling of matrix A

   D += A - B;    // Addition assignment
   D -= A + C;    // Subtraction assignment
   D *= A * B;    // Multiplication assignment
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
class UniformMatrix
   : public Expression< DenseMatrix< UniformMatrix<Type,SO,Tag>, SO > >
{
 public:
   //**Type definitions****************************************************************************
   using This     = UniformMatrix<Type,SO,Tag>;          //!< Type of this UniformMatrix instance.
   using BaseType = Expression< DenseMatrix<This,SO> >;  //!< Base type of this UniformMatrix instance.

   //! Result type for expression template evaluations.
   using ResultType = This;

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = UniformMatrix<Type,!SO,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = UniformMatrix<Type,!SO,Tag>;

   using ElementType   = Type;                      //!< Type of the matrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using TagType       = Tag;                       //!< Tag type of this UniformMatrix instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = const Type&;  //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = const Type*;  //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = UniformIterator<const Type,aligned>;  //!< Iterator over non-constant elements.
   using ConstIterator = UniformIterator<const Type,aligned>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a UniformMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using Other = UniformMatrix<NewType,SO,Tag>;  //!< The type of the other UniformMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a UniformMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using Other = UniformMatrix<Type,SO,Tag>;  //!< The type of the other UniformMatrix.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the matrix is involved
       in can be optimized via SIMD operations. In case the element type of the matrix is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = IsVectorizable_v<Type>;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the matrix can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr UniformMatrix() noexcept;
   constexpr UniformMatrix( size_t m, size_t n );
   constexpr UniformMatrix( size_t m, size_t n, const Type& init );

   template< typename MT, bool SO2 >
   inline UniformMatrix( const Matrix<MT,SO2>& m );

   inline UniformMatrix( const UniformMatrix& m ) = default;
   inline UniformMatrix( UniformMatrix&& m ) = default;

   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~UniformMatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline    ConstReference at( size_t i, size_t j ) const;
   constexpr ConstPointer   data  () const noexcept;
   constexpr ConstPointer   data  ( size_t i ) const noexcept;
   constexpr ConstIterator  begin ( size_t i ) const noexcept;
   constexpr ConstIterator  cbegin( size_t i ) const noexcept;
   constexpr ConstIterator  end   ( size_t i ) const noexcept;
   constexpr ConstIterator  cend  ( size_t i ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   constexpr UniformMatrix& operator=( const Type& rhs ) &;

   UniformMatrix& operator=( const UniformMatrix& ) & = default;
   UniformMatrix& operator=( UniformMatrix&& ) & = default;

   template< typename MT, bool SO2 > inline UniformMatrix& operator= ( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline UniformMatrix& operator+=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline UniformMatrix& operator-=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline UniformMatrix& operator%=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline UniformMatrix& operator*=( const Matrix<MT,SO2>& rhs ) &;

   template< typename ST >
   inline auto operator*=( ST rhs ) & -> EnableIf_t< IsScalar_v<ST>, UniformMatrix& >;

   template< typename ST >
   inline auto operator/=( ST rhs ) & -> EnableIf_t< IsScalar_v<ST>, UniformMatrix& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   constexpr size_t rows() const noexcept;
   constexpr size_t columns() const noexcept;
   constexpr size_t spacing() const noexcept;
   constexpr size_t capacity() const noexcept;
   constexpr size_t capacity( size_t i ) const noexcept;
   inline    size_t nonZeros() const;
   inline    size_t nonZeros( size_t i ) const;
   constexpr void   reset();
   constexpr void   clear();
   constexpr void   resize ( size_t m, size_t n, bool preserve=true );
   constexpr void   extend ( size_t m, size_t n, bool preserve=true );
   constexpr void   swap( UniformMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   constexpr UniformMatrix& transpose();
   constexpr UniformMatrix& ctranspose();

   template< typename Other > inline UniformMatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;      //!< The current number of rows of the matrix.
   size_t n_;      //!< The current number of columns of the matrix.
   Type   value_;  //!< The value of all elements of the uniform matrix.
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
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for UniformMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr UniformMatrix<Type,SO,Tag>::UniformMatrix() noexcept
   : m_    ( 0UL )  // The current number of rows of the matrix
   , n_    ( 0UL )  // The current number of columns of the matrix
   , value_()       // The value of all elements of the uniform matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr UniformMatrix<Type,SO,Tag>::UniformMatrix( size_t m, size_t n )
   : m_    ( m )  // The current number of rows of the matrix
   , n_    ( n )  // The current number of columns of the matrix
   , value_()     // The value of all elements of the uniform matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
constexpr UniformMatrix<Type,SO,Tag>::UniformMatrix( size_t m, size_t n, const Type& init )
   : m_    ( m )     // The current number of rows of the matrix
   , n_    ( n )     // The current number of columns of the matrix
   , value_( init )  // The value of all elements of the uniform matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
// \exception std::invalid_argument Invalid setup of uniform vector.
//
// The matrix is sized according to the given uniform matrix and initialized as a copy of this
// matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Tag type
template< typename MT     // Type of the foreign matrix
        , bool SO2 >      // Storage order of the foreign matrix
inline UniformMatrix<Type,SO,Tag>::UniformMatrix( const Matrix<MT,SO2>& m )
   : m_    ( (*m).rows()    )  // The current number of rows of the matrix
   , n_    ( (*m).columns() )  // The current number of columns of the matrix
   , value_()                  // The value of all elements of the uniform vector
{
   if( !IsUniform_v<MT> && !isUniform( *m ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of uniform matrix" );
   }

   if( m_ > 0UL && n_ > 0UL ) {
      value_ = (*m)(0UL,0UL);
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
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstReference
   UniformMatrix<Type,SO,Tag>::operator()( size_t i, size_t j ) const noexcept
{
   MAYBE_UNUSED( i, j );

   BLAZE_USER_ASSERT( i < m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < n_, "Invalid column access index" );

   return value_;
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
        , typename Tag >  // Type tag
inline typename UniformMatrix<Type,SO,Tag>::ConstReference
   UniformMatrix<Type,SO,Tag>::at( size_t i, size_t j ) const
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
// This function returns a pointer to the internal storage of the uniform matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The uniform matrix may
// use techniques such as padding to improve the alignment of the data. Whereas the number of
// elements within a row/column are given by the \c rows() and \c columns() member functions,
// respectively, the total number of elements including padding is given by the \c spacing()
// member function.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstPointer
   UniformMatrix<Type,SO,Tag>::data() const noexcept
{
   return &value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the matrix elements of row/column \a i.
//
// \param i The row/column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row/column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstPointer
   UniformMatrix<Type,SO,Tag>::data( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   return &value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// This function returns a row/column iterator to the first element of row/column \a i. In case
// the storage order is set to \a rowMajor the function returns an iterator to the first element
// of row \a i, in case the storage flag is set to \a columnMajor the function returns an iterator
// to the first element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstIterator
   UniformMatrix<Type,SO,Tag>::begin( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   return ConstIterator( &value_, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first element of row/column \a i.
//
// This function returns a row/column iterator to the first element of row/column \a i. In case
// the storage order is set to \a rowMajor the function returns an iterator to the first element
// of row \a i, in case the storage flag is set to \a columnMajor the function returns an iterator
// to the first element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstIterator
   UniformMatrix<Type,SO,Tag>::cbegin( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   return ConstIterator( &value_, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// This function returns an row/column iterator just past the last element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator just past
// the last element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator just past the last element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstIterator
   UniformMatrix<Type,SO,Tag>::end( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   return ConstIterator( &value_, SO ? m_ : n_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last element of row/column \a i.
//
// This function returns an row/column iterator just past the last element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator just past
// the last element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator just past the last element of column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr typename UniformMatrix<Type,SO,Tag>::ConstIterator
   UniformMatrix<Type,SO,Tag>::cend( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   return ConstIterator( &value_, SO ? m_ : n_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all matrix elements.
//
// \param rhs Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator=( const Type& rhs ) &
{
   value_ = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator=( const Matrix<MT,SO2>& rhs ) &
{
   using TT = decltype( trans( *this ) );
   using CT = decltype( ctrans( *this ) );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( !IsUniform_v<MT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment of uniform matrix" );
   }

   if( IsSame_v<MT,TT> && (*rhs).isAliased( this ) ) {
      transpose();
   }
   else if( IsSame_v<MT,CT> && (*rhs).isAliased( this ) ) {
      ctranspose();
   }
   else {
      m_ = (*rhs).rows();
      n_ = (*rhs).columns();

      if( (*rhs).rows() > 0UL && (*rhs).columns() > 0UL ) {
         value_ = (*rhs)(0UL,0UL);
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
inline UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator+=( const Matrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsUniform_v<MT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid addition assignment to uniform matrix" );
   }

   if( m_ > 0UL && n_ > 0UL ) {
      value_ += (*rhs)(0UL,0UL);
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
inline UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator-=( const Matrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsUniform_v<MT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subtraction assignment to uniform matrix" );
   }

   if( m_ > 0UL && n_ > 0UL ) {
      value_ -= (*rhs)(0UL,0UL);
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
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
inline UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator%=( const Matrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsUniform_v<MT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid Schur product assignment to uniform matrix" );
   }

   if( m_ > 0UL && n_ > 0UL ) {
      value_ *= (*rhs)(0UL,0UL);
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a matrix (\f$ A*=B \f$).
//
// \param rhs The right-hand side matrix for the multiplication.
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
inline UniformMatrix<Type,SO,Tag>&
   UniformMatrix<Type,SO,Tag>::operator*=( const Matrix<MT,SO2>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !IsUniform_v<MT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid multiplication assignment to uniform matrix" );
   }

   n_ = (*rhs).columns();

   if( m_ > 0UL && n_ > 0UL ) {
      value_ = ( value_ * (*rhs)(0UL,0UL) ) * Type( (*rhs).rows() );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a matrix and
//        a scalar value (\f$ A*=s \f$).
//
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename ST >   // Data type of the right-hand side scalar
inline auto UniformMatrix<Type,SO,Tag>::operator*=( ST scalar ) &
   -> EnableIf_t< IsScalar_v<ST>, UniformMatrix& >
{
   if( rows() > 0UL && columns() > 0UL ) {
      value_ *= scalar;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division between a matrix and a scalar value
//        (\f$ A/=s \f$).
//
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
template< typename ST >   // Data type of the right-hand side scalar
inline auto UniformMatrix<Type,SO,Tag>::operator/=( ST scalar ) &
   -> EnableIf_t< IsScalar_v<ST>, UniformMatrix& >
{
   if( rows() > 0UL && columns() > 0UL ) {
      value_ /= scalar;
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
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr size_t UniformMatrix<Type,SO,Tag>::rows() const noexcept
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
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr size_t UniformMatrix<Type,SO,Tag>::columns() const noexcept
{
   return n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr size_t UniformMatrix<Type,SO,Tag>::spacing() const noexcept
{
   return SO ? m_ : n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr size_t UniformMatrix<Type,SO,Tag>::capacity() const noexcept
{
   return m_ * n_;
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
constexpr size_t UniformMatrix<Type,SO,Tag>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );
   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );
   return SO ? m_ : n_;
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
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t UniformMatrix<Type,SO,Tag>::nonZeros() const
{
   if( m_ == 0UL || n_ == 0UL || isDefault<strict>( value_ ) )
      return 0UL;
   else
      return m_ * n_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column
// (i.e. the elements that compare unequal to their default value). In case the storage order
// is set to \a rowMajor the function returns the number of non-zero elements in row \a i, in
// case the storage flag is set to \a columnMajor the function returns the number of non-zero
// elements in column \a i.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline size_t UniformMatrix<Type,SO,Tag>::nonZeros( size_t i ) const
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( SO  || i < m_, "Invalid dense matrix row access index" );
   BLAZE_USER_ASSERT( !SO || i < n_, "Invalid dense matrix row access index" );

   const size_t tmp( SO ? m_ : n_ );
   if( tmp == 0UL || isDefault<strict>( value_ ) )
      return 0UL;
   else
      return tmp;
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
constexpr void UniformMatrix<Type,SO,Tag>::reset()
{
   using blaze::clear;

   clear( value_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the \f$ M \times N \f$ matrix.
//
// \return void
//
// After the clear() function, the size of the matrix is 0.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void UniformMatrix<Type,SO,Tag>::clear()
{
   m_ = 0UL;
   n_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. Note that this
// function may invalidate all existing views (submatrices, rows, columns, ...) on the matrix
// if it is used to shrink the matrix. Additionally, the resize operation potentially changes
// all matrix elements. In order to preserve the old matrix values, the \a preserve flag can
// be set to \a true.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
void constexpr UniformMatrix<Type,SO,Tag>::resize( size_t m, size_t n, bool preserve )
{
   MAYBE_UNUSED( preserve );

   m_  = m;
   n_  = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the matrix.
//
// \param m Number of additional rows.
// \param n Number of additional columns.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function increases the matrix size by \a m rows and \a n columns. Note that this function
// potentially changes all matrix elements. In order to preserve the old matrix values, the
// \a preserve flag can be set to \a true.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void UniformMatrix<Type,SO,Tag>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void UniformMatrix<Type,SO,Tag>::swap( UniformMatrix& m ) noexcept
{
   using std::swap;

   swap( m_, m.m_ );
   swap( n_, m.n_ );
   swap( value_, m.value_ );
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
constexpr UniformMatrix<Type,SO,Tag>& UniformMatrix<Type,SO,Tag>::transpose()
{
   using std::swap;

   swap( m_, n_ );

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
constexpr UniformMatrix<Type,SO,Tag>& UniformMatrix<Type,SO,Tag>::ctranspose()
{
   using std::swap;

   swap( m_, n_ );
   conjugate( value_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
//
// This function scales the matrix by applying the given scalar value \a scalar to each element
// of the matrix. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::UniformMatrix<int> A;
   // ... Resizing and initialization
   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline UniformMatrix<Type,SO,Tag>& UniformMatrix<Type,SO,Tag>::scale( const Other& scalar )
{
   if( m_ > 0UL && n_ > 0UL ) {
      value_ *= scalar;
   }

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
inline bool UniformMatrix<Type,SO,Tag>::canAlias( const Other* alias ) const noexcept
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
inline bool UniformMatrix<Type,SO,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix is properly aligned in memory.
//
// \return \a true in case the matrix is aligned, \a false if not.
//
// This function returns whether the matrix is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of each row/column of the matrix are guaranteed to conform
// to the alignment restrictions of the element type \a Type.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
inline bool UniformMatrix<Type,SO,Tag>::isAligned() const noexcept
{
   return true;
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
inline bool UniformMatrix<Type,SO,Tag>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense matrix. The row index
// must be smaller than the number of rows and the column index must be smaller then the number
// of columns. Additionally, the column index (in case of a row-major matrix) or the row index
// (in case of a column-major matrix) must be a multiple of the number of values inside the
// SIMD element. This function must \b NOT be called explicitly! It is used internally for the
// performance optimized evaluation of expression templates. Calling this function explicitly
// might result in erroneous results and/or in compilation errors.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformMatrix<Type,SO,Tag>::SIMDType
   UniformMatrix<Type,SO,Tag>::load( size_t i, size_t j ) const noexcept
{
   return loada( i, j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformMatrix<Type,SO,Tag>::SIMDType
   UniformMatrix<Type,SO,Tag>::loada( size_t i, size_t j ) const noexcept
{
   MAYBE_UNUSED( i, j );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( SO  || j + SIMDSIZE <= n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !SO || i + SIMDSIZE <= m_, "Invalid row access index" );

   return set( value_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformMatrix<Type,SO,Tag>::SIMDType
   UniformMatrix<Type,SO,Tag>::loadu( size_t i, size_t j ) const noexcept
{
   MAYBE_UNUSED( i, j );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( SO  || j + SIMDSIZE <= n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !SO || i + SIMDSIZE <= m_, "Invalid row access index" );

   return set( value_ );
}
//*************************************************************************************************








//=================================================================================================
//
//  UNIFORMMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name UniformMatrix operators */
//@{
template< typename Type, bool SO, typename Tag >
constexpr void reset( UniformMatrix<Type,SO,Tag>& m );

template< typename Type, bool SO, typename Tag >
constexpr void clear( UniformMatrix<Type,SO,Tag>& m );

template< RelaxationFlag RF, typename Type, bool SO, typename Tag >
constexpr bool isDefault( const UniformMatrix<Type,SO,Tag>& m );

template< typename Type, bool SO, typename Tag >
constexpr bool isIntact( const UniformMatrix<Type,SO,Tag>& m ) noexcept;

template< typename Type, bool SO, typename Tag >
constexpr void swap( UniformMatrix<Type,SO,Tag>& a, UniformMatrix<Type,SO,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given uniform matrix.
// \ingroup uniform_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void reset( UniformMatrix<Type,SO,Tag>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given dynamic matrix.
// \ingroup dynamic_matrix
//
// \param m The matrix to be cleared.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void clear( UniformMatrix<Type,SO,Tag>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given uniform matrix is in default state.
// \ingroup uniform_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the uniform matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   blaze::UniformMatrix<int> A;
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
constexpr bool isDefault( const UniformMatrix<Type,SO,Tag>& m )
{
   return ( m.rows() == 0UL && m.columns() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given uniform matrix are intact.
// \ingroup uniform_matrix
//
// \param m The uniform matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the uniform matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::UniformMatrix<int> A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr bool isIntact( const UniformMatrix<Type,SO,Tag>& m ) noexcept
{
   MAYBE_UNUSED( m );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two uniform matrices.
// \ingroup uniform_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
constexpr void swap( UniformMatrix<Type,SO,Tag>& a, UniformMatrix<Type,SO,Tag>& b ) noexcept
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
/*!\brief Creating a uniform matrix.
// \ingroup uniform_matrix
//
// \param m The number of rows of the uniform matrix.
// \param n The number of columns of the uniform matrix.
// \param init The initial value of the matrix elements.
// \return A uniform matrix of the given size.
//
// This function creates a uniform matrix of the given size. By default, the resulting uniform
// matrix is a row-major matrix, but it is possible to specify the storage order explicitly:

   \code
   using blaze::uniform;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Creates the uniform row-major matrix
   //    ( 1, 1, 1, 1, 1 )
   //    ( 1, 1, 1, 1, 1 )
   auto U1 = uniform<int>( 2UL, 5UL, 1 );

   // Creates the uniform row-major matrix
   //    ( 1.2, 1.2 )
   //    ( 1.2, 1.2 )
   //    ( 1.2, 1.2 )
   auto U2 = uniform<int,rowMajor>( 3UL, 2UL, 1.2 );

   // Creates the uniform column-major matrix
   //   ( 5U, 5U, 5U, 5U, 5U, 5U, 5U )
   //   ( 5U, 5U, 5U, 5U, 5U, 5U, 5U )
   auto U3 = uniform<int,columnMajor>( 2UL, 7UL, 5U );
   \endcode
*/
template< bool SO = defaultStorageOrder, typename T >
constexpr decltype(auto) uniform( size_t m, size_t n, T&& init )
{
   return UniformMatrix< RemoveCVRef_t<T>, SO >( m, n, std::forward<T>( init ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename Type, bool SO, typename Tag >
struct IsUniform< UniformMatrix<Type,SO,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename Type, bool SO, typename Tag >
struct IsAligned< UniformMatrix<Type,SO,Tag> >
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
template< typename Type, bool SO, typename Tag >
struct IsResizable< UniformMatrix<Type,SO,Tag> >
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
struct AddTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( IsDenseMatrix_v<T1> && IsDenseMatrix_v<T2>
                                ? ( IsSymmetric_v<T1> ^ IsSymmetric_v<T2>
                                    ? ( IsSymmetric_v<T1>
                                        ? SO2
                                        : SO1 )
                                    : SO1 && SO2 )
                                : ( IsDenseMatrix_v<T1>
                                    ? SO1
                                    : SO2 ) );

   using Type = UniformMatrix< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , SO
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
struct SubTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( IsDenseMatrix_v<T1> && IsDenseMatrix_v<T2>
                                ? ( IsSymmetric_v<T1> ^ IsSymmetric_v<T2>
                                    ? ( IsSymmetric_v<T1>
                                        ? SO2
                                        : SO1 )
                                    : SO1 && SO2 )
                                : ( IsDenseMatrix_v<T1>
                                    ? SO1
                                    : SO2 ) );

   using Type = UniformMatrix< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , SO
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
template< typename T1, typename T2 >
struct SchurTraitEval1< T1, T2
                      , EnableIf_t< IsMatrix_v<T1> &&
                                    IsMatrix_v<T2> &&
                                    ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                    !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( IsSymmetric_v<T1> ^ IsSymmetric_v<T2>
                                ? ( IsSymmetric_v<T1>
                                    ? SO2
                                    : SO1 )
                                : SO1 && SO2 );

   using Type = UniformMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
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
                                   IsUniform_v<T1> &&
                                   !IsZero_v<T1> &&
                                   IsScalar_v<T2> > >
{
   using Type = UniformMatrix< MultTrait_t< ElementType_t<T1>, T2 >
                             , StorageOrder_v<T1>
                             , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   IsUniform_v<T2> &&
                                   !IsZero_v<T2> > >
{
   using Type = UniformMatrix< MultTrait_t< T1, ElementType_t<T2> >
                             , StorageOrder_v<T2>
                             , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsColumnVector_v<T1> &&
                                   IsRowVector_v<T2> &&
                                   ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , false
                             , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = UniformMatrix< AddTrait_t<MultType,MultType>
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
                                   ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformMatrix< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , StorageOrder_v<T2>
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
struct DivTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsScalar_v<T2> &&
                                  IsUniform_v<T1> && !IsZero_v<T1> > >
{
   using Type = UniformMatrix< DivTrait_t< ElementType_t<T1>, T2 >
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
struct UnaryMapTraitEval1< T, OP
                         , EnableIf_t< IsMatrix_v<T> &&
                                       YieldsUniform_v<OP,T> &&
                                       !YieldsZero_v<OP,T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = UniformMatrix< EvaluateTrait_t<ElementType>
                             , StorageOrder_v<T>
                             , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< IsColumnVector_v<T1> &&
                                        IsRowVector_v<T2> &&
                                        YieldsUniform_v<OP,T1,T2> &&
                                        !YieldsZero_v<OP,T1,T2> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   using Type = UniformMatrix< EvaluateTrait_t<ElementType>
                             , false
                             , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
};

template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< IsMatrix_v<T1> &&
                                        IsMatrix_v<T2> &&
                                        YieldsUniform_v<OP,T1,T2> &&
                                        !YieldsZero_v<OP,T1,T2> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( IsDenseMatrix_v<T1> && IsDenseMatrix_v<T2>
                                ? ( IsSymmetric_v<T1> ^ IsSymmetric_v<T2>
                                    ? ( IsSymmetric_v<T1>
                                        ? SO2
                                        : SO1 )
                                    : SO1 && SO2 )
                                : ( IsDenseMatrix_v<T1>
                                    ? SO1
                                    : SO2 ) );

   using Type = UniformMatrix< EvaluateTrait_t<ElementType>
                             , SO
                             , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
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
template< typename T, size_t E >
struct ExpandTraitEval1< T, E
                       , EnableIf_t< IsVector_v<T> &&
                                     IsUniform_v<T> && !IsZero_v<T> > >
{
   using Type = UniformMatrix< ElementType_t<T>
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
struct RepeatTraitEval1< T, R0, R1, inf
                       , EnableIf_t< IsMatrix_v<T> &&
                                     IsUniform_v<T> && !IsZero_v<T> > >
{
   using Type = UniformMatrix< ElementType_t<T>
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
struct HighType< UniformMatrix<T1,SO,Tag>, UniformMatrix<T2,SO,Tag> >
{
   using Type = UniformMatrix< typename HighType<T1,T2>::Type, SO, Tag >;
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
struct LowType< UniformMatrix<T1,SO,Tag>, UniformMatrix<T2,SO,Tag> >
{
   using Type = UniformMatrix< typename LowType<T1,T2>::Type, SO, Tag >;
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
struct SubmatrixTraitEval1< MT, I, J, M, N
                          , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformMatrix< RemoveConst_t< ElementType_t<MT> >
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
struct RowsTraitEval1< MT, M
                     , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformMatrix< RemoveConst_t< ElementType_t<MT> >
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
struct ColumnsTraitEval1< MT, N
                        , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformMatrix< RemoveConst_t< ElementType_t<MT> >
                             , true
                             , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
