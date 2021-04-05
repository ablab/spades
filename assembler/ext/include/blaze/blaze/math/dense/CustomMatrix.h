//=================================================================================================
/*!
//  \file blaze/math/dense/CustomMatrix.h
//  \brief Header file for the implementation of a customizable matrix
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

#ifndef _BLAZE_MATH_DENSE_CUSTOMMATRIX_H_
#define _BLAZE_MATH_DENSE_CUSTOMMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <algorithm>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Diagonal.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/functors/Noop.h>
#include <blaze/math/functors/Clear.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/PaddingFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/CustomOppositeType.h>
#include <blaze/math/typetraits/CustomTransposeType.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsCustom.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Misalignment.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/typetraits/RemoveConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup custom_matrix CustomMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient implementation of a customizable matrix.
// \ingroup custom_matrix
//
// The CustomMatrix class template provides the functionality to represent an external array of
// elements of arbitrary type and a fixed size as a native \b Blaze dense matrix data structure.
// Thus in contrast to all other dense matrix types a custom matrix does not perform any kind
// of memory allocation by itself, but it is provided with an existing array of element during
// construction. A custom matrix can therefore be considered an alias to the existing array.
//
// The type of the elements, the properties of the given array of elements, the storage order,
// and the group tag of the matrix can be specified via the following five template parameters:

   \code
   namespace blaze {

   template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag >
   class CustomMatrix;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the matrix elements. CustomMatrix can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - AF  : specifies whether the represented, external arrays are properly aligned with
//          respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::aligned
//          or \c blaze::unaligned).
//  - PF  : specified whether the represented, external arrays are properly padded with
//          respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::padded
//          or \c blaze::unpadded).
//  - SO  : specifies the storage order (\c blaze::rowMajor, \c blaze::columnMajor) of the
//          matrix. The default value is \c blaze::defaultStorageOrder.
//  - Tag : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// The following examples give an impression of several possible types of custom matrices:

   \code
   using blaze::CustomMatrix;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a custom row-major matrix for unaligned, unpadded integer arrays
   using UnalignedUnpadded = CustomMatrix<int,unaligned,unpadded,rowMajor>;

   // Definition of a custom column-major matrix for unaligned but padded 'float' arrays
   using UnalignedPadded = CustomMatrix<float,unaligned,padded,columnMajor>;

   // Definition of a custom row-major matrix for aligned, unpadded 'double' arrays
   using AlignedUnpadded = CustomMatrix<double,aligned,unpadded,rowMajor>;

   // Definition of a custom column-major matrix for aligned, padded 'complex<double>' arrays
   using AlignedPadded = CustomMatrix<complex<double>,aligned,padded,columnMajor>;
   \endcode

// \n \section custommatrix_special_properties Special Properties of Custom Matrices
//
// In comparison with the remaining \b Blaze dense matrix types CustomMatrix has several special
// characteristics. All of these result from the fact that a custom matrix is not performing any
// kind of memory allocation, but instead is given an existing array of elements. The following
// sections discuss all of these characteristics:
//
//  -# <b>\ref custommatrix_memory_management</b>
//  -# <b>\ref custommatrix_copy_operations</b>
//  -# <b>\ref custommatrix_alignment</b>
//  -# <b>\ref custommatrix_padding</b>
//
// \n \subsection custommatrix_memory_management Memory Management
//
// The CustomMatrix class template acts as an adaptor for an existing array of elements. As such
// it provides everything that is required to use the array just like a native \b Blaze dense
// matrix data structure. However, this flexibility comes with the price that the user of a custom
// matrix is responsible for the resource management.
//
// The following examples give an impression of several possible custom matrices:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a 3x4 custom row-major matrix with unaligned, unpadded and externally
   // managed integer array. Note that the std::vector must be guaranteed to outlive the
   // custom matrix!
   std::vector<int> vec( 12UL );
   CustomMatrix<int,unaligned,unpadded> A( &vec[0], 3UL, 4UL );

   // Definition of a custom 8x12 matrix for an aligned and padded integer array of
   // capacity 128 (including 8 padding elements per row). Note that the std::unique_ptr
   // must be guaranteed to outlive the custom matrix!
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 128UL ) );
   CustomMatrix<int,aligned,padded> B( memory.get(), 8UL, 12UL, 16UL );
   \endcode

// \n \subsection custommatrix_copy_operations Copy Operations
//
// As with all dense matrices it is possible to copy construct a custom matrix:

   \code
   using blaze::CustomMatrix;
   using blaze::unaligned;
   using blaze::unpadded;

   using CustomType = CustomMatrix<int,unaligned,unpadded>;

   std::vector<int> vec( 6UL, 10 );    // Vector of 6 integers of the value 10
   CustomType A( &vec[0], 2UL, 3UL );  // Represent the std::vector as Blaze dense matrix
   a[1] = 20;                          // Also modifies the std::vector

   CustomType B( a );  // Creating a copy of vector a
   b[2] = 20;          // Also affects matrix A and the std::vector
   \endcode

// It is important to note that a custom matrix acts as a reference to the specified array. Thus
// the result of the copy constructor is a new custom matrix that is referencing and representing
// the same array as the original custom matrix.
//
// In contrast to copy construction, just as with references, copy assignment does not change
// which array is referenced by the custom matrices, but modifies the values of the array:

   \code
   std::vector<int> vec2( 6UL, 4 );     // Vector of 6 integers of the value 4
   CustomType C( &vec2[0], 2UL, 3UL );  // Represent the std::vector as Blaze dense matrix

   A = C;  // Copy assignment: Set all values of matrix A and B to 4.
   \endcode

// \n \subsection custommatrix_alignment Alignment
//
// In case the custom matrix is specified as \c aligned the passed array must adhere to some
// alignment restrictions based on the alignment requirements of the used data type and the
// used instruction set (SSE, AVX, ...). The restriction applies to the first element of each
// row/column: In case of a row-major matrix the first element of each row must be properly
// aligned, in case of a column-major matrix the first element of each column must be properly
// aligned. For instance, if a row-major matrix is used and AVX is active the first element of
// each row must be 32-bit aligned:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;
   using blaze::rowMajor;

   // Allocation of 32-bit aligned memory
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 40UL ) );

   CustomMatrix<int,aligned,padded,rowMajor> A( memory.get(), 5UL, 6UL, 8UL );
   \endcode

// In the example, the row-major matrix has six columns. However, since with AVX eight integer
// values are loaded together the matrix is padded with two additional elements. This guarantees
// that the first element of each row is 32-bit aligned. In case the alignment requirements are
// violated, a \c std::invalid_argument exception is thrown.
//
// \n \subsection custommatrix_padding Padding
//
// Adding padding elements to the end of an array can have a significant impact on performance.
// For instance, assuming that AVX is available, then two aligned, padded, 3x3 double precision
// matrices can be added via three SIMD addition operations:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;

   using CustomType = CustomMatrix<double,aligned,padded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 12UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 12UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 12UL ) );

   // Creating padded custom 3x3 matrix with an additional padding element in each row
   CustomType A( memory1.get(), 3UL, 3UL, 4UL );
   CustomType B( memory2.get(), 3UL, 3UL, 4UL );
   CustomType C( memory3.get(), 3UL, 3UL, 4UL );

   // ... Initialization

   C = A + B;  // AVX-based matrix addition
   \endcode

// In this example, maximum performance is possible. However, in case no padding elements are
// inserted a scalar addition has to be used:

   \code
   using blaze::CustomMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   using CustomType = CustomMatrix<double,aligned,unpadded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 9UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 9UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 9UL ) );

   // Creating unpadded custom 3x3 matrix
   CustomType A( memory1.get(), 3UL, 3UL );
   CustomType B( memory2.get(), 3UL, 3UL );
   CustomType C( memory3.get(), 3UL, 3UL );

   // ... Initialization

   C = A + B;  // Scalar matrix addition
   \endcode

// Note that the construction of padded and unpadded aligned matrices looks identical. However,
// in case of padded matrices, \b Blaze will zero initialize the padding element and use them
// in all computations in order to achieve maximum performance. In case of an unpadded matrix
// \b Blaze will ignore the elements with the downside that it is not possible to load a complete
// row to an AVX register, which makes it necessary to fall back to a scalar addition.
//
// The number of padding elements is required to be sufficient with respect to the available
// instruction set: In case of an aligned padded custom matrix the added padding elements must
// guarantee that the total number of elements in each row/column is a multiple of the SIMD
// vector width. In case of an unaligned padded matrix the number of padding elements can be
// greater or equal the number of padding elements of an aligned padded custom matrix. In case
// the padding is insufficient with respect to the available instruction set, a
// \c std::invalid_argument exception is thrown.
//
//
// \n \section custommatrix_arithmetic_operations Arithmetic Operations
//
// The use of custom matrices in arithmetic operations is designed to be as natural and intuitive
// as possible. All operations (addition, subtraction, multiplication, scaling, ...) can be
// expressed similar to a text book representation. Also, custom matrices can be combined with all
// other dense and sparse vectors and matrices. The following example gives an impression of the
// use of CustomMatrix:

   \code
   using blaze::CustomMatrix;
   using blaze::CompressedMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Non-initialized custom 2x3 row-major matrix. All given arrays are considered to be
   // unaligned and unpadded. The memory is managed via a 'std::vector'.
   std::vector<double> memory1( 6UL );
   CustomMatrix<double,unaligned,unpadded> A( memory1.data(), 2UL, 3UL );

   A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;  // Initialization of the first row
   A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;  // Initialization of the second row

   // Non-initialized custom 2x3 row-major matrix with padding elements. All given arrays are
   // required to be properly aligned and padded. The memory is managed via a 'std::unique_ptr'.
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 16UL ) );
   CustomMatrix<double,aligned,padded> B( memory2.get(), 2UL, 3UL, 8UL );

   B(0,0) = 1.0; B(0,1) = 3.0; B(0,2) = 5.0;    // Initialization of the first row
   B(1,0) = 2.0; B(1,1) = 4.0; B(1,2) = 6.0;    // Initialization of the second row

   CompressedMatrix<float> C( 2, 3 );        // Empty row-major sparse single precision matrix
   DynamicMatrix<float>    D( 3, 2, 4.0F );  // Directly, homogeneously initialized single precision 3x2 matrix

   DynamicMatrix<double,rowMajor>    E( A );  // Creation of a new row-major matrix as a copy of A
   DynamicMatrix<double,columnMajor> F;       // Creation of a default column-major matrix

   E = A + B;     // Matrix addition and assignment to a row-major matrix
   F = A - C;     // Matrix subtraction and assignment to a column-major matrix
   F = A * D;     // Matrix multiplication between two matrices of different element types

   A *= 2.0;      // In-place scaling of matrix A
   E  = 2.0 * B;  // Scaling of matrix B
   F  = D * 2.0;  // Scaling of matrix D

   E += A - B;    // Addition assignment
   E -= A + C;    // Subtraction assignment
   F *= A * D;    // Multiplication assignment
   \endcode
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
class CustomMatrix
   : public DenseMatrix< CustomMatrix<Type,AF,PF,SO,Tag,RT>, SO >
{
 public:
   //**Type definitions****************************************************************************
   using This     = CustomMatrix<Type,AF,PF,SO,Tag,RT>;  //!< Type of this CustomMatrix instance.
   using BaseType = DenseMatrix<This,SO>;                //!< Base type of this CustomMatrix instance.

   //! Result type for expression template evaluations.
   using ResultType = RT;

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = CustomOppositeType_t<RT>;

   //! Transpose type for expression template evaluations.
   using TransposeType = CustomTransposeType_t<RT>;

   using ElementType   = Type;                      //!< Type of the matrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using TagType       = Tag;                       //!< Tag type of this CustomMatrix instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = Type*;        //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = DenseIterator<Type,AF>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,AF>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CustomMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using RRT   = Rebind_t< RT, RemoveConst_t<NewType> >;  //!< The rebound result type.
      using Other = CustomMatrix<NewType,AF,PF,SO,Tag,RRT>;  //!< The type of the other CustomMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CustomMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using RRT   = Resize_t<RT,NewM,NewN>;               //!< The resized result type.
      using Other = CustomMatrix<Type,AF,PF,SO,Tag,RRT>;  //!< The type of the other CustomMatrix.
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
   inline CustomMatrix();
   inline CustomMatrix( Type* ptr, size_t m, size_t n );
   inline CustomMatrix( Type* ptr, size_t m, size_t n, size_t nn );
   inline CustomMatrix( const CustomMatrix& m );
   inline CustomMatrix( CustomMatrix&& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CustomMatrix();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j ) noexcept;
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
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
   inline CustomMatrix& operator=( const Type& set );
   inline CustomMatrix& operator=( initializer_list< initializer_list<Type> > list );

   template< typename Other, size_t Rows, size_t Cols >
   inline CustomMatrix& operator=( const Other (&array)[Rows][Cols] );

   template< typename Other, size_t Rows, size_t Cols >
   inline CustomMatrix& operator=( const std::array<std::array<Other,Cols>,Rows>& array );

   inline CustomMatrix& operator=( const CustomMatrix& rhs );
   inline CustomMatrix& operator=( CustomMatrix&& rhs ) noexcept;

   template< typename MT, bool SO2 > inline CustomMatrix& operator= ( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline CustomMatrix& operator+=( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline CustomMatrix& operator-=( const Matrix<MT,SO2>& rhs );
   template< typename MT, bool SO2 > inline CustomMatrix& operator%=( const Matrix<MT,SO2>& rhs );
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
   inline void   reset();
   inline void   reset( size_t i );
   inline void   clear();
   inline void   swap( CustomMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline CustomMatrix& transpose();
   inline CustomMatrix& ctranspose();

   template< typename Other > inline CustomMatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Resource management functions***************************************************************
   /*!\name Resource management functions */
   //@{
   inline void reset( Type* ptr, size_t m, size_t n );
   inline void reset( Type* ptr, size_t m, size_t n, size_t nn );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT::simdEnabled &&
        IsSIMDCombinable_v< Type, ElementType_t<MT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDAdd_v< Type, ElementType_t<MT> > &&
        !IsDiagonal_v<MT> );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDSub_v< Type, ElementType_t<MT> > &&
        !IsDiagonal_v<MT> );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDMult_v< Type, ElementType_t<MT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
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

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT >
   inline auto assign( const DenseMatrix<MT,SO>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT> >;

   template< typename MT >
   inline auto assign( const DenseMatrix<MT,SO>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT> >;

   template< typename MT > inline void assign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT >
   inline auto addAssign( const DenseMatrix<MT,SO>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT> >;

   template< typename MT >
   inline auto addAssign( const DenseMatrix<MT,SO>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT> >;

   template< typename MT > inline void addAssign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT >
   inline auto subAssign( const DenseMatrix<MT,SO>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT> >;

   template< typename MT >
   inline auto subAssign( const DenseMatrix<MT,SO>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT> >;

   template< typename MT > inline void subAssign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,!SO>& rhs );

   template< typename MT >
   inline auto schurAssign( const DenseMatrix<MT,SO>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT> >;

   template< typename MT >
   inline auto schurAssign( const DenseMatrix<MT,SO>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT> >;

   template< typename MT > inline void schurAssign( const DenseMatrix<MT,!SO>&  rhs );
   template< typename MT > inline void schurAssign( const SparseMatrix<MT,SO>&  rhs );
   template< typename MT > inline void schurAssign( const SparseMatrix<MT,!SO>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;   //!< The current number of rows of the matrix.
   size_t n_;   //!< The current number of columns of the matrix.
   size_t nn_;  //!< The number of elements between two rows.
   Type* v_;    //!< The custom array of elements.
                /*!< Access to the matrix elements is gained via the function call
                     operator. In case of row-major order the memory layout of the
                     elements is
                     \f[\left(\begin{array}{*{5}{c}}
                     0            & 1             & 2             & \cdots & N-1         \\
                     N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                     \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                     M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                     \end{array}\right)\f]. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
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
/*!\brief The default constructor for CustomMatrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::CustomMatrix()
   : m_ ( 0UL )      // The current number of rows of the matrix
   , n_ ( 0UL )      // The current number of columns of the matrix
   , nn_( 0UL )      // The number of elements between two rows
   , v_ ( nullptr )  // The custom array of elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This constructor creates an unpadded custom matrix of size \f$ m \times n \f$. The construction
// fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This constructor is \b NOT available for padded custom matrices!
// \note The custom matrix does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::CustomMatrix( Type* ptr, size_t m, size_t n )
   : m_ ( m )    // The current number of rows of the matrix
   , n_ ( n )    // The current number of columns of the matrix
   , nn_( n )    // The number of elements between two rows
   , v_ ( ptr )  // The custom array of elements
{
   BLAZE_STATIC_ASSERT( PF == unpadded );

   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && ( !checkAlignment( ptr ) || nn_ % SIMDSIZE != 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \param nn The total number of elements between two rows/columns.
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This constructor creates a custom matrix of size \f$ m \times n \f$. The construction fails
// if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...);
//  - ... the specified spacing \a nn is insufficient for the given data type \a Type and the
//    available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note The custom matrix does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::CustomMatrix( Type* ptr, size_t m, size_t n, size_t nn )
   : m_ ( m )    // The current number of rows of the matrix
   , n_ ( n )    // The current number of columns of the matrix
   , nn_( nn )   // The number of elements between two rows
   , v_ ( ptr )  // The custom array of elements
{
   using blaze::clear;

   using ClearFunctor = If_t< PF || !IsConst_v<Type>, Clear, Noop >;

   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && ( !checkAlignment( ptr ) || nn_ % SIMDSIZE != 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }

   if( PF && IsVectorizable_v<Type> && ( nn_ < nextMultiple<size_t>( n_, SIMDSIZE ) ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Insufficient capacity for padded matrix" );
   }

   if( PF && IsVectorizable_v<Type> ) {
      ClearFunctor clear;
      for( size_t i=0UL; i<m_; ++i ) {
         for( size_t j=n_; j<nn_; ++j )
            clear( v_[i*nn_+j] );
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for CustomMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor initializes the custom matrix as an exact copy of the given custom matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::CustomMatrix( const CustomMatrix& m )
   : m_ ( m.m_ )   // The current number of rows of the matrix
   , n_ ( m.n_ )   // The current number of columns of the matrix
   , nn_( m.nn_ )  // The number of elements between two rows
   , v_ ( m.v_ )   // The custom array of elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for CustomMatrix.
//
// \param m The matrix to be moved into this instance.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::CustomMatrix( CustomMatrix&& m ) noexcept
   : m_ ( m.m_ )   // The current number of rows of the matrix
   , n_ ( m.n_ )   // The current number of columns of the matrix
   , nn_( m.nn_ )  // The number of elements between two rows
   , v_ ( m.v_ )   // The custom array of elements
{
   m.m_  = 0UL;
   m.n_  = 0UL;
   m.nn_ = 0UL;
   m.v_  = nullptr;

   BLAZE_INTERNAL_ASSERT( m.data() == nullptr, "Invalid data reference detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for CustomMatrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>::~CustomMatrix()
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST( RT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( RT );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Reference
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator()( size_t i, size_t j ) noexcept
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i*nn_+j];
}
//*************************************************************************************************


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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstReference
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i*nn_+j];
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Reference
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::at( size_t i, size_t j )
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
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstReference
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::at( size_t i, size_t j ) const
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Pointer
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::data() noexcept
{
   return v_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstPointer
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::data() const noexcept
{
   return v_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Pointer
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::data( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_+i*nn_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstPointer
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::data( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_+i*nn_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Iterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::begin( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_+i*nn_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_+i*nn_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_+i*nn_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::Iterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::end( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_+i*nn_+n_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_+i*nn_+n_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_+i*nn_+n_ );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( const Type& rhs )
{
   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List assignment to all matrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to static matrix.
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   blaze::CustomMatrix<int,unaligned,unpadded,rowMajor> A( array, 3UL, 3UL );
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix elements are assigned the values from the given initializer list. Missing values
// are initialized as default (as e.g. the value 6 in the example). Note that in case the size
// of the top-level initializer list exceeds the number of rows or the size of any nested list
// exceeds the number of columns, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( initializer_list< initializer_list<Type> > list )
{
   if( list.size() != m_ || determineColumns( list ) > n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to custom matrix" );
   }

   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), v_+i*nn_ ),
                 v_+i*nn_+( PF ? nn_ : n_ ), Type() );
      ++i;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all matrix elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::CustomMatrix<int,unaligned,unpadded,rowMajor> A( array, 3UL, 3UL );
   A = init;
   \endcode

// The matrix is assigned the values from the given static array. Missing values are initialized
// with default values (as e.g. the value 6 in the example). Note that the size of the static
// array must match the size of the custom matrix. Otherwise a \a std::invalid_argument exception
// is thrown. Also note that after the assignment \a array will have the same entries as \a init.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Rows       // Number of rows of the static array
        , size_t Cols >     // Number of columns of the static array
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( const Other (&array)[Rows][Cols] )
{
   if( m_ != Rows || n_ != Cols ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Rows; ++i )
      for( size_t j=0UL; j<Cols; ++j )
         v_[i*nn_+j] = array[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all matrix elements.
//
// \param array The given std::array for the assignment.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } } };
   blaze::CustomMatrix<int,unaligned,unpadded,rowMajor> A( array, 3UL, 3UL );
   A = init;
   \endcode

// The matrix is assigned the values from the given std::array. Missing values are initialized
// with default values (as e.g. the value 6 in the example). Note that the size of the std::array
// must match the size of the custom matrix. Otherwise a \a std::invalid_argument exception is
// thrown. Also note that after the assignment \a array will have the same entries as \a init.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Rows       // Number of rows of the static array
        , size_t Cols >     // Number of columns of the static array
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( const std::array<std::array<Other,Cols>,Rows>& array )
{
   if( m_ != Rows || n_ != Cols ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Rows; ++i )
      for( size_t j=0UL; j<Cols; ++j )
         v_[i*nn_+j] = array[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for CustomMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The matrix is initialized as a copy of the given matrix. In case the current sizes of the two
// matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( const CustomMatrix& rhs )
{
   if( rhs.rows() != m_ || rhs.columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   smpAssign( *this, *rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for CustomMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( CustomMatrix&& rhs ) noexcept
{
   m_  = rhs.m_;
   n_  = rhs.n_;
   nn_ = rhs.nn_;
   v_  = rhs.v_;

   rhs.m_  = 0UL;
   rhs.n_  = 0UL;
   rhs.nn_ = 0UL;
   rhs.v_  = nullptr;

   BLAZE_INTERNAL_ASSERT( rhs.data() == nullptr, "Invalid data reference detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The matrix is initialized as a copy of the given matrix. In case the current sizes of the two
// matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator=( const Matrix<MT,SO2>& rhs )
{
   using TT = decltype( trans( *this ) );
   using CT = decltype( ctrans( *this ) );
   using IT = decltype( inv( *this ) );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsSame_v<MT,TT> && (*rhs).isAliased( this ) ) {
      transpose();
   }
   else if( IsSame_v<MT,CT> && (*rhs).isAliased( this ) ) {
      ctranspose();
   }
   else if( !IsSame_v<MT,IT> && (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseMatrix_v<MT> )
         reset();
      smpAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator+=( const Matrix<MT,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator-=( const Matrix<MT,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO2 >     // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::operator%=( const Matrix<MT,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpSchurAssign( *this, tmp );
   }
   else {
      smpSchurAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::rows() const noexcept
{
   return m_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::columns() const noexcept
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::spacing() const noexcept
{
   return nn_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::capacity() const noexcept
{
   return m_ * nn_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return nn_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         if( !isDefault<strict>( v_[i*nn_+j] ) )
            ++nonzeros;

   return nonzeros;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,SO,Tag,RT>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( i*nn_ + n_ );
   size_t nonzeros( 0UL );

   for( size_t j=i*nn_; j<jend; ++j )
      if( !isDefault<strict>( v_[j] ) )
         ++nonzeros;

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::reset()
{
   using blaze::clear;

   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         clear( v_[i*nn_+j] );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::reset( size_t i )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   for( size_t j=0UL; j<n_; ++j )
      clear( v_[i*nn_+j] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the \f$ M \times N \f$ matrix.
//
// \return void
//
// After the clear() function, the size of the matrix is 0.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::clear()
{
   m_  = 0UL;
   n_  = 0UL;
   nn_ = 0UL;
   v_  = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::swap( CustomMatrix& m ) noexcept
{
   using std::swap;

   swap( m_ , m.m_  );
   swap( n_ , m.n_  );
   swap( nn_, m.nn_ );
   swap( v_ , m.v_  );
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
// \exception std::logic_error Impossible transpose operation.
//
// In case the matrix is not a square matrix, a \a std::logic_error exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>& CustomMatrix<Type,AF,PF,SO,Tag,RT>::transpose()
{
   using std::swap;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Impossible transpose operation" );
   }

   for( size_t i=1UL; i<m_; ++i )
      for( size_t j=0UL; j<i; ++j )
         swap( v_[i*nn_+j], v_[j*nn_+i] );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the matrix.
//
// \return Reference to the transposed matrix.
// \exception std::logic_error Impossible transpose operation.
//
// In case the matrix is not a square matrix, a \a std::logic_error exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>& CustomMatrix<Type,AF,PF,SO,Tag,RT>::ctranspose()
{
   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Impossible transpose operation" );
   }

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<i; ++j ) {
         cswap( v_[i*nn_+j], v_[j*nn_+i] );
      }
      conjugate( v_[i*nn_+i] );
   }

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
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   CustomMatrix<int,unaligned,unpadded> A( ... );

   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the scalar value
inline CustomMatrix<Type,AF,PF,SO,Tag,RT>&
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] *= scalar;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  RESOURCE MANAGEMENT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Resets the custom matrix and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \return void
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This function resets the custom matrix to the given array of elements of size \f$ m \times n \f$.
// The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function is \b NOT available for padded custom matrices!
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom matrix referencing the array goes out of scope.
// \note The custom matrix does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::reset( Type* ptr, size_t m, size_t n )
{
   BLAZE_STATIC_ASSERT( PF == unpadded );

   CustomMatrix tmp( ptr, m, n );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resets the custom matrix and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \param nn The total number of elements between two rows/columns.
// \return void
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This function resets the custom matrix to the given array of elements of size \f$ m \times n \f$.
// The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom matrix referencing the array goes out of scope.
// \note The custom matrix does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::reset( Type* ptr, size_t m, size_t n, size_t nn )
{
   CustomMatrix tmp( ptr, m, n, nn );
   swap( tmp );
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
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomMatrix<Type,AF,PF,SO,Tag,RT>::canAlias( const Other* alias ) const noexcept
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
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomMatrix<Type,AF,PF,SO,Tag,RT>::isAliased( const Other* alias ) const noexcept
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomMatrix<Type,AF,PF,SO,Tag,RT>::isAligned() const noexcept
{
   return ( AF || ( checkAlignment( v_ ) && columns() % SIMDSIZE == 0UL ) );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomMatrix<Type,AF,PF,SO,Tag,RT>::canSMPAssign() const noexcept
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::load( size_t i, size_t j ) const noexcept
{
   if( AF && PF )
      return loada( i, j );
   else
      return loadu( i, j );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::loada( size_t i, size_t j ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= ( PF ? nn_ : n_ ), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !PF || j % SIMDSIZE == 0UL, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i*nn_+j ), "Invalid alignment detected" );

   return loada( v_+i*nn_+j );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,SO,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::loadu( size_t i, size_t j ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= ( PF ? nn_ : n_ ), "Invalid column access index" );

   return loadu( v_+i*nn_+j );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense matrix. The row index
// must be smaller than the number of rows and the column index must be smaller than the number
// of columns. Additionally, the column index (in case of a row-major matrix) or the row index
// (in case of a column-major matrix) must be a multiple of the number of values inside the
// SIMD element. This function must \b NOT be called explicitly! It is used internally for the
// performance optimized evaluation of expression templates. Calling this function explicitly
// might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( AF && PF )
      storea( i, j, value );
   else
      storeu( i, j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= ( PF ? nn_ : n_ ), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !PF || j % SIMDSIZE == 0UL, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i*nn_+j ), "Invalid alignment detected" );

   storea( v_+i*nn_+j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= ( PF ? nn_ : n_ ), "Invalid column access index" );

   storeu( v_+i*nn_+j, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense matrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of a
// row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,SO,Tag,RT>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= ( PF ? nn_ : n_ ), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !PF || j % SIMDSIZE == 0UL, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i*nn_+j ), "Invalid alignment detected" );

   stream( v_+i*nn_+j, value );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::assign( const DenseMatrix<MT,SO>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( n_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= n_, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         v_[i*nn_+j    ] = (*rhs)(i,j    );
         v_[i*nn_+j+1UL] = (*rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         v_[i*nn_+jpos] = (*rhs)(i,jpos);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::assign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   const size_t jpos( remainder ? prevMultiple( n_, SIMDSIZE ): n_ );
   BLAZE_INTERNAL_ASSERT( jpos <= n_, "Invalid end calculation" );

   if( AF && PF && useStreaming &&
       ( m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) ) && !(*rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<m_; ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT> right( (*rhs).begin(i) );

         for( ; j<jpos; j+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE, right += SIMDSIZE;
         }
         for( ; remainder && j<n_; ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<m_; ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT> right( (*rhs).begin(i) );

         for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<jpos; j+=SIMDSIZE ) {
            left.store( right.load() ); left+=SIMDSIZE, right+=SIMDSIZE;
         }
         for( ; remainder && j<n_; ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::assign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( min( m_, ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( min( n_, jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               v_[i*nn_+j] = (*rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::assign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i*nn_+element->index()] = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::assign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()*nn_+j] = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::addAssign( const DenseMatrix<MT,SO>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
   {
      if( IsDiagonal_v<MT> )
      {
         v_[i*nn_+i] += (*rhs)(i,i);
      }
      else
      {
         const size_t jbegin( ( IsUpper_v<MT> )
                              ?( IsStrictlyUpper_v<MT> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend  ( ( IsLower_v<MT> )
                              ?( IsStrictlyLower_v<MT> ? i : i+1UL )
                              :( n_ ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         size_t j( jbegin );

         for( ; (j+2UL) <= jend; j+=2UL ) {
            v_[i*nn_+j    ] += (*rhs)(i,j    );
            v_[i*nn_+j+1UL] += (*rhs)(i,j+1UL);
         }
         if( j < jend ) {
            v_[i*nn_+j] += (*rhs)(i,j);
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::addAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t i=0UL; i<m_; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( IsStrictlyLower_v<MT> ? i : i+1UL )
                           :( n_ ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( remainder ? prevMultiple( jend, SIMDSIZE ) : jend );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && j<jend; ++j ) {
         *left += *right; ++left; ++right;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::addAssign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( min( m_, ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block )
      {
         if( IsLower_v<MT> && ii < jj ) break;
         if( IsUpper_v<MT> && ii > jj ) continue;

         for( size_t i=ii; i<iend; ++i )
         {
            const size_t jbegin( ( IsUpper_v<MT> )
                                 ?( max( ( IsStrictlyUpper_v<MT> ? i+1UL : i ), jj ) )
                                 :( jj ) );
            const size_t jend  ( ( IsLower_v<MT> )
                                 ?( min( ( IsStrictlyLower_v<MT> ? i : i+1UL ), n_, jj+block ) )
                                 :( min( n_, jj+block ) ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; ++j ) {
               v_[i*nn_+j] += (*rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::addAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i*nn_+element->index()] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::addAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()*nn_+j] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::subAssign( const DenseMatrix<MT,SO>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
   {
      if( IsDiagonal_v<MT> )
      {
         v_[i*nn_+i] -= (*rhs)(i,i);
      }
      else
      {
         const size_t jbegin( ( IsUpper_v<MT> )
                              ?( IsStrictlyUpper_v<MT> ? i+1UL : i )
                              :( 0UL ) );
         const size_t jend  ( ( IsLower_v<MT> )
                              ?( IsStrictlyLower_v<MT> ? i : i+1UL )
                              :( n_ ) );
         BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

         size_t j( jbegin );

         for( ; (j+2UL) <= jend; j+=2UL ) {
            v_[i*nn_+j    ] -= (*rhs)(i,j    );
            v_[i*nn_+j+1UL] -= (*rhs)(i,j+1UL);
         }
         if( j < jend ) {
            v_[i*nn_+j] -= (*rhs)(i,j);
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::subAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t i=0UL; i<m_; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( IsStrictlyLower_v<MT> ? i : i+1UL )
                           :( n_ ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( remainder ? prevMultiple( jend, SIMDSIZE ) : jend );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && j<jend; ++j ) {
         *left -= *right; ++left; ++right;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::subAssign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( min( m_, ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block )
      {
         if( IsLower_v<MT> && ii < jj ) break;
         if( IsUpper_v<MT> && ii > jj ) continue;

         for( size_t i=ii; i<iend; ++i )
         {
            const size_t jbegin( ( IsUpper_v<MT> )
                                 ?( max( ( IsStrictlyUpper_v<MT> ? i+1UL : i ), jj ) )
                                 :( jj ) );
            const size_t jend  ( ( IsLower_v<MT> )
                                 ?( min( ( IsStrictlyLower_v<MT> ? i : i+1UL ), n_, jj+block ) )
                                 :( min( n_, jj+block ) ) );
            BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

            for( size_t j=jbegin; j<jend; ++j ) {
               v_[i*nn_+j] -= (*rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::subAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i*nn_+element->index()] -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::subAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()*nn_+j] -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::schurAssign( const DenseMatrix<MT,SO>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( n_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= n_, "Invalid end calculation" );

   for( size_t i=0UL; i<m_; ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         v_[i*nn_+j    ] *= (*rhs)(i,j    );
         v_[i*nn_+j+1UL] *= (*rhs)(i,j+1UL);
      }
      if( jpos < n_ ) {
         v_[i*nn_+jpos] *= (*rhs)(i,jpos);
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,SO,Tag,RT>::schurAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t i=0UL; i<m_; ++i )
   {
      const size_t jpos( remainder ? prevMultiple( n_, SIMDSIZE ): n_ );
      BLAZE_INTERNAL_ASSERT( jpos <= n_, "Invalid end calculation" );

      size_t j( 0UL );
      Iterator left( begin(i) );
      ConstIterator_t<MT> right( (*rhs).begin(i) );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && j<n_; ++j ) {
         *left *= *right; ++left; ++right;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::schurAssign( const DenseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<m_; ii+=block ) {
      const size_t iend( min( m_, ii+block ) );
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( min( n_, jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               v_[i*nn_+j] *= (*rhs)(i,j);
            }
         }
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::schurAssign( const SparseMatrix<MT,SO>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( size_t i=0UL; i<m_; ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i*nn_+element->index()] = tmp(i,element->index()) * element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,SO,Tag,RT>::schurAssign( const SparseMatrix<MT,!SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( size_t j=0UL; j<n_; ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()*nn_+j] = tmp(element->index(),j) * element->value();
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR COLUMN-MAJOR MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of CustomMatrix for column-major matrices.
// \ingroup custom_matrix
//
// This specialization of CustomMatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
class CustomMatrix<Type,AF,PF,true,Tag,RT>
   : public DenseMatrix< CustomMatrix<Type,AF,PF,true,Tag,RT>, true >
{
 public:
   //**Type definitions****************************************************************************
   using This     = CustomMatrix<Type,AF,PF,true,Tag,RT>;  //!< Type of this CustomMatrix instance.
   using BaseType = DenseMatrix<This,true>;                //!< Base type of this CustomMatrix instance.

   //! Result type for expression template evaluations.
   using ResultType = RT;

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = CustomOppositeType_t<RT>;

   //! Transpose type for expression template evaluations.
   using TransposeType = CustomTransposeType_t<RT>;

   using ElementType   = Type;                      //!< Type of the matrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using TagType       = Tag;                       //!< Tag type of this CustomMatrix instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = Type*;        //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = DenseIterator<Type,AF>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,AF>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CustomMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind {
      using RRT   = Rebind_t< RT, RemoveConst_t<NewType> >;    //!< The rebound result type.
      using Other = CustomMatrix<NewType,AF,PF,true,Tag,RRT>;  //!< The type of the other CustomMatrix.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CustomMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize {
      using RRT   = Resize_t<RT,NewM,NewN>;                 //!< The resized result type.
      using Other = CustomMatrix<Type,AF,PF,true,Tag,RRT>;  //!< The type of the other CustomMatrix.
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
   inline CustomMatrix();
   inline CustomMatrix( Type* ptr, size_t m, size_t n );
   inline CustomMatrix( Type* ptr, size_t m, size_t n, size_t mm );
   inline CustomMatrix( const CustomMatrix& m );
   inline CustomMatrix( CustomMatrix&& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CustomMatrix();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j ) noexcept;
   inline ConstReference operator()( size_t i, size_t j ) const noexcept;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
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
   inline CustomMatrix& operator=( const Type& set );
   inline CustomMatrix& operator=( initializer_list< initializer_list<Type> > list );

   template< typename Other, size_t Rows, size_t Cols >
   inline CustomMatrix& operator=( const Other (&array)[Rows][Cols] );

   template< typename Other, size_t Rows, size_t Cols >
   inline CustomMatrix& operator=( const std::array<std::array<Other,Cols>,Rows>& array );

   inline CustomMatrix& operator=( const CustomMatrix& rhs );
   inline CustomMatrix& operator=( CustomMatrix&& rhs ) noexcept;

   template< typename MT, bool SO > inline CustomMatrix& operator= ( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline CustomMatrix& operator+=( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline CustomMatrix& operator-=( const Matrix<MT,SO>& rhs );
   template< typename MT, bool SO > inline CustomMatrix& operator%=( const Matrix<MT,SO>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t rows() const noexcept;
   inline size_t columns() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t j ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t j ) const;
   inline void   reset();
   inline void   reset( size_t j );
   inline void   clear();
   inline void   swap( CustomMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline CustomMatrix& transpose();
   inline CustomMatrix& ctranspose();

   template< typename Other > inline CustomMatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Resource management functions***************************************************************
   /*!\name Resource management functions */
   //@{
   inline void reset( Type* ptr, size_t m, size_t n );
   inline void reset( Type* ptr, size_t m, size_t n, size_t mm );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT::simdEnabled &&
        IsSIMDCombinable_v< Type, ElementType_t<MT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDAdd_v< Type, ElementType_t<MT> > &&
        !IsDiagonal_v<MT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDSub_v< Type, ElementType_t<MT> > &&
        !IsDiagonal_v<MT> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT> &&
        HasSIMDMult_v< Type, ElementType_t<MT> > );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
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

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT >
   inline auto assign( const DenseMatrix<MT,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT> >;

   template< typename MT >
   inline auto assign( const DenseMatrix<MT,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT> >;

   template< typename MT > inline void assign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void assign( const SparseMatrix<MT,false>& rhs );

   template< typename MT >
   inline auto addAssign( const DenseMatrix<MT,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT> >;

   template< typename MT >
   inline auto addAssign( const DenseMatrix<MT,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT> >;

   template< typename MT > inline void addAssign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void addAssign( const SparseMatrix<MT,false>& rhs );

   template< typename MT >
   inline auto subAssign ( const DenseMatrix<MT,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT> >;

   template< typename MT >
   inline auto subAssign ( const DenseMatrix<MT,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT> >;

   template< typename MT > inline void subAssign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void subAssign( const SparseMatrix<MT,false>& rhs );

   template< typename MT >
   inline auto schurAssign ( const DenseMatrix<MT,true>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT> >;

   template< typename MT >
   inline auto schurAssign ( const DenseMatrix<MT,true>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT> >;

   template< typename MT > inline void schurAssign( const DenseMatrix<MT,false>&  rhs );
   template< typename MT > inline void schurAssign( const SparseMatrix<MT,true>&  rhs );
   template< typename MT > inline void schurAssign( const SparseMatrix<MT,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;   //!< The current number of rows of the matrix.
   size_t mm_;  //!< The number of elements between two columns.
   size_t n_;   //!< The current number of columns of the matrix.
   Type* v_;    //!< The custom array of elements.
                /*!< Access to the matrix elements is gained via the function
                     call operator. */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
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
/*!\brief The default constructor for CustomMatrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::CustomMatrix()
   : m_ ( 0UL )      // The current number of rows of the matrix
   , mm_( 0UL )      // The number of elements between two columns
   , n_ ( 0UL )      // The current number of columns of the matrix
   , v_ ( nullptr )  // The custom array of elements
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ m \times n \f$.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This constructor creates an unpadded custom matrix of size \f$ m \times n \f$. The construction
// fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This constructor is \b NOT available for padded custom matrices!
// \note The custom matrix does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::CustomMatrix( Type* ptr, size_t m, size_t n )
   : m_ ( m )    // The current number of rows of the matrix
   , mm_( m )    // The number of elements between two columns
   , n_ ( n )    // The current number of columns of the matrix
   , v_ ( ptr )  // The custom array of elements
{
   BLAZE_STATIC_ASSERT( PF == unpadded );

   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && ( !checkAlignment( ptr ) || mm_ % SIMDSIZE != 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ m \times n \f$.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param mm The total number of elements between two columns.
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This constructor creates a custom matrix of size \f$ m \times n \f$. The construction fails
// if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...);
//  - ... the specified spacing \a mm is insufficient for the given data type \a Type and the
//    available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note The custom matrix does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::CustomMatrix( Type* ptr, size_t m, size_t n, size_t mm )
   : m_ ( m )    // The current number of rows of the matrix
   , mm_( mm )   // The number of elements between two columns
   , n_ ( n )    // The current number of columns of the matrix
   , v_ ( ptr )  // The custom array of elements
{
   using blaze::clear;

   using ClearFunctor = If_t< PF || !IsConst_v<Type>, Clear, Noop >;

   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && ( !checkAlignment( ptr ) || mm_ % SIMDSIZE != 0UL ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }

   if( PF && IsVectorizable_v<Type> && ( mm_ < nextMultiple<size_t>( m_, SIMDSIZE ) ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Insufficient capacity for padded matrix" );
   }

   if( PF && IsVectorizable_v<Type> ) {
      ClearFunctor clear;
      for( size_t j=0UL; j<n_; ++j )
         for( size_t i=m_; i<mm_; ++i ) {
            clear( v_[i+j*mm_] );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for CustomMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::CustomMatrix( const CustomMatrix& m )
   : m_ ( m.m_ )   // The current number of rows of the matrix
   , mm_( m.mm_ )  // The number of elements between two columns
   , n_ ( m.n_ )   // The current number of columns of the matrix
   , v_ ( m.v_ )   // The custom array of elements
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The move constructor for CustomMatrix.
//
// \param m The matrix to be moved into this instance.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::CustomMatrix( CustomMatrix&& m ) noexcept
   : m_ ( m.m_ )   // The current number of rows of the matrix
   , mm_( m.mm_ )  // The number of elements between two columns
   , n_ ( m.n_ )   // The current number of columns of the matrix
   , v_ ( m.v_ )   // The custom array of elements
{
   m.m_  = 0UL;
   m.mm_ = 0UL;
   m.n_  = 0UL;
   m.v_  = nullptr;

   BLAZE_INTERNAL_ASSERT( m.data() == nullptr, "Invalid data reference detected" );
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
/*!\brief The destructor for CustomMatrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>::~CustomMatrix()
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST( RT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( RT );
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
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Reference
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator()( size_t i, size_t j ) noexcept
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i+j*mm_];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the matrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstReference
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator()( size_t i, size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( i<m_, "Invalid row access index"    );
   BLAZE_USER_ASSERT( j<n_, "Invalid column access index" );
   return v_[i+j*mm_];
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Reference
   CustomMatrix<Type,AF,PF,true,Tag,RT>::at( size_t i, size_t j )
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstReference
   CustomMatrix<Type,AF,PF,true,Tag,RT>::at( size_t i, size_t j ) const
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
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data. Whereas the number of
// elements within a column are given by the \c columns() member functions, the total number
// of elements including padding is given by the \c spacing() member function.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Pointer
   CustomMatrix<Type,AF,PF,true,Tag,RT>::data() noexcept
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dynamic matrix. Note that you
// can NOT assume that all matrix elements lie adjacent to each other! The dynamic matrix may
// use techniques such as padding to improve the alignment of the data. Whereas the number of
// elements within a column are given by the \c columns() member functions, the total number
// of elements including padding is given by the \c spacing() member function.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstPointer
   CustomMatrix<Type,AF,PF,true,Tag,RT>::data() const noexcept
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Pointer
   CustomMatrix<Type,AF,PF,true,Tag,RT>::data( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_+j*mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the matrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstPointer
   CustomMatrix<Type,AF,PF,true,Tag,RT>::data( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_+j*mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Iterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::begin( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_+j*mm_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::begin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_+j*mm_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of column \a j.
//
// \param j The column index.
// \return Iterator to the first element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::cbegin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_+j*mm_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::Iterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::end( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_+j*mm_+m_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::end( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_+j*mm_+m_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last element of column \a j.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomMatrix<Type,AF,PF,true,Tag,RT>::ConstIterator
   CustomMatrix<Type,AF,PF,true,Tag,RT>::cend( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_+j*mm_+m_ );
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
/*!\brief Homogenous assignment to all matrix elements.
//
// \param rhs Scalar value to be assigned to all matrix elements.
// \return Reference to the assigned matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( const Type& rhs )
{
   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         v_[i+j*mm_] = rhs;

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all matrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to static matrix.
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::rowMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   blaze::CustomMatrix<int,unaligned,unpadded,rowMajor> A( array, 3UL, 3UL );
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix elements are assigned the values from the given initializer list. Missing values
// are initialized as default (as e.g. the value 6 in the example). Note that in case the size
// of the top-level initializer list exceeds the number of rows or the size of any nested list
// exceeds the number of columns, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( initializer_list< initializer_list<Type> > list )
{
   using blaze::clear;

   if( list.size() != m_ || determineColumns( list ) > n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to custom matrix" );
   }

   size_t i( 0UL );

   for( const auto& rowList : list ) {
      size_t j( 0UL );
      for( const auto& element : rowList ) {
         v_[i+j*mm_] = element;
         ++j;
      }
      for( ; j<n_; ++j ) {
         clear( v_[i+j*mm_] );
      }
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array assignment to all matrix elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::columnMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::CustomMatrix<int,unaligned,unpadded,columnMajor> A( array, 3UL, 3UL );
   A = init;
   \endcode

// The matrix is assigned the values from the given static array. Missing values are initialized
// with default values (as e.g. the value 6 in the example). Note that the size of the static
// array must match the size of the custom matrix. Otherwise a \a std::invalid_argument exception
// is thrown. Also note that after the assignment \a array will have the same entries as \a init.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Rows       // Number of rows of the static array
        , size_t Cols >     // Number of columns of the static array
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( const Other (&array)[Rows][Cols] )
{
   if( m_ != Rows || n_ != Cols ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t j=0UL; j<Cols; ++j )
      for( size_t i=0UL; i<Rows; ++i )
         v_[i+j*mm_] = array[i][j];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array assignment to all matrix elements.
//
// \param array The given std::array for the assignment.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::unaligned;
   using blaze::unpadded;
   using blaze::columnMajor;

   const int array[9] = { 0, 0, 0,
                          0, 0, 0,
                          0, 0, 0 };
   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } } };
   blaze::CustomMatrix<int,unaligned,unpadded,columnMajor> A( array, 3UL, 3UL );
   A = init;
   \endcode

// The matrix is assigned the values from the given std::array. Missing values are initialized
// with default values (as e.g. the value 6 in the example). Note that the size of the std::array
// must match the size of the custom matrix. Otherwise a \a std::invalid_argument exception is
// thrown. Also note that after the assignment \a array will have the same entries as \a init.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Rows       // Number of rows of the static array
        , size_t Cols >     // Number of columns of the static array
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( const std::array<std::array<Other,Cols>,Rows>& array )
{
   if( m_ != Rows || n_ != Cols ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t j=0UL; j<Cols; ++j )
      for( size_t i=0UL; i<Rows; ++i )
         v_[i+j*mm_] = array[i][j];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for CustomMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The matrix is initialized as a copy of the given matrix. In case the current sizes of the two
// matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( const CustomMatrix& rhs )
{
   if( rhs.rows() != m_ || rhs.columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   smpAssign( *this, *rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Move assignment operator for CustomMatrix.
//
// \param rhs The matrix to be moved into this instance.
// \return Reference to the assigned matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( CustomMatrix&& rhs ) noexcept
{
   m_  = rhs.m_;
   mm_ = rhs.mm_;
   n_  = rhs.n_;
   v_  = rhs.v_;

   rhs.m_  = 0UL;
   rhs.mm_ = 0UL;
   rhs.n_  = 0UL;
   rhs.v_  = nullptr;

   BLAZE_INTERNAL_ASSERT( rhs.data() == nullptr, "Invalid data reference detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// The matrix is initialized as a copy of the given matrix. In case the current sizes of the two
// matrices don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator=( const Matrix<MT,SO>& rhs )
{
   using TT = decltype( trans( *this ) );
   using CT = decltype( ctrans( *this ) );
   using IT = decltype( inv( *this ) );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( IsSame_v<MT,TT> && (*rhs).isAliased( this ) ) {
      transpose();
   }
   else if( IsSame_v<MT,CT> && (*rhs).isAliased( this ) ) {
      ctranspose();
   }
   else if( !IsSame_v<MT,IT> && (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseMatrix_v<MT> )
         reset();
      smpAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator+=( const Matrix<MT,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, *rhs );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator-=( const Matrix<MT,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A\circ=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the matrix.
// \exception std::invalid_argument Matrix sizes do not match.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT    // Type of the right-hand side matrix
        , bool SO >      // Storage order of the right-hand side matrix
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::operator%=( const Matrix<MT,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( (*rhs).rows() != m_ || (*rhs).columns() != n_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<MT> tmp( *rhs );
      smpSchurAssign( *this, tmp );
   }
   else {
      smpSchurAssign( *this, *rhs );
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
/*!\brief Returns the current number of rows of the matrix.
//
// \return The number of rows of the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::rows() const noexcept
{
   return m_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current number of columns of the matrix.
//
// \return The number of columns of the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::columns() const noexcept
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two column, i.e. the total number
// of elements of a column.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::spacing() const noexcept
{
   return mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::capacity() const noexcept
{
   return mm_ * n_;
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   return mm_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the total number of non-zero elements in the matrix
//
// \return The number of non-zero elements in the dense matrix.
//
// This function returns the number of non-zero elements in the matrix (i.e. the elements that
// compare unequal to their default value). Note that the number of non-zero elements is always
// less than or equal to the total number of elements in the matrix.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::nonZeros() const
{
   size_t nonzeros( 0UL );

   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         if( !isDefault<strict>( v_[i+j*mm_] ) )
            ++nonzeros;

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
// This function returns the current number of non-zero elements in the specified column (i.e.
// the elements that compare unequal to their default value).
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomMatrix<Type,AF,PF,true,Tag,RT>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( j*mm_ + m_ );
   size_t nonzeros( 0UL );

   for( size_t i=j*mm_; i<iend; ++i )
      if( !isDefault<strict>( v_[i] ) )
         ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::reset()
{
   using blaze::clear;

   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         clear( v_[i+j*mm_] );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::reset( size_t j )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );
   for( size_t i=0UL; i<m_; ++i )
      clear( v_[i+j*mm_] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the \f$ M \times N \f$ matrix.
//
// \return void
//
// After the clear() function, the size of the matrix is 0.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::clear()
{
   m_  = 0UL;
   mm_ = 0UL;
   n_  = 0UL;
   v_  = nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Swapping the contents of two matrices.
//
// \param m The matrix to be swapped.
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::swap( CustomMatrix& m ) noexcept
{
   using std::swap;

   swap( m_ , m.m_  );
   swap( mm_, m.mm_ );
   swap( n_ , m.n_  );
   swap( v_ , m.v_  );
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
// \exception std::logic_error Impossible transpose operation.
//
// In case the matrix is not a square matrix, a \a std::logic_error exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>& CustomMatrix<Type,AF,PF,true,Tag,RT>::transpose()
{
   using std::swap;

   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Impossible transpose operation" );
   }

   for( size_t j=1UL; j<n_; ++j )
      for( size_t i=0UL; i<j; ++i )
         swap( v_[i+j*mm_], v_[j+i*mm_] );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the matrix.
//
// \return Reference to the transposed matrix.
// \exception std::logic_error Impossible transpose operation.
//
// In case the matrix is not a square matrix, a \a std::logic_error exception is thrown.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomMatrix<Type,AF,PF,true,Tag,RT>& CustomMatrix<Type,AF,PF,true,Tag,RT>::ctranspose()
{
   if( m_ != n_ ) {
      BLAZE_THROW_LOGIC_ERROR( "Impossible transpose operation" );
   }

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<j; ++i ) {
         cswap( v_[i+j*mm_], v_[j+i*mm_] );
      }
      conjugate( v_[j+j*mm_] );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the matrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the matrix scaling.
// \return Reference to the matrix.
//
// This function scales the matrix by applying the given scalar value \a scalar to each element
// of the matrix. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   CustomMatrix<int,unaligned,unpadded> A( ... );

   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the scalar value
inline CustomMatrix<Type,AF,PF,true,Tag,RT>&
   CustomMatrix<Type,AF,PF,true,Tag,RT>::scale( const Other& scalar )
{
   for( size_t j=0UL; j<n_; ++j )
      for( size_t i=0UL; i<m_; ++i )
         v_[i+j*mm_] *= scalar;

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RESOURCE MANAGEMENT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resets the custom matrix and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \return void
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This function resets the custom matrix to the given array of elements of size \f$ m \times n \f$.
// The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function is \b NOT available for padded custom matrices!
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom matrix referencing the array goes out of scope.
// \note The custom matrix does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::reset( Type* ptr, size_t m, size_t n )
{
   BLAZE_STATIC_ASSERT( PF == unpadded );

   CustomMatrix tmp( ptr, m, n );
   swap( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resets the custom matrix and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the matrix.
// \param m The number of rows of the array of elements.
// \param n The number of columns of the array of elements.
// \param mm The total number of elements between two columns.
// \return void
// \exception std::invalid_argument Invalid setup of custom matrix.
//
// This function resets the custom matrix to the given array of elements of size \f$ m \times n \f$.
// The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom matrix referencing the array goes out of scope.
// \note The custom matrix does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::reset( Type* ptr, size_t m, size_t n, size_t mm )
{
   CustomMatrix tmp( ptr, m, n, mm );
   swap( tmp );
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
// This function returns whether the given address can alias with the matrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomMatrix<Type,AF,PF,true,Tag,RT>::canAlias( const Other* alias ) const noexcept
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
// This function returns whether the given address is aliased with the matrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomMatrix<Type,AF,PF,true,Tag,RT>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the matrix is properly aligned in memory.
//
// \return \a true in case the matrix is aligned, \a false if not.
//
// This function returns whether the matrix is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of each column of the matrix are guaranteed to conform to
// the alignment restrictions of the element type \a Type.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomMatrix<Type,AF,PF,true,Tag,RT>::isAligned() const noexcept
{
   return ( AF || ( checkAlignment( v_ ) && rows() % SIMDSIZE == 0UL ) );
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomMatrix<Type,AF,PF,true,Tag,RT>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense matrix. The row index
// must be smaller than the number of rows and the column index must be smaller than the number
// of columns. Additionally, the row index must be a multiple of the number of values inside
// the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,true,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,true,Tag,RT>::load( size_t i, size_t j ) const noexcept
{
   if( AF && PF )
      return loada( i, j );
   else
      return loadu( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,true,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,true,Tag,RT>::loada( size_t i, size_t j ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= ( PF ? mm_ : m_ ), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !PF || i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i+j*mm_ ), "Invalid alignment detected" );

   return loada( v_+i+j*mm_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomMatrix<Type,AF,PF,true,Tag,RT>::SIMDType
   CustomMatrix<Type,AF,PF,true,Tag,RT>::loadu( size_t i, size_t j ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= ( PF ? mm_ : m_ ), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );

   return loadu( v_+i+j*mm_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense matrix. The row index
// must be smaller than the number of rows and the column index must be smaller then the number
// of columns. Additionally, the row index must be a multiple of the number of values inside the
// SIMD element. This function must \b NOT be called explicitly! It is used internally for the
// performance optimized evaluation of expression templates. Calling this function explicitly
// might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,true,Tag,RT>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( AF && PF )
      storea( i, j, value );
   else
      storeu( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,true,Tag,RT>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= ( PF ? mm_ : m_ ), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !PF || i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i+j*mm_ ), "Invalid alignment detected" );

   storea( v_+i+j*mm_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense matrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,true,Tag,RT>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= ( PF ? mm_ : m_ ), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );

   storeu( v_+i+j*mm_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the matrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense matrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the row index must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomMatrix<Type,AF,PF,true,Tag,RT>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= ( PF ? mm_ : m_ ), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !PF || i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+i+j*mm_ ), "Invalid alignment detected" );

   stream( v_+i+j*mm_, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::assign( const DenseMatrix<MT,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( m_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         v_[i    +j*mm_] = (*rhs)(i    ,j);
         v_[i+1UL+j*mm_] = (*rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         v_[ipos+j*mm_] = (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::assign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   const size_t ipos( remainder ? prevMultiple( m_, SIMDSIZE ) : m_ );
   BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

   if( AF && PF && useStreaming &&
       ( m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) ) && !(*rhs).isAliased( this ) )
   {
      for( size_t j=0UL; j<n_; ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT> right( (*rhs).begin(j) );

         for( ; i<ipos; i+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; remainder && i<m_; ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<n_; ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT> right( (*rhs).begin(j) );

         for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; remainder && i<m_; ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::assign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( min( n_, jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( min( m_, ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               v_[i+j*mm_] = (*rhs)(i,j);
            }
         }
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::assign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()+j*mm_] = element->value();
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
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::assign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).rows(); ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i+element->index()*mm_] = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::addAssign( const DenseMatrix<MT,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
   {
      if( IsDiagonal_v<MT> )
      {
         v_[j+j*mm_] += (*rhs)(j,j);
      }
      else
      {
         const size_t ibegin( ( IsLower_v<MT> )
                              ?( IsStrictlyLower_v<MT> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend  ( ( IsUpper_v<MT> )
                              ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                              :( m_ ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         size_t i( ibegin );

         for( ; (i+2UL) <= iend; i+=2UL ) {
            v_[i    +j*mm_] += (*rhs)(i    ,j);
            v_[i+1UL+j*mm_] += (*rhs)(i+1UL,j);
         }
         if( i < iend ) {
            v_[i+j*mm_] += (*rhs)(i,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::addAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t j=0UL; j<n_; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( m_ ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<iend; ++i ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::addAssign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( min( n_, jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block )
      {
         if( IsLower_v<MT> && ii < jj ) continue;
         if( IsUpper_v<MT> && ii > jj ) break;

         for( size_t j=jj; j<jend; ++j )
         {
            const size_t ibegin( ( IsLower_v<MT> )
                                 ?( max( ( IsStrictlyLower_v<MT> ? j+1UL : j ), ii ) )
                                 :( ii ) );
            const size_t iend  ( ( IsUpper_v<MT> )
                                 ?( min( ( IsStrictlyUpper_v<MT> ? j : j+1UL ), m_, ii+block ) )
                                 :( min( m_, ii+block ) ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            for( size_t i=ibegin; i<iend; ++i ) {
               v_[i+j*mm_] += (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::addAssign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()+j*mm_] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::addAssign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).rows(); ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i+element->index()*mm_] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::subAssign( const DenseMatrix<MT,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
   {
      if( IsDiagonal_v<MT> )
      {
         v_[j+j*mm_] -= (*rhs)(j,j);
      }
      else
      {
         const size_t ibegin( ( IsLower_v<MT> )
                              ?( IsStrictlyLower_v<MT> ? j+1UL : j )
                              :( 0UL ) );
         const size_t iend  ( ( IsUpper_v<MT> )
                              ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                              :( m_ ) );
         BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

         size_t i( ibegin );

         for( ; (i+2UL) <= iend; i+=2UL ) {
            v_[i  +j*mm_] -= (*rhs)(i  ,j);
            v_[i+1+j*mm_] -= (*rhs)(i+1,j);
         }
         if( i < iend ) {
            v_[i+j*mm_] -= (*rhs)(i,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a column-major
//        dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::subAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t j=0UL; j<n_; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( m_ ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( remainder ? prevMultiple( iend, SIMDSIZE ) : iend );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<iend; ++i ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::subAssign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( min( n_, jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block )
      {
         if( IsLower_v<MT> && ii < jj ) continue;
         if( IsUpper_v<MT> && ii > jj ) break;

         for( size_t j=jj; j<jend; ++j )
         {
            const size_t ibegin( ( IsLower_v<MT> )
                                 ?( max( ( IsStrictlyLower_v<MT> ? j+1UL : j ), ii ) )
                                 :( ii ) );
            const size_t iend  ( ( IsUpper_v<MT> )
                                 ?( min( ( IsStrictlyUpper_v<MT> ? j : j+1UL ), m_, ii+block ) )
                                 :( min( m_, ii+block ) ) );
            BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

            for( size_t i=ibegin; i<iend; ++i ) {
               v_[i+j*mm_] -= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::subAssign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()+j*mm_] -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::subAssign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).rows(); ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i+element->index()*mm_] -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::schurAssign( const DenseMatrix<MT,true>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( m_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         v_[i    +j*mm_] *= (*rhs)(i    ,j);
         v_[i+1UL+j*mm_] *= (*rhs)(i+1UL,j);
      }
      if( ipos < m_ ) {
         v_[ipos+j*mm_] *= (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a column-major
//        dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline auto CustomMatrix<Type,AF,PF,true,Tag,RT>::schurAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !PF || !IsPadded_v<MT> );

   for( size_t j=0UL; j<n_; ++j )
   {
      const size_t ipos( remainder ? prevMultiple( m_, SIMDSIZE ) : m_ );
      BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

      size_t i( 0UL );
      Iterator left( begin(j) );
      ConstIterator_t<MT> right( (*rhs).begin(j) );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; remainder && i<m_; ++i ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side dense matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::schurAssign( const DenseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<n_; jj+=block ) {
      const size_t jend( min( n_, jj+block ) );
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( min( m_, ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               v_[i+j*mm_] *= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::schurAssign( const SparseMatrix<MT,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( size_t j=0UL; j<(*rhs).columns(); ++j )
      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         v_[element->index()+j*mm_] = tmp(element->index(),j) * element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename MT >  // Type of the right-hand side sparse matrix
inline void CustomMatrix<Type,AF,PF,true,Tag,RT>::schurAssign( const SparseMatrix<MT,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( size_t i=0UL; i<(*rhs).rows(); ++i )
      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         v_[i+element->index()*mm_] = tmp(i,element->index()) * element->value();
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CUSTOMMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name CustomMatrix operators */
//@{
template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
void reset( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
void reset( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m, size_t i );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
void clear( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m );

template< RelaxationFlag RF, typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
bool isDefault( const CustomMatrix<Type,AF,PF,SO,Tag,RT>& m );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
bool isIntact( const CustomMatrix<Type,AF,PF,SO,Tag,RT>& m );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
void swap( CustomMatrix<Type,AF,PF,SO,Tag,RT>& a, CustomMatrix<Type,AF,PF,SO,Tag,RT>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given custom matrix.
// \ingroup custom_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void reset( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given custom matrix.
// \ingroup custom_matrix
//
// \param m The matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given custom matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void reset( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given custom matrix.
// \ingroup custom_matrix
//
// \param m The matrix to be cleared.
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void clear( CustomMatrix<Type,AF,PF,SO,Tag,RT>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given custom matrix is in default state.
// \ingroup custom_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the custom matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   using blaze::aligned;
   using blaze::padded;

   blaze::CustomMatrix<int,aligned,padded> A( ... );
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
        , AlignmentFlag AF   // Alignment flag
        , PaddingFlag PF     // Padding flag
        , bool SO            // Storage order
        , typename Tag       // Type tag
        , typename RT >      // Result type
inline bool isDefault( const CustomMatrix<Type,AF,PF,SO,Tag,RT>& m )
{
   return ( m.rows() == 0UL && m.columns() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given custom matrix are intact.
// \ingroup custom_matrix
//
// \param m The custom matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the custom matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::aligned;
   using blaze::padded;

   blaze::CustomMatrix<int,aligned,padded> A( ... );
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool isIntact( const CustomMatrix<Type,AF,PF,SO,Tag,RT>& m )
{
   return ( m.rows() * m.columns() <= m.capacity() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two custom matrices.
// \ingroup custom_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type     // Data type of the matrix
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool SO           // Storage order
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void swap( CustomMatrix<Type,AF,PF,SO,Tag,RT>& a, CustomMatrix<Type,AF,PF,SO,Tag,RT>& b ) noexcept
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
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
struct HasConstDataAccess< CustomMatrix<T,AF,PF,SO,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
struct HasMutableDataAccess< CustomMatrix<T,AF,PF,SO,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCUSTOM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
struct IsCustom< CustomMatrix<T,AF,PF,SO,Tag,RT> >
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
template< typename T, PaddingFlag PF, bool SO, typename Tag, typename RT >
struct IsAligned< CustomMatrix<T,aligned,PF,SO,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool SO, typename Tag, typename RT >
struct IsContiguous< CustomMatrix<T,AF,PF,SO,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************





//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, AlignmentFlag AF, bool SO, typename Tag, typename RT >
struct IsPadded< CustomMatrix<T,AF,padded,SO,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
