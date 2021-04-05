//=================================================================================================
/*!
//  \file blaze/math/dense/DynamicMatrix.h
//  \brief Header file for the implementation of a dynamic MxN matrix
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

#ifndef _BLAZE_MATH_DENSE_DYNAMICMATRIX_H_
#define _BLAZE_MATH_DENSE_DYNAMICMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <algorithm>
#include <memory>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/Diagonal.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
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
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/DynamicAllocator.h>
#include <blaze/math/typetraits/GetAllocator.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/NoUniqueAddress.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Restrict.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Algorithms.h>
#include <blaze/util/AlignedAllocator.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsSame.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/typetraits/RemoveConst.h>
#include <blaze/util/typetraits/RemoveCV.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup dynamic_matrix DynamicMatrix
// \ingroup dense_matrix
*/
/*!\brief Efficient implementation of a dynamic \f$ M \times N \f$ matrix.
// \ingroup dynamic_matrix
//
// The DynamicMatrix class template is the representation of an arbitrary sized matrix with
// \f$ M \times N \f$ dynamically allocated elements of arbitrary type. The type of the elements,
// the storage order, the type of the allocator, and the group tag of the matrix can be specified
// via the four template parameters:

   \code
   namespace blaze {

   template< typename Type, bool SO, typename Alloc, typename Tag >
   class DynamicMatrix;

   } // namespace blaze
   \endcode

//  - Type : specifies the type of the matrix elements. DynamicMatrix can be used with any
//           non-cv-qualified, non-reference, non-pointer element type.
//  - SO   : specifies the storage order (blaze::rowMajor, blaze::columnMajor) of the matrix.
//           The default value is blaze::defaultStorageOrder.
//  - Alloc: specifies the type of allocator used to allocate dynamic memory. The default type
//           of allocator is \c blaze::AlignedAllocator.
//  - Tag  : optional type parameter to tag the matrix. The default type is \c blaze::Group0.
//           See \ref grouping_tagging for details.
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

// The use of DynamicMatrix is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of row-major and
// column-major dense and sparse matrices with fitting element types. The following example gives
// an impression of the use of DynamicMatrix:

   \code
   using blaze::DynamicMatrix;
   using blaze::CompressedMatrix;
   using blaze::rowMajor;
   using blaze::columnMajor;

   DynamicMatrix<double,rowMajor> A( 2, 3 );  // Default constructed, non-initialized, row-major 2x3 matrix
   A(0,0) = 1.0; A(0,1) = 2.0; A(0,2) = 3.0;  // Initialization of the first row
   A(1,0) = 4.0; A(1,1) = 5.0; A(1,2) = 6.0;  // Initialization of the second row

   DynamicMatrix<float,columnMajor> B( 2, 3 );  // Default constructed column-major single precision 2x3 matrix
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
class DynamicMatrix
   : public DenseMatrix< DynamicMatrix<Type,SO,Alloc,Tag>, SO >
{
 private:
   //**********************************************************************************************
   //! Compilation switch for the choice of alignment.
   static constexpr AlignmentFlag align = ( usePadding ? aligned : unaligned );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This       = DynamicMatrix<Type,SO,Alloc,Tag>;  //!< Type of this DynamicMatrix instance.
   using BaseType   = DenseMatrix<This,SO>;              //!< Base type of this DynamicMatrix instance.
   using ResultType = This;                              //!< Result type for expression template evaluations.

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = DynamicMatrix<Type,!SO,Alloc,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = DynamicMatrix<Type,!SO,Alloc,Tag>;

   using ElementType   = Type;                      //!< Type of the matrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using AllocatorType = AlignedAllocator<Type>;    //!< Allocator type of this DynamicMatrix instance.
   using TagType       = Tag;                       //!< Tag type of this DynamicMatrix instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = Type*;        //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = DenseIterator<Type,align>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,align>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a DynamicMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind
   {
      //! The new type of allocator.
      using NewAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<NewType>;

      //! The type of the other DynamicMatrix.
      using Other = DynamicMatrix<NewType,SO,NewAlloc,Tag>;
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a DynamicMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize
   {
      using Other = DynamicMatrix<Type,SO,Alloc,Tag>;  //!< The type of the other DynamicMatrix.
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
   inline DynamicMatrix( const Alloc& alloc = Alloc{} ) noexcept;
   inline DynamicMatrix( size_t m, size_t n, const Alloc& alloc = Alloc{} );
   inline DynamicMatrix( size_t m, size_t n, const Type& init, const Alloc& alloc = Alloc{} );
   inline DynamicMatrix( initializer_list< initializer_list<Type> > list, const Alloc& alloc = Alloc{} );

   template< typename Other >
   inline DynamicMatrix( size_t m, size_t n, const Other* array, const Alloc& alloc = Alloc{} );

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix( const Other (&array)[Rows][Cols], const Alloc& alloc = Alloc{} );

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix( const std::array<std::array<Other,Cols>,Rows>& array, const Alloc& alloc = Alloc{} );

   inline DynamicMatrix( const DynamicMatrix& m );
   inline DynamicMatrix( DynamicMatrix&& m ) noexcept;

   template< typename MT, bool SO2 >
   inline DynamicMatrix( const Matrix<MT,SO2>& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~DynamicMatrix();
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
   inline DynamicMatrix& operator=( const Type& rhs ) &;
   inline DynamicMatrix& operator=( initializer_list< initializer_list<Type> > list ) &;

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix& operator=( const Other (&array)[Rows][Cols] ) &;

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix& operator=( const std::array<std::array<Other,Cols>,Rows>& array ) &;

   inline DynamicMatrix& operator=( const DynamicMatrix& rhs ) &;
   inline DynamicMatrix& operator=( DynamicMatrix&& rhs ) & noexcept;

   template< typename MT, bool SO2 > inline DynamicMatrix& operator= ( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline DynamicMatrix& operator+=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline DynamicMatrix& operator-=( const Matrix<MT,SO2>& rhs ) &;
   template< typename MT, bool SO2 > inline DynamicMatrix& operator%=( const Matrix<MT,SO2>& rhs ) &;
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
          void   resize ( size_t m, size_t n, bool preserve=true );
   inline void   extend ( size_t m, size_t n, bool preserve=true );
   inline void   reserve( size_t elements );
   inline void   shrinkToFit();
   inline void   swap( DynamicMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline DynamicMatrix& transpose();
   inline DynamicMatrix& ctranspose();

   template< typename Other > inline DynamicMatrix& scale( const Other& scalar );
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

   //**********************************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Debugging functions*************************************************************************
   /*!\name Debugging functions */
   //@{
   inline bool isIntact() const noexcept;
   //@}
   //**********************************************************************************************

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
   //**Uninitialized struct definition*************************************************************
   /*!\brief Definition of the nested auxiliary struct Uninitialized.
   */
   struct Uninitialized {};
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline DynamicMatrix( size_t m, size_t n, size_t nn, const Alloc& alloc, Uninitialized );
   inline DynamicMatrix( size_t m, size_t n, size_t nn, size_t capa, const Alloc& alloc, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t addPadding( size_t value ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;         //!< The current number of rows of the matrix.
   size_t n_;         //!< The current number of columns of the matrix.
   size_t nn_;        //!< The alignment adjusted number of columns.
   size_t capacity_;  //!< The maximum capacity of the matrix.

   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                             /*!< Access to the matrix elements is gained via the function call
                                  operator. In case of row-major order the memory layout of the
                                  elements is
                                  \f[\left(\begin{array}{*{5}{c}}
                                  0            & 1             & 2             & \cdots & N-1         \\
                                  N            & N+1           & N+2           & \cdots & 2 \cdot N-1 \\
                                  \vdots       & \vdots        & \vdots        & \ddots & \vdots      \\
                                  M \cdot N-N  & M \cdot N-N+1 & M \cdot N-N+2 & \cdots & M \cdot N-1 \\
                                  \end{array}\right)\f]. */

   BLAZE_NO_UNIQUE_ADDRESS Alloc alloc_;  //!< The allocator of the matrix.
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
//  DEDUCTION GUIDES
//
//=================================================================================================

//*************************************************************************************************
#if BLAZE_CPP17_MODE

template< typename Type >
DynamicMatrix( size_t, size_t, Type* ) -> DynamicMatrix< RemoveCV_t<Type> >;

template< typename Type, size_t M, size_t N >
DynamicMatrix( Type (&)[M][N] ) -> DynamicMatrix< RemoveCV_t<Type> >;

template< typename Type, size_t M, size_t N >
DynamicMatrix( std::array<std::array<Type,N>,M> ) -> DynamicMatrix<Type>;

#endif
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The (default) constructor for DynamicMatrix.
//
// \param alloc Allocator for all memory allocations of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( const Alloc& alloc ) noexcept
   : m_       ( 0UL )      // The current number of rows of the matrix
   , n_       ( 0UL )      // The current number of columns of the matrix
   , nn_      ( 0UL )      // The alignment adjusted number of columns
   , capacity_( 0UL )      // The maximum capacity of the matrix
   , v_       ( nullptr )  // The matrix elements
   , alloc_   ( alloc )    // The allocator of the matrix
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for DynamicMatrix.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nn The alignment adjusted number of columns.
// \param alloc Allocator for all memory allocations of this matrix.
// \exception std::bad_alloc Allocation failed.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, size_t nn, const Alloc& alloc, Uninitialized )
   : DynamicMatrix( m, n, nn, m*nn, alloc, Uninitialized{} )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for DynamicMatrix.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param nn The alignment adjusted number of columns.
// \param capa The initial capacity of the matrix.
// \param alloc Allocator for all memory allocations of this matrix.
// \exception std::bad_alloc Allocation failed.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, size_t nn, size_t capa, const Alloc& alloc, Uninitialized )
   : m_       ( m )        // The current number of rows of the matrix
   , n_       ( n )        // The current number of columns of the matrix
   , nn_      ( nn )       // The alignment adjusted number of columns
   , capacity_( capa )     // The maximum capacity of the matrix
   , v_       ( nullptr )  // The matrix elements
   , alloc_   ( alloc )    // The allocator of the matrix
{
   v_ = alloc_.allocate( capacity_ );

   if( !checkAlignment( v_ ) ) {
      alloc_.deallocate( v_, capacity_ );
      BLAZE_THROW_BAD_ALLOC;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. For built-in types no initialization
//        is performed!
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param alloc Allocator for all memory allocations of this matrix.
//
// \note This constructor is only responsible to allocate the required dynamic memory. For
// built-in types no initialization of the elements is performed!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Alloc& alloc )
   : DynamicMatrix( m, n, addPadding(n), alloc, Uninitialized{} )
{
   using blaze::clear;

   blaze::uninitialized_default_construct_n( v_, capacity_ );

   if( IsVectorizable_v<Type> && IsBuiltin_v<Type> ) {
      for( size_t i=0UL; i<m_; ++i ) {
         for( size_t j=n_; j<nn_; ++j ) {
            clear( v_[i*nn_+j] );
         }
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
// \param alloc Allocator for all memory allocations of this matrix.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Type& init, const Alloc& alloc )
   : DynamicMatrix( m, n, alloc )
{
   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n_; ++j ) {
         v_[i*nn_+j] = init;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all matrix elements.
//
// \param list The initializer list.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor provides the option to explicitly initialize the elements of the matrix by
// means of an initializer list:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<int,rowMajor> A{ { 1, 2, 3 },
                                         { 4, 5 },
                                         { 7, 8, 9 } };
   \endcode

// The matrix is sized according to the size of the initializer list and all its elements are
// (copy) assigned the values of the given initializer list. Missing values are initialized as
// default (as e.g. the value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( initializer_list< initializer_list<Type> > list, const Alloc& alloc )
   : DynamicMatrix( list.size(), determineColumns( list ), alloc )
{
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), begin(i) ), end(i), Type() );
      ++i;
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param array Dynamic array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a dynamic array:

   \code
   using blaze::rowMajor;

   int* array = new int[20];
   // ... Initialization of the dynamic array
   blaze::DynamicMatrix<int,rowMajor> v( 4UL, 5UL, array );
   delete[] array;
   \endcode

// The matrix is sized according to the given size of the array and initialized with the values
// from the given array. Note that it is expected that the given \a array has at least \a m by
// \a n elements. Providing an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the initialization array
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Other* array, const Alloc& alloc )
   : DynamicMatrix( m, n, alloc )
{
   for( size_t i=0UL; i<m; ++i ) {
      for( size_t j=0UL; j<n; ++j ) {
         v_[i*nn_+j] = array[i*n+j];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all matrix elements.
//
// \param array Static array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a static array:

   \code
   using blaze::rowMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::DynamicMatrix<int,rowMajor> A( init );
   \endcode

// The matrix is sized according to the size of the static array and initialized with the values
// from the given static array. Missing values are initialized with default values (as e.g. the
// value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the static array
        , size_t Rows     // Number of rows of the static array
        , size_t Cols >   // Number of columns of the static array
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( const Other (&array)[Rows][Cols], const Alloc& alloc )
   : DynamicMatrix( Rows, Cols, alloc )
{
   for( size_t i=0UL; i<Rows; ++i ) {
      for( size_t j=0UL; j<Cols; ++j ) {
         v_[i*nn_+j] = array[i][j];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of all matrix elements from the given std::array.
//
// \param array The given std::array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a std::array:

   \code
   using blaze::rowMajor;

   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } } };
   blaze::DynamicMatrix<int,rowMajor> A( init );
   \endcode

// The matrix is sized according to the size of the std::array and initialized with the values
// from the given std::array. Missing values are initialized with default values (as e.g. the
// value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the std::array
        , size_t Rows     // Number of rows of the std::array
        , size_t Cols >   // Number of columns of the std::array
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( const std::array<std::array<Other,Cols>,Rows>& array, const Alloc& alloc )
   : DynamicMatrix( Rows, Cols, alloc )
{
   for( size_t i=0UL; i<Rows; ++i ) {
      for( size_t j=0UL; j<Cols; ++j ) {
         v_[i*nn_+j] = array[i][j];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for DynamicMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( const DynamicMatrix& m )
   : DynamicMatrix( m.m_, m.n_ )
{
   BLAZE_INTERNAL_ASSERT( capacity_ <= m.capacity_, "Invalid capacity estimation" );

   smpAssign( *this, m );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for DynamicMatrix.
//
// \param m The matrix to be move into this instance.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( DynamicMatrix&& m ) noexcept
   : m_       ( m.m_        )  // The current number of rows of the matrix
   , n_       ( m.n_        )  // The current number of columns of the matrix
   , nn_      ( m.nn_       )  // The alignment adjusted number of columns
   , capacity_( m.capacity_ )  // The maximum capacity of the matrix
   , v_       ( m.v_        )  // The matrix elements
{
   m.m_        = 0UL;
   m.n_        = 0UL;
   m.nn_       = 0UL;
   m.capacity_ = 0UL;
   m.v_        = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign matrix
        , bool SO2 >      // Storage order of the foreign matrix
inline DynamicMatrix<Type,SO,Alloc,Tag>::DynamicMatrix( const Matrix<MT,SO2>& m )
   : DynamicMatrix( (*m).rows(), (*m).columns() )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( IsSparseMatrix_v<MT> && IsBuiltin_v<Type> ) {
      reset();
   }

   smpAssign( *this, *m );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for DynamicMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>::~DynamicMatrix()
{
   blaze::destroy_n( v_, capacity_ );
   alloc_.deallocate( v_, capacity_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Reference
   DynamicMatrix<Type,SO,Alloc,Tag>::operator()( size_t i, size_t j ) noexcept
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstReference
   DynamicMatrix<Type,SO,Alloc,Tag>::operator()( size_t i, size_t j ) const noexcept
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Reference
   DynamicMatrix<Type,SO,Alloc,Tag>::at( size_t i, size_t j )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstReference
   DynamicMatrix<Type,SO,Alloc,Tag>::at( size_t i, size_t j ) const
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
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Pointer
   DynamicMatrix<Type,SO,Alloc,Tag>::data() noexcept
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstPointer
   DynamicMatrix<Type,SO,Alloc,Tag>::data() const noexcept
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Pointer
   DynamicMatrix<Type,SO,Alloc,Tag>::data( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_ + i*nn_;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstPointer
   DynamicMatrix<Type,SO,Alloc,Tag>::data( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return v_ + i*nn_;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Iterator
   DynamicMatrix<Type,SO,Alloc,Tag>::begin( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*nn_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,SO,Alloc,Tag>::begin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,SO,Alloc,Tag>::cbegin( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::Iterator
   DynamicMatrix<Type,SO,Alloc,Tag>::end( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return Iterator( v_ + i*nn_ + n_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,SO,Alloc,Tag>::end( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ + n_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,SO,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,SO,Alloc,Tag>::cend( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < m_, "Invalid dense matrix row access index" );
   return ConstIterator( v_ + i*nn_ + n_ );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( const Type& rhs ) &
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
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<int,rowMajor> A;
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix is resized according to the given initializer list and all its elements are
// (copy) assigned the values from the given initializer list. Missing values are initialized
// as default (as e.g. the value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( initializer_list< initializer_list<Type> > list ) &
{
   resize( list.size(), determineColumns( list ), false );

   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), v_+i*nn_ ), v_+(i+1UL)*nn_, Type() );
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
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::rowMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::DynamicMatrix<int,rowMajor> A;
   A = init;
   \endcode

// The matrix is resized according to the size of the static array and assigned the values of the
// given static array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the static array
        , size_t Rows     // Number of rows of the static array
        , size_t Cols >   // Number of columns of the static array
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( const Other (&array)[Rows][Cols] ) &
{
   resize( Rows, Cols, false );

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
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::rowMajor;

   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } } };
   blaze::DynamicMatrix<int,rowMajor> A;
   A = init;
   \endcode

// The matrix is resized according to the size of the std::array and assigned the values of the
// given std::array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the std::array
        , size_t Rows     // Number of rows of the std::array
        , size_t Cols >   // Number of columns of the std::array
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( const std::array<std::array<Other,Cols>,Rows>& array ) &
{
   resize( Rows, Cols, false );

   for( size_t i=0UL; i<Rows; ++i )
      for( size_t j=0UL; j<Cols; ++j )
         v_[i*nn_+j] = array[i][j];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for DynamicMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( const DynamicMatrix& rhs ) &
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );
   smpAssign( *this, *rhs );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for DynamicMatrix.
//
// \param rhs The matrix to be moved into this instance.
// \return Reference to the assigned matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( DynamicMatrix&& rhs ) & noexcept
{
   blaze::destroy_n( v_, capacity_ );
   alloc_.deallocate( v_, capacity_ );

   m_        = rhs.m_;
   n_        = rhs.n_;
   nn_       = rhs.nn_;
   capacity_ = rhs.capacity_;
   v_        = rhs.v_;

   rhs.m_        = 0UL;
   rhs.n_        = 0UL;
   rhs.nn_       = 0UL;
   rhs.capacity_ = 0UL;
   rhs.v_        = nullptr;

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator=( const Matrix<MT,SO2>& rhs ) &
{
   using TT = decltype( trans( *this ) );
   using CT = decltype( ctrans( *this ) );
   using IT = decltype( inv( *this ) );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( IsSame_v<MT,TT> && (*rhs).isAliased( this ) ) {
      transpose();
   }
   else if( IsSame_v<MT,CT> && (*rhs).isAliased( this ) ) {
      ctranspose();
   }
   else if( !IsSame_v<MT,IT> && (*rhs).canAlias( this ) ) {
      DynamicMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      if( IsSparseMatrix_v<MT> )
         reset();
      smpAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator+=( const Matrix<MT,SO2>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator-=( const Matrix<MT,SO2>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO2 >      // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::operator%=( const Matrix<MT,SO2>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::rows() const noexcept
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::columns() const noexcept
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::spacing() const noexcept
{
   return nn_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the matrix.
//
// \return The capacity of the matrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::capacity() const noexcept
{
   return capacity_;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::capacity( size_t i ) const noexcept
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::nonZeros() const
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::nonZeros( size_t i ) const
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::reset()
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::reset( size_t i )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::clear()
{
   resize( 0UL, 0UL, false );
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
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Note that this function may invalidate all existing views (submatrices, rows, columns,
// ...) on the matrix if it is used to shrink the matrix. Additionally, the resize operation
// potentially changes all matrix elements. In order to preserve the old matrix values, the
// \a preserve flag can be set to \a true. However, new matrix elements of built-in type are
// not initialized!
//
// The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
// \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
void DynamicMatrix<Type,SO,Alloc,Tag>::resize( size_t m, size_t n, bool preserve )
{
   using blaze::clear;
   using blaze::min;

   if( m == m_ && n == n_ ) return;

   const size_t nn( addPadding( n ) );

   if( preserve )
   {
      const size_t min_m( min( m, m_ ) );
      const size_t min_n( min( n, n_ ) );

      DynamicMatrix tmp( m, n, nn, Alloc{}, Uninitialized{} );

      for( size_t i=0UL; i<min_m; ++i ) {
         blaze::uninitialized_transfer( v_+i*nn_, v_+i*nn_+min_n, tmp.v_+i*nn );
         blaze::uninitialized_default_construct( tmp.v_+i*nn+min_n, tmp.v_+i*nn+nn );
      }
      blaze::uninitialized_default_construct( tmp.v_+min_m*nn, tmp.v_+m*nn );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }
   else if( m*nn > capacity_ )
   {
      DynamicMatrix tmp( m, n, nn, Alloc{}, Uninitialized{} );

      blaze::uninitialized_default_construct( tmp.v_, tmp.v_+tmp.capacity_ );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }

   if( IsVectorizable_v<Type> ) {
      for( size_t i=0UL; i<m; ++i )
         for( size_t j=n; j<nn; ++j )
            clear( v_[i*nn+j] );
   }

   m_  = m;
   n_  = n;
   nn_ = nn;

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
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
// This function increases the matrix size by \a m rows and \a n columns. During this operation,
// new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
// this function potentially changes all matrix elements. In order to preserve the old matrix
// values, the \a preserve flag can be set to \a true. However, new matrix elements of built-in
// type are not initialized!
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the matrix.
//
// \param elements The new minimum capacity of the dense matrix.
// \return void
//
// This function increases the capacity of the dense matrix to at least \a elements elements.
// The current values of the matrix elements are preserved.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::reserve( size_t elements )
{
   using blaze::clear;

   if( elements > capacity_ )
   {
      DynamicMatrix tmp( m_, n_, nn_, elements, Alloc{}, Uninitialized{} );

      blaze::uninitialized_transfer( v_, v_+capacity_, tmp.v_ );
      blaze::uninitialized_value_construct( tmp.v_+capacity_, tmp.v_+elements );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the matrix by removing unused capacity. Please note
// that due to padding the capacity might not be reduced exactly to rows() times columns().
// Please also note that in case a reallocation occurs, all iterators (including end() iterators),
// all pointers and references to elements of this matrix are invalidated.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::shrinkToFit()
{
   if( ( m_ * nn_ ) < capacity_ ) {
      DynamicMatrix( *this ).swap( *this );
   }
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,SO,Alloc,Tag>::swap( DynamicMatrix& m ) noexcept
{
   using std::swap;

   swap( m_ , m.m_  );
   swap( n_ , m.n_  );
   swap( nn_, m.nn_ );
   swap( capacity_, m.capacity_ );
   swap( v_ , m.v_  );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Add the necessary amount of padding to the given value.
//
// \param value The value to be padded.
// \return The padded value.
//
// This function increments the given \a value by the necessary amount of padding based on the
// vector's data type \a Type.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,SO,Alloc,Tag>::addPadding( size_t value ) const noexcept
{
   if( usePadding && IsVectorizable_v<Type> )
      return nextMultiple<size_t>( value, SIMDSIZE );
   else return value;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>& DynamicMatrix<Type,SO,Alloc,Tag>::transpose()
{
   using std::swap;

   constexpr size_t block( BLOCK_SIZE );

   if( m_ == n_ )
   {
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( min( ii+block, m_ ) );
         for( size_t jj=0UL; jj<=ii; jj+=block ) {
            for( size_t i=ii; i<iend; ++i ) {
               const size_t jend( min( jj+block, n_, i ) );
               for( size_t j=jj; j<jend; ++j ) {
                  swap( v_[i*nn_+j], v_[j*nn_+i] );
               }
            }
         }
      }
   }
   else
   {
      DynamicMatrix tmp( trans(*this) );
      this->swap( tmp );
   }

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,SO,Alloc,Tag>& DynamicMatrix<Type,SO,Alloc,Tag>::ctranspose()
{
   constexpr size_t block( BLOCK_SIZE );

   if( m_ == n_ )
   {
      for( size_t ii=0UL; ii<m_; ii+=block ) {
         const size_t iend( min( ii+block, m_ ) );
         for( size_t jj=0UL; jj<ii; jj+=block ) {
            const size_t jend( min( jj+block, n_ ) );
            for( size_t i=ii; i<iend; ++i ) {
               for( size_t j=jj; j<jend; ++j ) {
                  cswap( v_[i*nn_+j], v_[j*nn_+i] );
               }
            }
         }
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=ii; j<i; ++j ) {
               cswap( v_[i*nn_+j], v_[j*nn_+i] );
            }
            conjugate( v_[i*nn_+i] );
         }
      }
   }
   else
   {
      DynamicMatrix tmp( ctrans(*this) );
      swap( tmp );
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
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , bool SO           // Storage order
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline DynamicMatrix<Type,SO,Alloc,Tag>&
   DynamicMatrix<Type,SO,Alloc,Tag>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<m_; ++i )
      for( size_t j=0UL; j<n_; ++j )
         v_[i*nn_+j] *= scalar;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  DEBUGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the invariants of the dynamic matrix are intact.
//
// \return \a true in case the dynamic matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dynamic matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,SO,Alloc,Tag>::isIntact() const noexcept
{
   if( n_ > nn_ || m_ * nn_ > capacity_ )
      return false;

   if( IsVectorizable_v<Type> ) {
      for( size_t i=0UL; i<m_; ++i ) {
         for( size_t j=n_; j<nn_; ++j ) {
            if( !isDefault<strict>( v_[i*nn_+j] ) )
               return false;
         }
      }
   }

   return true;
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
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,SO,Alloc,Tag>::canAlias( const Other* alias ) const noexcept
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
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,SO,Alloc,Tag>::isAliased( const Other* alias ) const noexcept
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,SO,Alloc,Tag>::isAligned() const noexcept
{
   return ( usePadding || columns() % SIMDSIZE == 0UL );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,SO,Alloc,Tag>::canSMPAssign() const noexcept
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,SO,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,SO,Alloc,Tag>::load( size_t i, size_t j ) const noexcept
{
   if( usePadding )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,SO,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,SO,Alloc,Tag>::loada( size_t i, size_t j ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || j % SIMDSIZE == 0UL, "Invalid column access index" );
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,SO,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,SO,Alloc,Tag>::loadu( size_t i, size_t j ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= nn_, "Invalid column access index" );

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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,SO,Alloc,Tag>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( usePadding )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,SO,Alloc,Tag>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || j % SIMDSIZE == 0UL, "Invalid column access index" );
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,SO,Alloc,Tag>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= nn_, "Invalid column access index" );

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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,SO,Alloc,Tag>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < n_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= nn_, "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || j % SIMDSIZE == 0UL, "Invalid column access index" );
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::assign( const DenseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::assign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

   const size_t jpos( remainder ? prevMultiple( n_, SIMDSIZE ) : n_ );
   BLAZE_INTERNAL_ASSERT( jpos <= n_, "Invalid end calculation" );

   if( usePadding && useStreaming &&
       ( m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<m_; ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT> right( (*rhs).begin(i) );

         for( ; j<jpos; j+=SIMDSIZE, left+=SIMDSIZE, right+=SIMDSIZE ) {
            left.stream( right.load() );
         }
         for( ; remainder && j<n_; ++j, ++left, ++right ) {
            *left = *right;
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::assign( const DenseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::assign( const SparseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::assign( const SparseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::addAssign( const DenseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::addAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::addAssign( const DenseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::addAssign( const SparseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::addAssign( const SparseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::subAssign( const DenseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::subAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::subAssign( const DenseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::subAssign( const SparseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::subAssign( const SparseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::schurAssign( const DenseMatrix<MT,SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,SO,Alloc,Tag>::schurAssign( const DenseMatrix<MT,SO>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

   for( size_t i=0UL; i<m_; ++i )
   {
      const size_t jpos( remainder ? prevMultiple( n_, SIMDSIZE ) : n_ );
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::schurAssign( const DenseMatrix<MT,!SO>& rhs )
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::schurAssign( const SparseMatrix<MT,SO>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
   {
      size_t j( 0UL );

      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( v_[i*nn_+j] );
         v_[i*nn_+j] *= element->value();
         ++j;
      }

      for( ; j<n_; ++j ) {
         reset( v_[i*nn_+j] );
      }
   }
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
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,SO,Alloc,Tag>::schurAssign( const SparseMatrix<MT,!SO>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
   {
      size_t i( 0UL );

      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( v_[i*nn_+j] );
         v_[i*nn_+j] *= element->value();
         ++i;
      }

      for( ; i<m_; ++i ) {
         reset( v_[i*nn_+j] );
      }
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
/*!\brief Specialization of DynamicMatrix for column-major matrices.
// \ingroup dynamic_matrix
//
// This specialization of DynamicMatrix adapts the class template to the requirements of
// column-major matrices.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
class DynamicMatrix<Type,true,Alloc,Tag>
   : public DenseMatrix< DynamicMatrix<Type,true,Alloc,Tag>, true >
{
 private:
   //**********************************************************************************************
   //! Compilation switch for the choice of alignment.
   static constexpr AlignmentFlag align = ( usePadding ? aligned : unaligned );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This       = DynamicMatrix<Type,true,Alloc,Tag>;  //!< Type of this DynamicMatrix instance.
   using BaseType   = DenseMatrix<This,true>;              //!< Base type of this DynamicMatrix instance.
   using ResultType = This;                                //!< Result type for expression template evaluations.

   //! Result type with opposite storage order for expression template evaluations.
   using OppositeType = DynamicMatrix<Type,false,Alloc,Tag>;

   //! Transpose type for expression template evaluations.
   using TransposeType = DynamicMatrix<Type,false,Alloc,Tag>;

   using ElementType   = Type;                      //!< Type of the matrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the matrix elements.
   using AllocatorType = AlignedAllocator<Type>;    //!< Allocator type of this DynamicMatrix instance.
   using TagType       = Tag;                       //!< Tag type of this DynamicMatrix instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using CompositeType = const This&;               //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant matrix value.
   using ConstReference = const Type&;  //!< Reference to a constant matrix value.
   using Pointer        = Type*;        //!< Pointer to a non-constant matrix value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant matrix value.

   using Iterator      = DenseIterator<Type,align>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,align>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a DynamicMatrix with different data/element type.
   */
   template< typename NewType >  // Data type of the other matrix
   struct Rebind
   {
      //! The new type of allocator.
      using NewAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<NewType>;

      //! The type of the other DynamicMatrix.
      using Other = DynamicMatrix<NewType,true,NewAlloc,Tag>;
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a DynamicMatrix with different fixed dimensions.
   */
   template< size_t NewM    // Number of rows of the other matrix
           , size_t NewN >  // Number of columns of the other matrix
   struct Resize
   {
      using Other = DynamicMatrix<Type,true,Alloc,Tag>;  //!< The type of the other DynamicMatrix.
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
   inline DynamicMatrix( const Alloc& alloc = Alloc{} ) noexcept;
   inline DynamicMatrix( size_t m, size_t n, const Alloc& alloc = Alloc{} );
   inline DynamicMatrix( size_t m, size_t n, const Type& init, const Alloc& alloc = Alloc{} );
   inline DynamicMatrix( initializer_list< initializer_list<Type> > list, const Alloc& alloc = Alloc{} );

   template< typename Other >
   inline DynamicMatrix( size_t m, size_t n, const Other* array, const Alloc& alloc = Alloc{} );

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix( const Other (&array)[Rows][Cols], const Alloc& alloc = Alloc{} );

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix( const std::array<std::array<Other,Cols>,Rows>& array, const Alloc& alloc = Alloc{} );

   inline DynamicMatrix( const DynamicMatrix& m );
   inline DynamicMatrix( DynamicMatrix&& m );

   template< typename MT, bool SO >
   inline DynamicMatrix( const Matrix<MT,SO>& m );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~DynamicMatrix();
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
   inline DynamicMatrix& operator=( const Type& rhs ) &;
   inline DynamicMatrix& operator=( initializer_list< initializer_list<Type> > list ) &;

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix& operator=( const Other (&array)[Rows][Cols] ) &;

   template< typename Other, size_t Rows, size_t Cols >
   inline DynamicMatrix& operator=( const std::array<std::array<Other,Cols>,Rows>& array ) &;

   inline DynamicMatrix& operator=( const DynamicMatrix& rhs ) &;
   inline DynamicMatrix& operator=( DynamicMatrix&& rhs ) &;

   template< typename MT, bool SO > inline DynamicMatrix& operator= ( const Matrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline DynamicMatrix& operator+=( const Matrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline DynamicMatrix& operator-=( const Matrix<MT,SO>& rhs ) &;
   template< typename MT, bool SO > inline DynamicMatrix& operator%=( const Matrix<MT,SO>& rhs ) &;
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
          void   resize ( size_t m, size_t n, bool preserve=true );
   inline void   extend ( size_t m, size_t n, bool preserve=true );
   inline void   reserve( size_t elements );
   inline void   shrinkToFit();
   inline void   swap( DynamicMatrix& m ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline DynamicMatrix& transpose();
   inline DynamicMatrix& ctranspose();

   template< typename Other > inline DynamicMatrix& scale( const Other& scalar );
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

   //**********************************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Debugging functions*************************************************************************
   /*!\name Debugging functions */
   //@{
   inline bool isIntact() const noexcept;
   //@}
   //**********************************************************************************************

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
   //**Uninitialized struct definition*************************************************************
   /*!\brief Definition of the nested auxiliary struct Uninitialized.
   */
   struct Uninitialized {};
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline DynamicMatrix( size_t m, size_t mm, size_t n, const Alloc& alloc, Uninitialized );
   inline DynamicMatrix( size_t m, size_t mm, size_t n, size_t capa, const Alloc& alloc, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t addPadding( size_t minRows ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t m_;         //!< The current number of rows of the matrix.
   size_t mm_;        //!< The alignment adjusted number of rows.
   size_t n_;         //!< The current number of columns of the matrix.
   size_t capacity_;  //!< The maximum capacity of the matrix.

   Type* BLAZE_RESTRICT v_;  //!< The dynamically allocated matrix elements.
                             /*!< Access to the matrix elements is gained via the function
                                  call operator. */

   BLAZE_NO_UNIQUE_ADDRESS Alloc alloc_;  //!< The allocator of the matrix.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
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
/*!\brief The (default) constructor for DynamicMatrix.
//
// \param alloc Allocator for all memory allocations of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( const Alloc& alloc ) noexcept
   : m_       ( 0UL )      // The current number of rows of the matrix
   , mm_      ( 0UL )      // The alignment adjusted number of rows
   , n_       ( 0UL )      // The current number of columns of the matrix
   , capacity_( 0UL )      // The maximum capacity of the matrix
   , v_       ( nullptr )  // The matrix elements
   , alloc_   ( alloc )    // The allocator of the matrix
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for DynamicMatrix.
//
// \param m The number of rows of the matrix.
// \param mm The alignment adjusted number of rows.
// \param n The number of columns of the matrix.
// \param alloc Allocator for all memory allocations of this matrix.
// \exception std::bad_alloc Allocation failed.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( size_t m, size_t mm, size_t n, const Alloc& alloc, Uninitialized )
   : DynamicMatrix( m, mm, n, mm*n, alloc, Uninitialized{} )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for DynamicMatrix.
//
// \param m The number of rows of the matrix.
// \param mm The alignment adjusted number of rows.
// \param n The number of columns of the matrix.
// \param capa The initial capacity of the matrix.
// \param alloc Allocator for all memory allocations of this matrix.
// \exception std::bad_alloc Allocation failed.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( size_t m, size_t mm, size_t n, size_t capa, const Alloc& alloc, Uninitialized )
   : m_       ( m )        // The current number of rows of the matrix
   , mm_      ( mm )       // The alignment adjusted number of rows
   , n_       ( n )        // The current number of columns of the matrix
   , capacity_( capa )     // The maximum capacity of the matrix
   , v_       ( nullptr )  // The matrix elements
   , alloc_   ( alloc )    // The allocator of the matrix
{
   v_ = alloc_.allocate( capacity_ );

   if( !checkAlignment( v_ ) ) {
      alloc_.deallocate( v_, capacity_ );
      BLAZE_THROW_BAD_ALLOC;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a matrix of size \f$ m \times n \f$. For built-in types no initialization
//        is performed!
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param alloc Allocator for all memory allocations of this matrix.
//
// \note This constructor is only responsible to allocate the required dynamic memory. For
// built-in types no initialization of the elements is performed!
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Alloc& alloc )
   : DynamicMatrix( m, addPadding(m), n, alloc, Uninitialized{} )
{
   using blaze::clear;

   blaze::uninitialized_default_construct_n( v_, capacity_ );

   if( IsVectorizable_v<Type> && IsBuiltin_v<Type> ) {
      for( size_t j=0UL; j<n_; ++j ) {
         for( size_t i=m_; i<mm_; ++i ) {
            clear( v_[i+j*mm_] );
         }
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a homogenous initialization of all \f$ m \times n \f$ matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param init The initial value of the matrix elements.
// \param alloc Allocator for all memory allocations of this matrix.
//
// All matrix elements are initialized with the specified value.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Type& init, const Alloc& alloc )
   : DynamicMatrix( m, n, alloc )
{
   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; i<m_; ++i ) {
         v_[i+j*mm_] = init;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List initialization of all matrix elements.
//
// \param list The initializer list.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor provides the option to explicitly initialize the elements of the matrix by
// means of an initializer list:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<int,columnMajor> A{ { 1, 2, 3 },
                                            { 4, 5 },
                                            { 7, 8, 9 } };
   \endcode

// The matrix is sized according to the size of the initializer list and all its elements are
// (copy) assigned the values of the given initializer list. Missing values are initialized as
// default (as e.g. the value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( initializer_list< initializer_list<Type> > list, const Alloc& alloc )
   : DynamicMatrix( list.size(), determineColumns( list ), alloc )
{
   using blaze::clear;

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

   BLAZE_INTERNAL_ASSERT( i == m_, "Invalid number of elements detected" );
   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array initialization of all matrix elements.
//
// \param m The number of rows of the matrix.
// \param n The number of columns of the matrix.
// \param array Dynamic array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a dynamic array:

   \code
   using blaze::columnMajor;

   int* array = new int[20];
   // ... Initialization of the dynamic array
   blaze::DynamicMatrix<int,columnMajor> v( array, 5UL, 4UL );
   delete[] array;
   \endcode

// The matrix is sized according to the given size of the array and initialized with the values
// from the given array. Note that it is expected that the given \a array has at least \a m by
// \a n elements. Providing an array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the matrix
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the initialization array
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( size_t m, size_t n, const Other* array, const Alloc& alloc )
   : DynamicMatrix( m, n, alloc )
{
   for( size_t j=0UL; j<n; ++j ) {
      for( size_t i=0UL; i<m; ++i ) {
         v_[i+j*mm_] = array[i+j*m];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array initialization of all matrix elements.
//
// \param array Static array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a static array:

   \code
   using blaze::columnMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::DynamicMatrix<int,columnMajor> A( init );
   \endcode

// The matrix is sized according to the size of the static array and initialized with the values
// from the given static array. Missing values are initialized with default values (as e.g. the
// value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the static array
        , size_t Rows     // Number of rows of the static array
        , size_t Cols >   // Number of columns of the static array
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( const Other (&array)[Rows][Cols], const Alloc& alloc )
   : DynamicMatrix( Rows, Cols, alloc )
{
   for( size_t j=0UL; j<Cols; ++j ) {
      for( size_t i=0UL; i<Rows; ++i ) {
         v_[i+j*mm_] = array[i][j];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Initialization of all matrix elements from the given std::array.
//
// \param array The given std::array for the initialization.
// \param alloc Allocator for all memory allocations of this matrix.
//
// This constructor offers the option to directly initialize the elements of the matrix with
// a std::array:

   \code
   using blaze::columnMajor;

   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } } };
   blaze::DynamicMatrix<int,columnMajor> A( init );
   \endcode

// The matrix is sized according to the size of the std::array and initialized with the values
// from the given std::array. Missing values are initialized with default values (as e.g. the
// value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the std::array
        , size_t Rows     // Number of rows of the std::array
        , size_t Cols >   // Number of columns of the std::array
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( const std::array<std::array<Other,Cols>,Rows>& array, const Alloc& alloc )
   : DynamicMatrix( Rows, Cols, alloc )
{
   for( size_t j=0UL; j<Cols; ++j ) {
      for( size_t i=0UL; i<Rows; ++i ) {
         v_[i+j*mm_] = array[i][j];
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for DynamicMatrix.
//
// \param m Matrix to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( const DynamicMatrix& m )
   : DynamicMatrix( m.m_, m.n_ )
{
   BLAZE_INTERNAL_ASSERT( capacity_ <= m.capacity_, "Invalid capacity estimation" );

   smpAssign( *this, m );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The move constructor for DynamicMatrix.
//
// \param m The matrix to be moved into this instance.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( DynamicMatrix&& m )
   : m_       ( m.m_        )  // The current number of rows of the matrix
   , mm_      ( m.mm_       )  // The alignment adjusted number of rows
   , n_       ( m.n_        )  // The current number of columns of the matrix
   , capacity_( m.capacity_ )  // The maximum capacity of the matrix
   , v_       ( m.v_        )  // The matrix elements
{
   m.m_        = 0UL;
   m.mm_       = 0UL;
   m.n_        = 0UL;
   m.capacity_ = 0UL;
   m.v_        = nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Conversion constructor from different matrices.
//
// \param m Matrix to be copied.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the foreign matrix
        , bool SO >       // Storage order of the foreign matrix
inline DynamicMatrix<Type,true,Alloc,Tag>::DynamicMatrix( const Matrix<MT,SO>& m )
   : DynamicMatrix( (*m).rows(), (*m).columns() )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( IsSparseMatrix_v<MT> && IsBuiltin_v<Type> ) {
      reset();
   }

   smpAssign( *this, *m );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
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
/*!\brief The destructor for DynamicMatrix.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>::~DynamicMatrix()
{
   blaze::destroy_n( v_, capacity_ );
   alloc_.deallocate( v_, capacity_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Reference
   DynamicMatrix<Type,true,Alloc,Tag>::operator()( size_t i, size_t j ) noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstReference
   DynamicMatrix<Type,true,Alloc,Tag>::operator()( size_t i, size_t j ) const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Reference
   DynamicMatrix<Type,true,Alloc,Tag>::at( size_t i, size_t j )
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstReference
   DynamicMatrix<Type,true,Alloc,Tag>::at( size_t i, size_t j ) const
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Pointer
   DynamicMatrix<Type,true,Alloc,Tag>::data() noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstPointer
   DynamicMatrix<Type,true,Alloc,Tag>::data() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Pointer
   DynamicMatrix<Type,true,Alloc,Tag>::data( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_ + j*mm_;
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstPointer
   DynamicMatrix<Type,true,Alloc,Tag>::data( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return v_ + j*mm_;
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Iterator
   DynamicMatrix<Type,true,Alloc,Tag>::begin( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*mm_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,true,Alloc,Tag>::begin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,true,Alloc,Tag>::cbegin( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::Iterator
   DynamicMatrix<Type,true,Alloc,Tag>::end( size_t j ) noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return Iterator( v_ + j*mm_ + m_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,true,Alloc,Tag>::end( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ + m_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline typename DynamicMatrix<Type,true,Alloc,Tag>::ConstIterator
   DynamicMatrix<Type,true,Alloc,Tag>::cend( size_t j ) const noexcept
{
   BLAZE_USER_ASSERT( j < n_, "Invalid dense matrix column access index" );
   return ConstIterator( v_ + j*mm_ + m_ );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( const Type& rhs ) &
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
//
// This assignment operator offers the option to directly assign to all elements of the matrix
// by means of an initializer list:

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<int,columnMajor> A;
   A = { { 1, 2, 3 },
         { 4, 5 },
         { 7, 8, 9 } };
   \endcode

// The matrix is resized according to the given initializer list and all its elements are
// (copy) assigned the values from the given initializer list. Missing values are initialized
// as default (as e.g. the value 6 in the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( initializer_list< initializer_list<Type> > list ) &
{
   using blaze::clear;

   resize( list.size(), determineColumns( list ), false );

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
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::columnMajor;

   const int init[3][3] = { { 1, 2, 3 },
                            { 4, 5 },
                            { 7, 8, 9 } };
   blaze::DynamicMatrix<int,columnMajor> A;
   A = init;
   \endcode

// The matrix is resized according to the size of the static array and assigned the values of the
// given static array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the static array
        , size_t Rows     // Number of rows of the static array
        , size_t Cols >   // Number of columns of the static array
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( const Other (&array)[Rows][Cols] ) &
{
   resize( Rows, Cols, false );

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
//
// This assignment operator offers the option to directly set all elements of the matrix:

   \code
   using blaze::columnMajor;

   const std::array<std::array<int,3UL>,3UL> init{ { { 1, 2, 3 },
                                                     { 4, 5 },
                                                     { 7, 8, 9 } };
   blaze::DynamicMatrix<int,columnMajor> A;
   A = init;
   \endcode

// The matrix is resized according to the size of the std::array and assigned the values of the
// given std::array. Missing values are initialized with default values (as e.g. the value 6 in
// the example).
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename Other  // Data type of the std::array
        , size_t Rows     // Number of rows of the std::array
        , size_t Cols >   // Number of columns of the std::array
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( const std::array<std::array<Other,Cols>,Rows>& array ) &
{
   resize( Rows, Cols, false );

   for( size_t j=0UL; j<Cols; ++j )
      for( size_t i=0UL; i<Rows; ++i )
         v_[i+j*mm_] = array[i][j];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for DynamicMatrix.
//
// \param rhs Matrix to be copied.
// \return Reference to the assigned matrix.
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( const DynamicMatrix& rhs ) &
{
   if( &rhs == this ) return *this;

   resize( rhs.m_, rhs.n_, false );
   smpAssign( *this, *rhs );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Move assignment operator for DynamicMatrix.
//
// \param rhs The matrix to be moved into this instance.
// \return Reference to the assigned matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( DynamicMatrix&& rhs ) &
{
   blaze::destroy_n( v_, capacity_ );
   alloc_.deallocate( v_, capacity_ );

   m_        = rhs.m_;
   mm_       = rhs.mm_;
   n_        = rhs.n_;
   capacity_ = rhs.capacity_;
   v_        = rhs.v_;

   rhs.m_        = 0UL;
   rhs.mm_       = 0UL;
   rhs.n_        = 0UL;
   rhs.capacity_ = 0UL;
   rhs.v_        = nullptr;

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
//
// The matrix is resized according to the given \f$ M \times N \f$ matrix and initialized as a
// copy of this matrix.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT      // Type of the right-hand side matrix
        , bool SO >        // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator=( const Matrix<MT,SO>& rhs ) &
{
   using TT = decltype( trans( *this ) );
   using CT = decltype( ctrans( *this ) );
   using IT = decltype( inv( *this ) );

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<MT> );

   if( IsSame_v<MT,TT> && (*rhs).isAliased( this ) ) {
      transpose();
   }
   else if( IsSame_v<MT,CT> && (*rhs).isAliased( this ) ) {
      ctranspose();
   }
   else if( !IsSame_v<MT,IT> && (*rhs).canAlias( this ) ) {
      DynamicMatrix tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).rows(), (*rhs).columns(), false );
      if( IsSparseMatrix_v<MT> )
         reset();
      smpAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO >       // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator+=( const Matrix<MT,SO>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO >       // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator-=( const Matrix<MT,SO>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT     // Type of the right-hand side matrix
        , bool SO >       // Storage order of the right-hand side matrix
inline DynamicMatrix<Type,true,Alloc,Tag>&
   DynamicMatrix<Type,true,Alloc,Tag>::operator%=( const Matrix<MT,SO>& rhs ) &
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

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::rows() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::columns() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::spacing() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::capacity() const noexcept
{
   return capacity_;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::capacity( size_t j ) const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::nonZeros() const
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::nonZeros( size_t j ) const
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::reset()
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::reset( size_t j )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::clear()
{
   resize( 0UL, 0UL, false );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Changing the size of the matrix.
//
// \param m The new number of rows of the matrix.
// \param n The new number of columns of the matrix.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function resizes the matrix using the given size to \f$ m \times n \f$. During this
// operation, new dynamic memory may be allocated in case the capacity of the matrix is too
// small. Note that this function may invalidate all existing views (submatrices, rows, columns,
// ...) on the matrix if it is used to shrink the matrix. Additionally, the resize operation
// potentially changes all matrix elements. In order to preserve the old matrix values, the
// \a preserve flag can be set to \a true. However, new matrix elements of built-in type are
// not initialized!
//
// The following example illustrates the resize operation of a \f$ 2 \times 4 \f$ matrix to a
// \f$ 4 \times 2 \f$ matrix. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & 3 & 4 \\
                              5 & 6 & 7 & 8 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              5 & 6 \\
                              x & x \\
                              x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
void DynamicMatrix<Type,true,Alloc,Tag>::resize( size_t m, size_t n, bool preserve )
{
   using blaze::clear;
   using blaze::min;

   if( m == m_ && n == n_ ) return;

   const size_t mm( addPadding( m ) );

   if( preserve )
   {
      const size_t min_m( min( m, m_ ) );
      const size_t min_n( min( n, n_ ) );

      DynamicMatrix tmp( m, mm, n, Alloc{}, Uninitialized{} );

      for( size_t j=0UL; j<min_n; ++j ) {
         blaze::uninitialized_transfer( v_+j*mm_, v_+j*mm_+min_m, tmp.v_+j*mm );
         blaze::uninitialized_default_construct( tmp.v_+j*mm+min_m, tmp.v_+j*mm+mm );
      }
      blaze::uninitialized_default_construct( tmp.v_+min_n*mm, tmp.v_+mm*n );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }
   else if( mm*n > capacity_ )
   {
      DynamicMatrix tmp( m, mm, n, Alloc{}, Uninitialized{} );

      blaze::uninitialized_default_construct( tmp.v_, tmp.v_+tmp.capacity_ );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }

   if( IsVectorizable_v<Type> ) {
      for( size_t j=0UL; j<n; ++j )
         for( size_t i=m; i<mm; ++i )
            clear( v_[i+j*mm] );
   }

   m_  = m;
   mm_ = mm;
   n_  = n;

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Extending the size of the matrix.
//
// \param m Number of additional rows.
// \param n Number of additional columns.
// \param preserve \a true if the old values of the matrix should be preserved, \a false if not.
// \return void
//
// This function increases the matrix size by \a m rows and \a n columns. During this operation,
// new dynamic memory may be allocated in case the capacity of the matrix is too small. Therefore
// this function potentially changes all matrix elements. In order to preserve the old matrix
// values, the \a preserve flag can be set to \a true. However, new matrix elements of built-in
// type are not initialized!
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::extend( size_t m, size_t n, bool preserve )
{
   resize( m_+m, n_+n, preserve );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the matrix.
//
// \param elements The new minimum capacity of the dense matrix.
// \return void
//
// This function increases the capacity of the dense matrix to at least \a elements elements.
// The current values of the matrix elements are preserved.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::reserve( size_t elements )
{
   using blaze::clear;

   if( elements > capacity_ )
   {
      DynamicMatrix tmp( m_, mm_, n_, elements, Alloc{}, Uninitialized{} );

      blaze::uninitialized_transfer( v_, v_+capacity_, tmp.v_ );
      blaze::uninitialized_value_construct( tmp.v_+capacity_, tmp.v_+elements );

      std::swap( capacity_, tmp.capacity_ );
      std::swap( v_, tmp.v_ );
   }
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
// that due to padding the capacity might not be reduced exactly to rows() times columns().
// Please also note that in case a reallocation occurs, all iterators (including end() iterators),
// all pointers and references to elements of this matrix are invalidated.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::shrinkToFit()
{
   if( ( mm_ * n_ ) < capacity_ ) {
      DynamicMatrix( *this ).swap( *this );
   }
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void DynamicMatrix<Type,true,Alloc,Tag>::swap( DynamicMatrix& m ) noexcept
{
   using std::swap;

   swap( m_ , m.m_  );
   swap( mm_, m.mm_ );
   swap( n_ , m.n_  );
   swap( capacity_, m.capacity_ );
   swap( v_ , m.v_  );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Add the necessary amount of padding to the given value.
//
// \param value The value to be padded.
// \return The padded value.
//
// This function increments the given \a value by the necessary amount of padding based on the
// vector's data type \a Type.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline size_t DynamicMatrix<Type,true,Alloc,Tag>::addPadding( size_t values ) const noexcept
{
   if( usePadding && IsVectorizable_v<Type> )
      return nextMultiple<size_t>( values, SIMDSIZE );
   else return values;
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>& DynamicMatrix<Type,true,Alloc,Tag>::transpose()
{
   using std::swap;

   constexpr size_t block( BLOCK_SIZE );

   if( m_ == n_ )
   {
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( min( jj+block, n_ ) );
         for( size_t ii=0UL; ii<=jj; ii+=block ) {
            for( size_t j=jj; j<jend; ++j ) {
               const size_t iend( min( ii+block, m_, j ) );
               for( size_t i=ii; i<iend; ++i ) {
                  swap( v_[i+j*mm_], v_[j+i*mm_] );
               }
            }
         }
      }
   }
   else
   {
      DynamicMatrix tmp( trans(*this) );
      this->swap( tmp );
   }

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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline DynamicMatrix<Type,true,Alloc,Tag>& DynamicMatrix<Type,true,Alloc,Tag>::ctranspose()
{
   constexpr size_t block( BLOCK_SIZE );

   if( m_ == n_ )
   {
      for( size_t jj=0UL; jj<n_; jj+=block ) {
         const size_t jend( min( jj+block, n_ ) );
         for( size_t ii=0UL; ii<jj; ii+=block ) {
            const size_t iend( min( ii+block, m_ ) );
            for( size_t j=jj; j<jend; ++j ) {
               for( size_t i=ii; i<iend; ++i ) {
                  cswap( v_[i+j*mm_], v_[j+i*mm_] );
               }
            }
         }
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=jj; i<j; ++i ) {
               cswap( v_[i+j*mm_], v_[j+i*mm_] );
            }
            conjugate( v_[j+j*mm_] );
         }
      }
   }
   else
   {
      DynamicMatrix tmp( ctrans(*this) );
      this->swap( tmp );
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
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   A *= 4;        // Scaling of the matrix
   A.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the matrix
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline DynamicMatrix<Type,true,Alloc,Tag>& DynamicMatrix<Type,true,Alloc,Tag>::scale( const Other& scalar )
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
//  DEBUGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the dynamic matrix are intact.
//
// \return \a true in case the dynamic matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dynamic matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false.
*/
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,true,Alloc,Tag>::isIntact() const noexcept
{
   if( m_ > mm_ || mm_ * n_ > capacity_ )
      return false;

   if( IsVectorizable_v<Type> ) {
      for( size_t j=0UL; j<n_; ++j ) {
         for( size_t i=m_; i<mm_; ++i ) {
            if( !isDefault<strict>( v_[i+j*mm_] ) )
               return false;
         }
      }
   }

   return true;
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
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,true,Alloc,Tag>::canAlias( const Other* alias ) const noexcept
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
        , typename Alloc    // Type of the allocator
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool DynamicMatrix<Type,true,Alloc,Tag>::isAliased( const Other* alias ) const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,true,Alloc,Tag>::isAligned() const noexcept
{
   return ( usePadding || rows() % SIMDSIZE == 0UL );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool DynamicMatrix<Type,true,Alloc,Tag>::canSMPAssign() const noexcept
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,true,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,true,Alloc,Tag>::load( size_t i, size_t j ) const noexcept
{
   if( usePadding )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,true,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,true,Alloc,Tag>::loada( size_t i, size_t j ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= mm_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || i % SIMDSIZE == 0UL, "Invalid row access index" );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename DynamicMatrix<Type,true,Alloc,Tag>::SIMDType
   DynamicMatrix<Type,true,Alloc,Tag>::loadu( size_t i, size_t j ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= mm_, "Invalid row access index" );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,true,Alloc,Tag>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( usePadding )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,true,Alloc,Tag>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= mm_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || i % SIMDSIZE == 0UL, "Invalid row access index" );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,true,Alloc,Tag>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= mm_, "Invalid row access index" );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE void
   DynamicMatrix<Type,true,Alloc,Tag>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( i < m_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= mm_, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( !usePadding || i % SIMDSIZE == 0UL, "Invalid row access index" );
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::assign( const DenseMatrix<MT,true>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::assign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

   const size_t ipos( remainder ? prevMultiple( m_, SIMDSIZE ) : m_ );
   BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

   if( usePadding && useStreaming &&
       ( m_*n_ > ( cacheSize / ( sizeof(Type) * 3UL ) ) ) &&
       !(*rhs).isAliased( this ) )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::assign( const DenseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::assign( const SparseMatrix<MT,true>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::assign( const SparseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::addAssign( const DenseMatrix<MT,true>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::addAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::addAssign( const DenseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::addAssign( const SparseMatrix<MT,true>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::addAssign( const SparseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::subAssign( const DenseMatrix<MT,true>& rhs )
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
/*!\brief SIMD optimized implementation of the subtraction assignment of a column-major dense matrix.
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::subAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_DIAGONAL_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::subAssign( const DenseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::subAssign( const SparseMatrix<MT,true>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::subAssign( const SparseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::schurAssign( const DenseMatrix<MT,true>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( m_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= m_, "Invalid end calculation" );

   for( size_t j=0UL; j<n_; ++j ) {
      for( size_t i=0UL; (i+2UL) <= ipos; i+=2UL ) {
         v_[i  +j*mm_] *= (*rhs)(i  ,j);
         v_[i+1+j*mm_] *= (*rhs)(i+1,j);
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
/*!\brief SIMD optimized implementation of the Schur product assignment of a column-major dense matrix.
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline auto DynamicMatrix<Type,true,Alloc,Tag>::schurAssign( const DenseMatrix<MT,true>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   constexpr bool remainder( !usePadding || !IsPadded_v<MT> );

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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side dense matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::schurAssign( const DenseMatrix<MT,false>& rhs )
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::schurAssign( const SparseMatrix<MT,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<n_; ++j )
   {
      size_t i( 0UL );

      for( auto element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( v_[i+j*mm_] );
         v_[i+j*mm_] *= element->value();
         ++i;
      }

      for( ; i<m_; ++i ) {
         reset( v_[i+j*mm_] );
      }
   }
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
template< typename Type   // Data type of the matrix
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
template< typename MT >   // Type of the right-hand side sparse matrix
inline void DynamicMatrix<Type,true,Alloc,Tag>::schurAssign( const SparseMatrix<MT,false>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );

   BLAZE_INTERNAL_ASSERT( m_ == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( n_ == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<m_; ++i )
   {
      size_t j( 0UL );

      for( auto element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( v_[i+j*mm_] );
         v_[i+j*mm_] *= element->value();
         ++j;
      }

      for( ; j<n_; ++j ) {
         reset( v_[i+j*mm_] );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  DYNAMICMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DynamicMatrix operators */
//@{
template< typename Type, bool SO, typename Alloc, typename Tag >
void reset( DynamicMatrix<Type,SO,Alloc,Tag>& m );

template< typename Type, bool SO, typename Alloc, typename Tag >
void reset( DynamicMatrix<Type,SO,Alloc,Tag>& m, size_t i );

template< typename Type, bool SO, typename Alloc, typename Tag >
void clear( DynamicMatrix<Type,SO,Alloc,Tag>& m );

template< RelaxationFlag RF, typename Type, bool SO, typename Alloc, typename Tag >
bool isDefault( const DynamicMatrix<Type,SO,Alloc,Tag>& m );

template< typename Type, bool SO, typename Alloc, typename Tag >
bool isIntact( const DynamicMatrix<Type,SO,Alloc,Tag>& m ) noexcept;

template< typename Type, bool SO, typename Alloc, typename Tag >
void swap( DynamicMatrix<Type,SO,Alloc,Tag>& a, DynamicMatrix<Type,SO,Alloc,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given dynamic matrix.
// \ingroup dynamic_matrix
//
// \param m The matrix to be resetted.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void reset( DynamicMatrix<Type,SO,Alloc,Tag>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset the specified row/column of the given dynamic matrix.
// \ingroup dynamic_matrix
//
// \param m The matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given dynamic matrix to
// their default value. In case the given matrix is a \a rowMajor matrix the function resets the
// values in row \a i, if it is a \a columnMajor matrix the function resets the values in column
// \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void reset( DynamicMatrix<Type,SO,Alloc,Tag>& m, size_t i )
{
   m.reset( i );
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
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void clear( DynamicMatrix<Type,SO,Alloc,Tag>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dynamic matrix is in default state.
// \ingroup dynamic_matrix
//
// \param m The matrix to be tested for its default state.
// \return \a true in case the given matrix's rows and columns are zero, \a false otherwise.
//
// This function checks whether the dynamic matrix is in default (constructed) state, i.e. if
// it's number of rows and columns is 0. In case it is in default state, the function returns
// \a true, else it will return \a false. The following example demonstrates the use of the
// \a isDefault() function:

   \code
   blaze::DynamicMatrix<int> A;
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
        , typename Alloc     // Type of the allocator
        , typename Tag >     // Type tag
inline bool isDefault( const DynamicMatrix<Type,SO,Alloc,Tag>& m )
{
   return ( m.rows() == 0UL && m.columns() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given dynamic matrix are intact.
// \ingroup dynamic_matrix
//
// \param m The dynamic matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the dynamic matrix are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::DynamicMatrix<int> A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline bool isIntact( const DynamicMatrix<Type,SO,Alloc,Tag>& m ) noexcept
{
   return m.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two dynamic matrices.
// \ingroup dynamic_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Alloc  // Type of the allocator
        , typename Tag >  // Type tag
inline void swap( DynamicMatrix<Type,SO,Alloc,Tag>& a, DynamicMatrix<Type,SO,Alloc,Tag>& b ) noexcept
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct HasConstDataAccess< DynamicMatrix<T,SO,Alloc,Tag> >
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct HasMutableDataAccess< DynamicMatrix<T,SO,Alloc,Tag> >
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct IsAligned< DynamicMatrix<T,SO,Alloc,Tag> >
   : public BoolConstant<usePadding>
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct IsContiguous< DynamicMatrix<T,SO,Alloc,Tag> >
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct IsPadded< DynamicMatrix<T,SO,Alloc,Tag> >
   : public BoolConstant<usePadding>
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct IsResizable< DynamicMatrix<T,SO,Alloc,Tag> >
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
template< typename T, bool SO, typename Alloc, typename Tag >
struct IsShrinkable< DynamicMatrix<T,SO,Alloc,Tag> >
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
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( IsDenseMatrix_v<T1> || IsDenseMatrix_v<T2> ) &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) > >
{
   using ET = AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >;

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

   using Type = DynamicMatrix< ET
                             , SO
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
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
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( IsDenseMatrix_v<T1> || IsDenseMatrix_v<T2> ) &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) > >
{
   using ET = SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >;

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

   using Type = DynamicMatrix< ET
                             , SO
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
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
struct SchurTraitEval2< T1, T2
                      , EnableIf_t< IsDenseMatrix_v<T1> &&
                                    IsDenseMatrix_v<T2> &&
                                    ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                    ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                    ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                    ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                    ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                    ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                    ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) &&
                                    ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) > >
{
   using ET = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;

   static constexpr bool SO1 = StorageOrder_v<T1>;
   static constexpr bool SO2 = StorageOrder_v<T2>;

   static constexpr bool SO = ( IsSymmetric_v<T1> ^ IsSymmetric_v<T2>
                                ? ( IsSymmetric_v<T1>
                                    ? SO2
                                    : SO1 )
                                : SO1 && SO2 );

   using Type = DynamicMatrix< ET
                             , SO
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
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
                     , EnableIf_t< IsDenseMatrix_v<T1> &&
                                   IsScalar_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                   ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                   ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) > >
{
   using ET = MultTrait_t< ElementType_t<T1>, T2 >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<T1>
                             , DynamicAllocator_t< ET, GetAllocator_t<T1> >
                             , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsDenseMatrix_v<T2> &&
                                   ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                   ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                   ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) > >
{
   using ET = MultTrait_t< T1, ElementType_t<T2> >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<T2>
                             , DynamicAllocator_t< ET, GetAllocator_t<T2> >
                             , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsDenseVector_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   IsColumnVector_v<T1> &&
                                   IsRowVector_v<T2> &&
                                   ( ( Size_v<T1,0UL> == DefaultSize_v ) ||
                                     ( Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                   ( ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) ||
                                     ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ) > >
{
   using ET = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;

   using Type = DynamicMatrix< ET
                             , false
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
                             , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsDenseMatrix_v<T1> || IsDenseMatrix_v<T2> ) &&
                                   ( ( Size_v<T1,0UL> == DefaultSize_v &&
                                       ( !IsSquare_v<T1> || Size_v<T2,0UL> == DefaultSize_v ) ) ||
                                     ( Size_v<T2,1UL> == DefaultSize_v &&
                                       ( !IsSquare_v<T2> || Size_v<T1,1UL> == DefaultSize_v ) ) ) &&
                                   ( ( MaxSize_v<T1,0UL> == DefaultMaxSize_v &&
                                       ( !IsSquare_v<T1> || MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ) ||
                                     ( MaxSize_v<T2,1UL> == DefaultMaxSize_v &&
                                       ( !IsSquare_v<T2> || MaxSize_v<T1,1UL> == DefaultMaxSize_v ) ) ) > >
{
   using M1 = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using M2 = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;
   using ET = AddTrait_t<M1,M1>;

   using Type = DynamicMatrix< ET
                             , ( IsSparseMatrix_v<T1> ? StorageOrder_v<T2> : StorageOrder_v<T1> )
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
                             , AddTrait_t<M2,M2> >;
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
struct KronTraitEval2< T1, T2
                     , EnableIf_t< IsDenseMatrix_v<T1> &&
                                   IsDenseMatrix_v<T2> &&
                                   ( ( Size_v<T1,0UL> == DefaultSize_v ) ||
                                     ( Size_v<T2,0UL> == DefaultSize_v ) ||
                                     ( Size_v<T1,1UL> == DefaultSize_v ) ||
                                     ( Size_v<T2,1UL> == DefaultSize_v ) ) &&
                                   ( ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) ||
                                     ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ||
                                     ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) ||
                                     ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) ) > >
{
   using ET = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<T2>
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
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
                    , EnableIf_t< IsDenseMatrix_v<T1> &&
                                  IsScalar_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                  ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) > >
{
   using ET = DivTrait_t< ElementType_t<T1>, T2 >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<T1>
                             , DynamicAllocator_t< ET, GetAllocator_t<T1> >
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
                         , EnableIf_t< IsDenseMatrix_v<T> &&
                                       ( Size_v<T,0UL> == DefaultSize_v ||
                                         Size_v<T,1UL> == DefaultSize_v ) &&
                                       ( MaxSize_v<T,0UL> == DefaultMaxSize_v ||
                                         MaxSize_v<T,1UL> == DefaultMaxSize_v ) > >
{
   using ET =
      EvaluateTrait_t< decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) ) >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<T>
                             , DynamicAllocator_t< ET, GetAllocator_t<T> >
                             , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval2< T1, T2, OP
                          , EnableIf_t< IsColumnVector_v<T1> &&
                                        IsRowVector_v<T2> &&
                                        ( ( Size_v<T1,0UL> == DefaultSize_v ) ||
                                          ( Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                        ( ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) ||
                                          ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) ) > >
{
   using ET =
      EvaluateTrait_t< decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) ) >;

   using Type = DynamicMatrix< ET
                             , false
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
                             , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
};

template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval2< T1, T2, OP
                          , EnableIf_t< IsMatrix_v<T1> &&
                                        IsMatrix_v<T2> &&
                                        ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                        ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                        ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                        ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                        ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                        ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                        ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) &&
                                        ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) > >
{
   using ET =
      EvaluateTrait_t< decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) ) >;

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

   using Type = DynamicMatrix< ET
                             , SO
                             , DynamicAllocator_t< ET, GetAllocator_t<T1>, GetAllocator_t<T2> >
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
struct ExpandTraitEval2< T, E
                       , EnableIf_t< IsDenseVector_v<T> &&
                                     ( ( E == inf ) ||
                                       ( ( Size_v<T,0UL> == DefaultSize_v ) &&
                                         ( MaxSize_v<T,0UL> == DefaultMaxSize_v ) ) ) > >
{
   using Type = DynamicMatrix< ElementType_t<T>
                             , ( IsColumnVector_v<T> ? columnMajor : rowMajor )
                             , DynamicAllocator_t< ElementType_t<T>, GetAllocator_t<T> >
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
                       , EnableIf_t< IsDenseMatrix_v<T> &&
                                     ( ( R0 == inf && R1 == inf ) ||
                                       ( ( Size_v<T,0UL> == DefaultSize_v ) &&
                                         ( MaxSize_v<T,0UL> == DefaultMaxSize_v ) ) ||
                                       ( ( Size_v<T,1UL> == DefaultSize_v ) &&
                                         ( MaxSize_v<T,1UL> == DefaultMaxSize_v ) ) ) > >
{
   using Type = DynamicMatrix< ElementType_t<T>
                             , StorageOrder_v<T>
                             , DynamicAllocator_t< ElementType_t<T>, GetAllocator_t<T> >
                             , TagType_t<T> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SOLVETRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SolveTraitEval2< T1, T2
                      , EnableIf_t< IsDenseMatrix_v<T1> &&
                                    IsDenseMatrix_v<T2> &&
                                    ( ( ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                        ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                        ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                        ( MaxSize_v<T1,0UL> == DefaultMaxSize_v ) &&
                                        ( MaxSize_v<T2,0UL> == DefaultMaxSize_v ) &&
                                        ( MaxSize_v<T1,1UL> == DefaultMaxSize_v ) ) ||
                                      ( ( Size_v<T2,1UL> == DefaultSize_v ) &&
                                        ( MaxSize_v<T2,1UL> == DefaultMaxSize_v ) ) ) > >
{
   using Type = DynamicMatrix< ElementType_t<T2>
                             , StorageOrder_v<T2>
                             , DynamicAllocator_t< ElementType_t<T2>, GetAllocator_t<T2> >
                             , TagType_t<T2> >;
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
template< typename T1, bool SO, typename Alloc, typename Tag, typename T2 >
struct HighType< DynamicMatrix<T1,SO,Alloc,Tag>, DynamicMatrix<T2,SO,Alloc,Tag> >
{
   using Type = DynamicMatrix< typename HighType<T1,T2>::Type, SO, Alloc, Tag >;
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
template< typename T1, bool SO, typename Alloc, typename Tag, typename T2 >
struct LowType< DynamicMatrix<T1,SO,Alloc,Tag>, DynamicMatrix<T2,SO,Alloc,Tag> >
{
   using Type = DynamicMatrix< typename LowType<T1,T2>::Type, SO, Alloc, Tag >;
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
template< typename MT >
struct SubmatrixTraitEval2< MT, inf, inf, inf, inf
                          , EnableIf_t< IsDenseMatrix_v<MT> &&
                                        ( Size_v<MT,0UL> == DefaultSize_v ||
                                          Size_v<MT,1UL> == DefaultSize_v ) &&
                                        ( MaxSize_v<MT,0UL> == DefaultMaxSize_v ||
                                          MaxSize_v<MT,1UL> == DefaultMaxSize_v ) > >
{
   using ET = RemoveConst_t< ElementType_t<MT> >;

   using Type = DynamicMatrix< ET
                             , StorageOrder_v<MT>
                             , DynamicAllocator_t< ET, GetAllocator_t<MT> >
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
                     , EnableIf_t< IsDenseMatrix_v<MT> &&
                                   ( M == 0UL || Size_v<MT,1UL> == DefaultSize_v ) &&
                                   ( M == 0UL || MaxSize_v<MT,1UL> == DefaultMaxSize_v ) > >
{
   using ET = RemoveConst_t< ElementType_t<MT> >;

   using Type = DynamicMatrix< ET
                             , false
                             , DynamicAllocator_t< ET, GetAllocator_t<MT> >
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
                        , EnableIf_t< IsDenseMatrix_v<MT> &&
                                      ( N == 0UL || Size_v<MT,0UL> == DefaultSize_v ) &&
                                      ( N == 0UL || MaxSize_v<MT,0UL> == DefaultMaxSize_v ) > >
{
   using ET = RemoveConst_t< ElementType_t<MT> >;

   using Type = DynamicMatrix< ET
                             , true
                             , DynamicAllocator_t< ET, GetAllocator_t<MT> >
                             , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
