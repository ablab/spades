//=================================================================================================
/*!
//  \file blaze/math/dense/CustomVector.h
//  \brief Header file for the implementation of a customizable vector
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

#ifndef _BLAZE_MATH_DENSE_CUSTOMVECTOR_H_
#define _BLAZE_MATH_DENSE_CUSTOMVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <array>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/DenseIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/PaddingFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/typetraits/CustomTransposeType.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsCustom.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
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
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsVectorizable.h>
#include <blaze/util/typetraits/RemoveConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup custom_vector CustomVector
// \ingroup dense_vector
*/
/*!\brief Efficient implementation of a customizable vector.
// \ingroup custom_vector
//
// \section customvector_general General
//
// The CustomVector class template provides the functionality to represent an external array of
// elements of arbitrary type and a fixed size as a native \b Blaze dense vector data structure.
// Thus in contrast to all other dense vector types a custom vector does not perform any kind
// of memory allocation by itself, but it is provided with an existing array of element during
// construction. A custom vector can therefore be considered an alias to the existing array.
//
// The type of the elements, the properties of the given array of elements, the transpose flag,
// and the group tag of the vector can be specified via the following five template parameters:

   \code
   namespace blaze {

   template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag >
   class CustomVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. CustomVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - AF  : specifies whether the represented, external arrays are properly aligned with
//          respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::aligned
//          or \c blaze::unaligned).
//  - PF  : specified whether the represented, external arrays are properly padded with
//          respect to the available instruction set (SSE, AVX, ...) or not (\c blaze::padded
//          or \c blaze::unpadded).
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// The following examples give an impression of several possible types of custom vectors:

   \code
   using blaze::CustomVector;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a custom column vector for unaligned, unpadded integer arrays
   using UnalignedUnpadded = CustomVector<int,unaligned,unpadded,columnVector>;

   // Definition of a custom column vector for unaligned but padded 'float' arrays
   using UnalignedPadded = CustomVector<float,unaligned,padded,columnVector>;

   // Definition of a custom row vector for aligned, unpadded 'double' arrays
   using AlignedUnpadded = CustomVector<double,aligned,unpadded,rowVector>;

   // Definition of a custom row vector for aligned, padded 'complex<double>' arrays
   using AlignedPadded = CustomVector<complex<double>,aligned,padded,rowVector>;
   \endcode

// \n \section customvector_special_properties Special Properties of Custom Vectors
//
// In comparison with the remaining \b Blaze dense vector types CustomVector has several special
// characteristics. All of these result from the fact that a custom vector is not performing any
// kind of memory allocation, but instead is given an existing array of elements. The following
// sections discuss all of these characteristics:
//
//  -# <b>\ref customvector_memory_management</b>
//  -# <b>\ref customvector_copy_operations</b>
//  -# <b>\ref customvector_alignment</b>
//  -# <b>\ref customvector_padding</b>
//
// \n \subsection customvector_memory_management Memory Management
//
// The CustomVector class template acts as an adaptor for an existing array of elements. As such
// it provides everything that is required to use the array just like a native \b Blaze dense
// vector data structure. However, this flexibility comes with the price that the user of a custom
// vector is responsible for the resource management.
//
// The following examples give an impression of several possible types of custom vectors:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Definition of a 3-dimensional custom vector with unaligned, unpadded and externally
   // managed integer array. Note that the std::vector must be guaranteed to outlive the
   // custom vector!
   std::vector<int> vec( 3UL );
   CustomVector<int,unaligned,unpadded> a( &vec[0], 3UL );

   // Definition of a custom vector with size 3 and capacity 16 with aligned, padded and
   // externally managed integer array. Note that the std::unique_ptr must be guaranteed
   // to outlive the custom vector!
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 16UL ) );
   CustomVector<int,aligned,padded> b( memory.get(), 3UL, 16UL );
   \endcode

// \n \subsection customvector_copy_operations Copy Operations
//
// As with all dense vectors it is possible to copy construct a custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaligned;
   using blaze::unpadded;

   using CustomType = CustomVector<int,unaligned,unpadded>;

   std::vector<int> vec( 5UL, 10 );  // Vector of 5 integers of the value 10
   CustomType a( &vec[0], 5UL );     // Represent the std::vector as Blaze dense vector
   a[1] = 20;                        // Also modifies the std::vector

   CustomType b( a );  // Creating a copy of vector a
   b[2] = 20;          // Also affects vector a and the std::vector
   \endcode

// It is important to note that a custom vector acts as a reference to the specified array. Thus
// the result of the copy constructor is a new custom vector that is referencing and representing
// the same array as the original custom vector.
//
// In contrast to copy construction, just as with references, copy assignment does not change
// which array is referenced by the custom vector, but modifies the values of the array:

   \code
   std::vector<int> vec2( 5UL, 4 );  // Vector of 5 integers of the value 4
   CustomType c( &vec2[0], 5UL );    // Represent the std::vector as Blaze dense vector

   a = c;  // Copy assignment: Set all values of vector a and b to 4.
   \endcode

// \n \subsection customvector_alignment Alignment
//
// In case the custom vector is specified as \c aligned the passed array must be guaranteed to
// be aligned according to the requirements of the used instruction set (SSE, AVX, ...). For
// instance, if AVX is active an array of integers must be 32-bit aligned:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   // Allocation of 32-bit aligned memory
   std::unique_ptr<int[],Deallocate> memory( allocate<int>( 5UL ) );

   CustomVector<int,aligned,unpadded> a( memory.get(), 5UL );
   \endcode

// In case the alignment requirements are violated, a \c std::invalid_argument exception is
// thrown.
//
// \n \subsection customvector_padding Padding
//
// Adding padding elements to the end of an array can have a significant impact on performance.
// For instance, assuming that AVX is available, then two aligned, padded, 3-dimensional vectors
// of double precision values can be added via a single SIMD addition operations:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::padded;

   using CustomType = CustomVector<double,aligned,padded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 4UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 4UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 4UL ) );

   // Creating padded custom vectors of size 3 and a capacity of 4
   CustomType a( memory1.get(), 3UL, 4UL );
   CustomType b( memory2.get(), 3UL, 4UL );
   CustomType c( memory3.get(), 3UL, 4UL );

   // ... Initialization

   c = a + b;  // AVX-based vector addition
   \endcode

// In this example, maximum performance is possible. However, in case no padding elements are
// inserted, a scalar addition has to be used:

   \code
   using blaze::CustomVector;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unpadded;

   using CustomType = CustomVector<double,aligned,unpadded>;

   std::unique_ptr<double[],Deallocate> memory1( allocate<double>( 3UL ) );
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 3UL ) );
   std::unique_ptr<double[],Deallocate> memory3( allocate<double>( 3UL ) );

   // Creating unpadded custom vector of size 3
   CustomType a( allocate<double>( 3UL ), 3UL );
   CustomType b( allocate<double>( 3UL ), 3UL );
   CustomType c( allocate<double>( 3UL ), 3UL );

   // ... Initialization

   c = a + b;  // Scalar vector addition
   \endcode

// Note the different number of constructor parameters for unpadded and padded custom vectors:
// In contrast to unpadded vectors, where during the construction only the size of the array
// has to be specified, during the construction of a padded custom vector it is additionally
// necessary to explicitly specify the capacity of the array.
//
// The number of padding elements is required to be sufficient with respect to the available
// instruction set: In case of an aligned padded custom vector the added padding elements must
// guarantee that the capacity is greater or equal than the size and a multiple of the SIMD vector
// width. In case of unaligned padded vectors the number of padding elements can be greater or
// equal the number of padding elements of an aligned padded custom vector. In case the padding
// is insufficient with respect to the available instruction set, a \c std::invalid_argument
// exception is thrown.
//
// Please also note that \b Blaze will zero initialize the padding elements in order to achieve
// maximum performance!
//
//
// \n \section customvector_arithmetic_operations Arithmetic Operations
//
// The use of custom vectors in arithmetic operations is designed to be as natural and intuitive
// as possible. All operations (addition, subtraction, multiplication, scaling, ...) can be
// expressed similar to a text book representation. Also, custom vectors can be combined with all
// other dense and sparse vectors and matrices. The following example gives an impression of the
// use of CustomVector:

   \code
   using blaze::CustomVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;
   using blaze::Deallocate;
   using blaze::allocate;
   using blaze::aligned;
   using blaze::unaligned;
   using blaze::padded;
   using blaze::unpadded;

   // Non-initialized custom column vector of size 2. All given arrays are considered to be
   // unaligned and unpadded. The memory is managed via a 'std::vector'.
   std::vector<double> memory1( 2UL );
   CustomVector<double,unaligned,unpadded> a( memory1.data(), 2UL );

   a[0] = 1.0;  // Initialization of the first element
   a[1] = 2.0;  // Initialization of the second element

   // Non-initialized custom column vector of size 2 and capacity 4. All given arrays are required
   // to be properly aligned and padded. The memory is managed via a 'std::unique_ptr'.
   std::unique_ptr<double[],Deallocate> memory2( allocate<double>( 4UL ) );
   CustomVector<double,aligned,padded> b( memory2.get(), 2UL, 4UL );

   b = 2.0;  // Homogeneous initialization of all elements

   CompressedVector<float> c( 2 );  // Empty sparse single precision vector
   DynamicVector<double>   d;       // Default constructed dynamic vector
   DynamicMatrix<double>   A;       // Default constructed row-major matrix

   d = a + b;  // Vector addition between custom vectors of equal element type
   d = a - c;  // Vector subtraction between a dense and sparse vector with different element types
   d = a * b;  // Component-wise vector multiplication

   a *= 2.0;      // In-place scaling of vector
   d  = a * 2.0;  // Scaling of vector a
   d  = 2.0 * a;  // Scaling of vector a

   d += a - b;  // Addition assignment
   d -= a + c;  // Subtraction assignment
   d *= a * b;  // Multiplication assignment

   double scalar = trans( a ) * b;  // Scalar/dot/inner product between two vectors

   A = a * trans( b );  // Outer product between two vectors
   \endcode
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
class CustomVector
   : public DenseVector< CustomVector<Type,AF,PF,TF,Tag,RT>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This     = CustomVector<Type,AF,PF,TF,Tag,RT>;  //!< Type of this CustomVector instance.
   using BaseType = DenseVector<This,TF>;                //!< Base type of this CustomVector instance.

   //! Result type for expression template evaluations.
   using ResultType = RT;

   //! Transpose type for expression template evaluations.
   using TransposeType = CustomTransposeType_t<RT>;

   using ElementType   = Type;                      //!< Type of the vector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using TagType       = Tag;                       //!< Tag type of this CustomVector instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations
   using CompositeType = const CustomVector&;       //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = Type*;        //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = DenseIterator<Type,AF>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,AF>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CustomVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using RRT   = Rebind_t< RT, RemoveConst_t<NewType> >;  //!< The rebound result type.
      using Other = CustomVector<NewType,AF,PF,TF,Tag,RRT>;  //!< The type of the other CustomVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CustomVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using RRT   = Resize_t<RT,NewN>;                    //!< The resized result type.
      using Other = CustomVector<Type,AF,PF,TF,Tag,RRT>;  //!< The type of the other CustomVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operations. In case the element type of the vector is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = IsVectorizable_v<Type>;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline CustomVector();
   inline CustomVector( Type* ptr, size_t n );
   inline CustomVector( Type* ptr, size_t n, size_t nn );
   inline CustomVector( const CustomVector& v );
   inline CustomVector( CustomVector&& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CustomVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin () noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline Iterator       end   () noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline CustomVector& operator=( const Type& rhs );
   inline CustomVector& operator=( initializer_list<Type> list );

   template< typename Other, size_t Dim >
   inline CustomVector& operator=( const Other (&array)[Dim] );

   template< typename Other, size_t Dim >
   inline CustomVector& operator=( const std::array<Other,Dim>& array );

   inline CustomVector& operator=( const CustomVector& rhs );
   inline CustomVector& operator=( CustomVector&& rhs ) noexcept;

   template< typename VT > inline CustomVector& operator= ( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator+=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator-=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator*=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator/=( const DenseVector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator%=( const Vector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   clear();
   inline void   swap( CustomVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline CustomVector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Resource management functions***************************************************************
   /*!\name Resource management functions */
   //@{
   inline void reset( Type* ptr, size_t n );
   inline void reset( Type* ptr, size_t n, size_t nn );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDAdd_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDSub_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDMult_v< Type, ElementType_t<VT> > );
   /*! \endcond */
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDDiv_v< Type, ElementType_t<VT> > );
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

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;  //!< The size/dimension of the custom vector.
   Type* v_;      //!< The custom array of elements.
                  /*!< Access to the array of elements is gained via the
                       subscript operator. The order of the elements is
                       \f[\left(\begin{array}{*{5}{c}}
                       0 & 1 & 2 & \cdots & N-1 \\
                       \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( Type );
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
/*!\brief The default constructor for CustomVector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::CustomVector()
   : size_( 0UL )      // The size/dimension of the vector
   , v_   ( nullptr )  // The custom array of elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an unpadded custom vector of size \a n.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This constructor creates an unpadded custom vector of size \a n. The construction fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This constructor is \b NOT available for padded custom vectors!
// \note The custom vector does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::CustomVector( Type* ptr, size_t n )
   : size_( n )    // The size/dimension of the vector
   , v_   ( ptr )  // The custom array of elements
{
   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && !checkAlignment( ptr ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a padded custom vector of size \a n and capacity \a nn.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \param nn The maximum size of the given array.
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This constructor creates a padded custom vector of size \a n and capacity \a nn. The
// construction fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...);
//  - ... the specified capacity \a nn is insufficient for the given data type \a Type and the
//    available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This constructor is \b NOT available for unpadded custom vectors!
// \note The custom vector does \b NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::CustomVector( Type* ptr, size_t n, size_t nn )
   : size_( 0UL )      // The size/dimension of the vector
   , v_   ( nullptr )  // The custom array of elements
{
   BLAZE_STATIC_ASSERT( PF == padded );

   MAYBE_UNUSED( ptr, n, nn );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for CustomVector.
//
// \param v Vector to be copied.
//
// The copy constructor initializes the custom vector as an exact copy of the given custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::CustomVector( const CustomVector& v )
   : size_( v.size_ )  // The size/dimension of the vector
   , v_   ( v.v_ )     // The custom array of elements
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for CustomVector.
//
// \param v The vector to be moved into this instance.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::CustomVector( CustomVector&& v ) noexcept
   : size_( v.size_ )             // The size/dimension of the vector
   , v_   ( v.v_ )                // The custom array of elements
{
   v.size_ = 0UL;
   v.v_    = nullptr;

   BLAZE_INTERNAL_ASSERT( v.data() == nullptr, "Invalid data reference detected" );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for CustomVector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>::~CustomVector()
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<RT> );
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
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::Reference
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstReference
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::Reference
   CustomVector<Type,AF,PF,TF,Tag,RT>::at( size_t index )
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstReference
   CustomVector<Type,AF,PF,TF,Tag,RT>::at( size_t index ) const
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::Pointer
   CustomVector<Type,AF,PF,TF,Tag,RT>::data() noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstPointer
   CustomVector<Type,AF,PF,TF,Tag,RT>::data() const noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::Iterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::begin() noexcept
{
   return Iterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::begin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::cbegin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::Iterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::end() noexcept
{
   return Iterator( v_+size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::end() const noexcept
{
   return ConstIterator( v_+size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,PF,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,PF,TF,Tag,RT>::cend() const noexcept
{
   return ConstIterator( v_+size_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Homogenous assignment to all vector elements.
//
// \param rhs Scalar value to be assigned to all vector elements.
// \return Reference to the assigned vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( const Type& rhs )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] = rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to custom vector.
//
// This assignment operator offers the option to directly assign to all elements of the vector
// by means of an initializer list:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   const int array[4] = { 1, 2, 3, 4 };

   CustomVector<double,unaligned,unpadded> v( array, 4UL );
   v = { 5, 6, 7 };
   \endcode

// The vector elements are assigned the values from the given initializer list. Missing values
// are reset to their default state. Note that in case the size of the initializer list exceeds
// the size of the vector, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( initializer_list<Type> list )
{
   if( list.size() > size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to custom vector" );
   }

   std::fill( std::copy( list.begin(), list.end(), v_ ), v_+size_, Type() );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all vector elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the vector. The
// following example demonstrates this by means of an unaligned, unpadded custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   const int array[4] = { 1, 2, 3, 4 };
   const int init[4]  = { 5, 6, 7 };

   CustomVector<double,unaligned,unpadded> v( array, 4UL );
   v = init;
   \endcode

// The vector is assigned the values from the given static array. Missing values are initialized
// with default values (as e.g. the fourth element in the example). Note that the size of the
// static array must match the size of the custom vector. Otherwise a \a std::invalid_argument
// exception is thrown. Also note that after the assignment \a array will have the same entries
// as \a init.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Dim >      // Dimension of the static array
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( const Other (&array)[Dim] )
{
   if( size_ != Dim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all vector elements.
//
// \param array The given std::array for the assignment.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid std::array size.
//
// This assignment operator offers the option to directly set all elements of the vector. The
// following example demonstrates this by means of an unaligned, unpadded custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   const int array[4] = { 1, 2, 3, 4 };
   const std::array<int,4UL> init{ 5, 6, 7 };

   CustomVector<double,unaligned,unpadded> v( array, 4UL );
   v = init;
   \endcode

// The vector is assigned the values from the given std::array. Missing values are initialized
// with default values (as e.g. the fourth element in the example). Note that the size of the
// std::array must match the size of the custom vector. Otherwise a \a std::invalid_argument
// exception is thrown. Also note that after the assignment \a array will have the same entries
// as \a init.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the std::array
        , size_t Dim >      // Dimension of the std::array
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( const std::array<Other,Dim>& array )
{
   if( size_ != Dim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for CustomVector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// The vector is initialized as a copy of the given vector. In case the current sizes of the two
// vectors don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( const CustomVector& rhs )
{
   if( rhs.size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   smpAssign( *this, *rhs );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for CustomVector.
//
// \param rhs The vector to be moved into this instance.
// \return Reference to the assigned vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( CustomVector&& rhs ) noexcept
{
   size_ = rhs.size_;
   v_    = rhs.v_;

   rhs.size_ = 0UL;
   rhs.v_    = nullptr;

   BLAZE_INTERNAL_ASSERT( rhs.data() == nullptr, "Invalid data reference detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// The vector is initialized as a copy of the given vector. In case the current sizes of the two
// vectors don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator+=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      smpAddAssign( *this, tmp );
   }
   else {
      smpAddAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator-=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      smpSubAssign( *this, tmp );
   }
   else {
      smpSubAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator*=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( MultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( IsSparseVector_v<VT> || (*rhs).canAlias( this ) ) {
      const MultType tmp( *this * (*rhs) );
      if( IsSparseVector_v<MultType> )
         reset();
      smpAssign( *this, tmp );
   }
   else {
      smpMultAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator/=( const DenseVector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( DivType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( DivType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const DivType tmp( *this / (*rhs) );
      smpAssign( *this, tmp );
   }
   else {
      smpDivAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the vector.
// \exception std::invalid_argument Invalid vector size for cross product.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::operator%=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size_ != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );
   assign( *this, tmp );

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,PF,TF,Tag,RT>::size() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the minimum capacity of the vector.
//
// \return The minimum capacity of the vector.
//
// This function returns the minimum capacity of the vector, which corresponds to the current
// size plus padding.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,PF,TF,Tag,RT>::spacing() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
//
// \return The maximum capacity of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,PF,TF,Tag,RT>::capacity() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the vector.
//
// \return The number of non-zero elements in the vector.
//
// This function returns the number of non-zero elements in the vector (i.e. the elements that
// compare unequal to their default value). Note that the number of non-zero elements is always
// less than or equal to the current size of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,PF,TF,Tag,RT>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( size_t i=0UL; i<size_; ++i ) {
      if( !isDefault<strict>( v_[i] ) )
         ++nonzeros;
   }

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::reset()
{
   using blaze::clear;
   for( size_t i=0UL; i<size_; ++i )
      clear( v_[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the vector to its default state.
//
// \return void
//
// This function clears the vector to its default state. In case the vector has been passed the
// responsibility to manage the given array, it disposes the resource via the specified deleter.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::clear()
{
   size_ = 0UL;
   v_ = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
//
// \param v The vector to be swapped.
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::swap( CustomVector& v ) noexcept
{
   using std::swap;

   swap( size_, v.size_ );
   swap( v_, v.v_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scaling of the vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the vector.
//
// This function scales the vector by applying the given scalar value \a scalar to each element
// of the vector. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::unpadded;

   CustomVector<double,unaligned,unpadded> v( ... );

   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the scalar value
inline CustomVector<Type,AF,PF,TF,Tag,RT>&
   CustomVector<Type,AF,PF,TF,Tag,RT>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  RESOURCE MANAGEMENT FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Resets the custom vector and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \return void
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This function resets the custom vector to the given array of elements of size \a n. The
// function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...).
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function is \b NOT available for padded custom vectors!
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom vector referencing the array goes out of scope.
// \note The custom vector does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::reset( Type* ptr, size_t n )
{
   CustomVector tmp( ptr, n );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resets the custom vector and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \param nn The maximum size of the given array.
// \return void
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This function resets the custom vector to the given array of elements of size \a n and
// capacity \a nn. The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...);
//  - ... the specified capacity \a nn is insufficient for the given data type \a Type and
//    the available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note This function is \a NOT available for unpadded custom vectors!
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom vector referencing the array goes out of scope.
// \note The custom vector does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::reset( Type* ptr, size_t n, size_t nn )
{
   BLAZE_STATIC_ASSERT( PF == padded );

   MAYBE_UNUSED( ptr, n, nn );
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the vector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomVector<Type,AF,PF,TF,Tag,RT>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address is aliased with the vector. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomVector<Type,AF,PF,TF,Tag,RT>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector is properly aligned in memory.
//
// \return \a true in case the vector is aligned, \a false if not.
//
// This function returns whether the vector is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of the vector are guaranteed to conform to the alignment
// restrictions of the element type \a Type.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomVector<Type,AF,PF,TF,Tag,RT>::isAligned() const noexcept
{
   return ( AF || checkAlignment( v_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector can be used in SMP assignments.
//
// \return \a true in case the vector can be used in SMP assignments, \a false if not.
//
// This function returns whether the vector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// vector).
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomVector<Type,AF,PF,TF,Tag,RT>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,PF,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,PF,TF,Tag,RT>::load( size_t index ) const noexcept
{
   if( AF )
      return loada( index );
   else
      return loadu( index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,PF,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,PF,TF,Tag,RT>::loada( size_t index ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( !AF || index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   return loada( v_+index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,PF,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,PF,TF,Tag,RT>::loadu( size_t index ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index< size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );

   return loadu( v_+index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,PF,TF,Tag,RT>::store( size_t index, const SIMDType& value ) noexcept
{
   if( AF )
      storea( index, value );
   else
      storeu( index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,PF,TF,Tag,RT>::storea( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( !AF || index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   storea( v_+index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense vector.
// The index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,PF,TF,Tag,RT>::storeu( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );

   storeu( v_+index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned, non-temporal store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense vector. The index must be smaller than the number of vector elements and it must be
// a multiple of the number of values inside the SIMD element. This function must \b NOT be
// called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,PF,TF,Tag,RT>::stream( size_t index, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( !AF || index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   stream( v_+index, value );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::assign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] = (*rhs)[i    ];
      v_[i+1UL] = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] = (*rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::assign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   if( AF && useStreaming && size_ > ( cacheSize/( sizeof(Type) * 3UL ) ) && !(*rhs).isAliased( this ) )
   {
      size_t i( 0UL );

      for( ; i<ipos; i+=SIMDSIZE ) {
         stream( i, (*rhs).load(i) );
      }
      for( ; i<size_; ++i ) {
         v_[i] = (*rhs)[i];
      }
   }
   else
   {
      const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
      BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

      size_t i( 0UL );
      ConstIterator_t<VT> it( (*rhs).begin() );

      for( ; i<i4way; i+=SIMDSIZE*4UL ) {
         store( i             , it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE    , it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE*2UL, it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE*3UL, it.load() ); it += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
         store( i, it.load() );
      }
      for( ; i<size_; ++i, ++it ) {
         v_[i] = *it;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] = element->value();
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::addAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] += (*rhs)[i    ];
      v_[i+1UL] += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] += (*rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::addAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) + it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) + it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      v_[i] += *it;
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] += element->value();
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::subAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] -= (*rhs)[i    ];
      v_[i+1UL] -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] -= (*rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::subAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) - it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) - it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      v_[i] -= *it;
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] -= element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::multAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] *= (*rhs)[i    ];
      v_[i+1UL] *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] *= (*rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::multAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) * it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) * it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      v_[i] *= *it;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,PF,TF,Tag,RT>::multAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] = tmp[element->index()] * element->value();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::divAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] /= (*rhs)[i    ];
      v_[i+1UL] /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] /= (*rhs)[ipos];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,PF,TF,Tag,RT>::divAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) / it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) / it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      v_[i] /= *it;
   }
}
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR PADDED VECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of CustomVector for padded vectors.
// \ingroup custom_vector
//
// This specialization of CustomVector adapts the class template to the requirements of padded
// vectors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
class CustomVector<Type,AF,padded,TF,Tag,RT>
   : public DenseVector< CustomVector<Type,AF,padded,TF,Tag,RT>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This     = CustomVector<Type,AF,padded,TF,Tag,RT>;  //!< Type of this CustomVector instance.
   using BaseType = DenseVector<This,TF>;                    //!< Base type of this CustomVector instance.

   //! Result type for expression template evaluations.
   using ResultType = RT;

   //! Transpose type for expression template evaluations.
   using TransposeType = CustomTransposeType_t<RT>;

   using ElementType   = Type;                      //!< Type of the vector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using TagType       = Tag;                       //!< Tag type of this CustomVector instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations
   using CompositeType = const CustomVector&;       //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = Type*;        //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = DenseIterator<Type,AF>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,AF>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CustomVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using RRT   = Rebind_t< RT, RemoveConst_t<NewType> >;      //!< The rebound result type.
      using Other = CustomVector<NewType,AF,padded,TF,Tag,RRT>;  //!< The type of the other CustomVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CustomVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using RRT   = Resize_t<RT,NewN>;                        //!< The resized result type.
      using Other = CustomVector<Type,AF,padded,TF,Tag,RRT>;  //!< The type of the other CustomVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operations. In case the element type of the vector is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = IsVectorizable_v<Type>;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline CustomVector();
   inline CustomVector( Type* ptr, size_t n, size_t nn );
   inline CustomVector( const CustomVector& v );
   inline CustomVector( CustomVector&& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CustomVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin () noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline Iterator       end   () noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline CustomVector& operator=( const Type& rhs );
   inline CustomVector& operator=( initializer_list<Type> list );

   template< typename Other, size_t Dim >
   inline CustomVector& operator=( const Other (&array)[Dim] );

   template< typename Other, size_t Dim >
   inline CustomVector& operator=( const std::array<Other,Dim>& array );

   inline CustomVector& operator=( const CustomVector& rhs );
   inline CustomVector& operator=( CustomVector&& rhs ) noexcept;

   template< typename VT > inline CustomVector& operator= ( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator+=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator-=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator*=( const Vector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator/=( const DenseVector<VT,TF>& rhs );
   template< typename VT > inline CustomVector& operator%=( const Vector<VT,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   clear();
   inline void   swap( CustomVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline CustomVector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Resource management functions***************************************************************
   /*!\name Resource management functions */
   //@{
   inline void reset( Type* ptr, size_t n, size_t nn );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT::simdEnabled &&
        IsSIMDCombinable_v< Type, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDAdd_v< Type, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDSub_v< Type, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDMult_v< Type, ElementType_t<VT> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT> &&
        HasSIMDDiv_v< Type, ElementType_t<VT> > );
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

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT >
   inline auto assign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT> >;

   template< typename VT > inline void assign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT >
   inline auto addAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT> >;

   template< typename VT > inline void addAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT >
   inline auto subAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT> >;

   template< typename VT > inline void subAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT >
   inline auto multAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT> >;

   template< typename VT > inline void multAssign( const SparseVector<VT,TF>& rhs );

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT> >;

   template< typename VT >
   inline auto divAssign( const DenseVector<VT,TF>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;      //!< The size/dimension of the custom vector.
   size_t capacity_;  //!< The maximum capacity of the custom vector.
   Type* v_;          //!< The custom array of elements.
                      /*!< Access to the array of elements is gained via the
                           subscript operator. The order of the elements is
                           \f[\left(\begin{array}{*{5}{c}}
                           0 & 1 & 2 & \cdots & N-1 \\
                           \end{array}\right)\f] */
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( Type );
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
/*!\brief The default constructor for CustomVector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>::CustomVector()
   : size_    ( 0UL )      // The size/dimension of the vector
   , capacity_( 0UL )      // The maximum capacity of the vector
   , v_       ( nullptr )  // The custom array of elements
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for a padded custom vector of size \a n and capacity \a nn.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \param nn The maximum size of the given array.
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This constructor creates a padded custom vector of size \a n and capacity \a nn. The
// construction of the vector fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the alignment flag \a AF is set to \a aligned, but the passed pointer is not properly
//    aligned according to the available instruction set (SSE, AVX, ...);
//  - ... the specified capacity \a nn is insufficient for the given data type \a Type and
//    the available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note The custom vector does NOT take responsibility for the given array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>::CustomVector( Type* ptr, size_t n, size_t nn )
   : size_    ( n )    // The size/dimension of the vector
   , capacity_( nn )   // The maximum capacity of the vector
   , v_       ( ptr )  // The custom array of elements
{
   using blaze::clear;

   if( ptr == nullptr ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array of elements" );
   }

   if( AF && !checkAlignment( ptr ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid alignment detected" );
   }

   if( IsVectorizable_v<Type> && capacity_ < nextMultiple<size_t>( size_, SIMDSIZE ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Insufficient capacity for padded vector" );
   }

   if( IsVectorizable_v<Type> ) {
      for( size_t i=size_; i<capacity_; ++i )
         clear( v_[i] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The copy constructor for CustomVector.
//
// \param v Vector to be copied.
//
// The copy constructor initializes the custom vector as an exact copy of the given custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>::CustomVector( const CustomVector& v )
   : size_    ( v.size_ )      // The size/dimension of the vector
   , capacity_( v.capacity_ )  // The maximum capacity of the vector
   , v_       ( v.v_ )         // The custom array of elements
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The move constructor for CustomVector.
//
// \param v The vector to be moved into this instance.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>::CustomVector( CustomVector&& v ) noexcept
   : size_    ( v.size_ )             // The size/dimension of the vector
   , capacity_( v.capacity_ )         // The maximum capacity of the vector
   , v_       ( v.v_ )                // The custom array of elements
{
   v.size_     = 0UL;
   v.capacity_ = 0UL;
   v.v_        = nullptr;

   BLAZE_INTERNAL_ASSERT( v.data() == nullptr, "Invalid data reference detected" );
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
/*!\brief The destructor for CustomVector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>::~CustomVector()
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( RT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( RT, TF );
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<RT> );
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
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::Reference
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstReference
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::Reference
   CustomVector<Type,AF,padded,TF,Tag,RT>::at( size_t index )
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstReference
   CustomVector<Type,AF,padded,TF,Tag,RT>::at( size_t index ) const
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::Pointer
   CustomVector<Type,AF,padded,TF,Tag,RT>::data() noexcept
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstPointer
   CustomVector<Type,AF,padded,TF,Tag,RT>::data() const noexcept
{
   return v_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::Iterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::begin() noexcept
{
   return Iterator( v_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::begin() const noexcept
{
   return ConstIterator( v_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the custom vector.
//
// \return Iterator to the first element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::cbegin() const noexcept
{
   return ConstIterator( v_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::Iterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::end() noexcept
{
   return Iterator( v_+size_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::end() const noexcept
{
   return ConstIterator( v_+size_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the custom vector.
//
// \return Iterator just past the last element of the custom vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline typename CustomVector<Type,AF,padded,TF,Tag,RT>::ConstIterator
   CustomVector<Type,AF,padded,TF,Tag,RT>::cend() const noexcept
{
   return ConstIterator( v_+size_ );
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
/*!\brief Homogenous assignment to all vector elements.
//
// \param rhs Scalar value to be assigned to all vector elements.
// \return Reference to the assigned vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( const Type& rhs )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] = rhs;
   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to custom vector.
//
// This assignment operator offers the option to directly assign to all elements of the vector
// by means of an initializer list:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::padded;

   const int array[8] = { 1, 2, 3, 4, 0, 0, 0, 0 };

   CustomVector<double,unaligned,padded> v( array, 4UL, 8UL );
   v = { 5, 6, 7 };
   \endcode

// The vector elements are assigned the values from the given initializer list. Missing values
// are reset to their default state. Note that in case the size of the initializer list exceeds
// the size of the vector, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( initializer_list<Type> list )
{
   if( list.size() > size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to custom vector" );
   }

   std::fill( std::copy( list.begin(), list.end(), v_ ), v_+capacity_, Type() );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array assignment to all vector elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid array size.
//
// This assignment operator offers the option to directly set all elements of the vector. The
// following example demonstrates this by means of an unaligned, padded custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::padded;

   const int array[8] = { 1, 2, 3, 4, 0, 0, 0, 0 };
   const int init[4]  = { 5, 6, 7 };

   CustomVector<double,unaligned,padded> v( array, 4UL, 8UL );
   v = init;
   \endcode

// The vector is assigned the values from the given static array. Missing values are initialized
// with default values (as e.g. the fourth element in the example). Note that the size of the
// static array must match the size of the custom vector. Otherwise a \a std::invalid_argument
// exception is thrown. Also note that after the assignment \a array will have the same entries
// as \a init.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the static array
        , size_t Dim >      // Dimension of the static array
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( const Other (&array)[Dim] )
{
   if( size_ != Dim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Array assignment to all vector elements.
//
// \param array The given std::array for the assignment.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid std::array size.
//
// This assignment operator offers the option to directly set all elements of the vector. The
// following example demonstrates this by means of an unaligned, padded custom vector:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::padded;

   const int array[8] = { 1, 2, 3, 4, 0, 0, 0, 0 };
   const std::array<int,4UL> init{ 5, 6, 7 };

   CustomVector<double,unaligned,padded> v( array, 4UL, 8UL );
   v = init;
   \endcode

// The vector is assigned the values from the given std::array. Missing values are initialized
// with default values (as e.g. the fourth element in the example). Note that the size of the
// std::array must match the size of the custom vector. Otherwise a \a std::invalid_argument
// exception is thrown. Also note that after the assignment \a array will have the same entries
// as \a init.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other    // Data type of the std::array
        , size_t Dim >      // Dimension of the std::array
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( const std::array<Other,Dim>& array )
{
   if( size_ != Dim ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid array size" );
   }

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for CustomVector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// The vector is initialized as a copy of the given vector. In case the current sizes of the two
// vectors don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( const CustomVector& rhs )
{
   if( rhs.size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   smpAssign( *this, *rhs );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Move assignment operator for CustomVector.
//
// \param rhs The vector to be moved into this instance.
// \return Reference to the assigned vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( CustomVector&& rhs ) noexcept
{
   size_     = rhs.size_;
   capacity_ = rhs.capacity_;
   v_        = rhs.v_;

   rhs.size_     = 0UL;
   rhs.capacity_ = 0UL;
   rhs.v_        = nullptr;

   BLAZE_INTERNAL_ASSERT( rhs.data() == nullptr, "Invalid data reference detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// The vector is initialized as a copy of the given vector. In case the current sizes of the two
// vectors don't match, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      smpAssign( *this, tmp );
   }
   else {
      if( IsSparseVector_v<VT> )
         reset();
      smpAssign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator+=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
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
/*!\brief Subtraction assignment operator for the subtraction of a vector
//        (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator-=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
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
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator*=( const Vector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( MultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( IsSparseVector_v<VT> || (*rhs).canAlias( this ) ) {
      const MultType tmp( *this * (*rhs) );
      if( IsSparseVector_v<MultType> )
         reset();
      smpAssign( *this, tmp );
   }
   else {
      smpMultAssign( *this, *rhs );
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator/=( const DenseVector<VT,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( DivType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( DivType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const DivType tmp( *this / (*rhs) );
      smpAssign( *this, tmp );
   }
   else {
      smpDivAssign( *this, *rhs );
   }

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
// \return Reference to the vector.
// \exception std::invalid_argument Invalid vector size for cross product.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side vector
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::operator%=( const Vector<VT,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size_ != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );
   assign( *this, tmp );

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
/*!\brief Returns the size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,padded,TF,Tag,RT>::size() const noexcept
{
   return size_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the vector.
//
// \return The minimum capacity of the vector.
//
// This function returns the minimum capacity of the vector, which corresponds to the current
// size plus padding.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,padded,TF,Tag,RT>::spacing() const noexcept
{
   return capacity_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the vector.
//
// \return The maximum capacity of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,padded,TF,Tag,RT>::capacity() const noexcept
{
   return capacity_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the vector.
//
// \return The number of non-zero elements in the vector.
//
// This function returns the number of non-zero elements in the vector (i.e. the elements that
// compare unequal to their default value). Note that the number of non-zero elements is always
// less than or equal to the current size of the vector.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline size_t CustomVector<Type,AF,padded,TF,Tag,RT>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( size_t i=0UL; i<size_; ++i ) {
      if( !isDefault<strict>( v_[i] ) )
         ++nonzeros;
   }

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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::reset()
{
   using blaze::clear;
   for( size_t i=0UL; i<size_; ++i )
      clear( v_[i] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the vector to its default state.
//
// \return void
//
// This function clears the vector to its default state. In case the vector has been passed the
// responsibility to manage the given array, it disposes the resource via the specified deleter.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::clear()
{
   size_     = 0UL;
   capacity_ = 0UL;
   v_        = nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Swapping the contents of two vectors.
//
// \param v The vector to be swapped.
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::swap( CustomVector& v ) noexcept
{
   using std::swap;

   swap( size_, v.size_ );
   swap( capacity_, v.capacity_ );
   swap( v_, v.v_ );
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
/*!\brief Scaling of the vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the vector.
//
// This function scales the vector by applying the given scalar value \a scalar to each element
// of the vector. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   using blaze::CustomVector;
   using blaze::unaliged;
   using blaze::padded;

   CustomVector<double,unaligned,padded> v( ... );

   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the scalar value
inline CustomVector<Type,AF,padded,TF,Tag,RT>&
   CustomVector<Type,AF,padded,TF,Tag,RT>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<size_; ++i )
      v_[i] *= scalar;
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
/*!\brief Resets the custom vector and replaces the array of elements with the given array.
//
// \param ptr The array of elements to be used by the vector.
// \param n The number of array elements to be used by the custom vector.
// \param nn The maximum size of the given array.
// \return void
// \exception std::invalid_argument Invalid setup of custom vector.
//
// This function resets the custom vector to the given array of elements of size \a n and capacity
// \a nn. The function fails if ...
//
//  - ... the passed pointer is \c nullptr;
//  - ... the specified capacity \a nn is insufficient for the given data type \a Type and
//    the available instruction set.
//
// In all failure cases a \a std::invalid_argument exception is thrown.
//
// \note In case a deleter was specified, the previously referenced array will only be destroyed
//       when the last custom vector referencing the array goes out of scope.
// \note The custom vector does NOT take responsibility for the new array of elements!
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::reset( Type* ptr, size_t n, size_t nn )
{
   CustomVector tmp( ptr, n, nn );
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
/*!\brief Returns whether the vector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomVector<Type,AF,padded,TF,Tag,RT>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the vector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address is aliased with the vector. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename Other >  // Data type of the foreign expression
inline bool CustomVector<Type,AF,padded,TF,Tag,RT>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the vector is properly aligned in memory.
//
// \return \a true in case the vector is aligned, \a false if not.
//
// This function returns whether the vector is guaranteed to be properly aligned in memory, i.e.
// whether the beginning and the end of the vector are guaranteed to conform to the alignment
// restrictions of the element type \a Type.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomVector<Type,AF,padded,TF,Tag,RT>::isAligned() const noexcept
{
   return ( AF || checkAlignment( v_ ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the vector can be used in SMP assignments.
//
// \return \a true in case the vector can be used in SMP assignments, \a false if not.
//
// This function returns whether the vector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// vector).
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool CustomVector<Type,AF,padded,TF,Tag,RT>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,padded,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,padded,TF,Tag,RT>::load( size_t index ) const noexcept
{
   if( AF )
      return loada( index );
   else
      return loadu( index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,padded,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,padded,TF,Tag,RT>::loada( size_t index ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( !AF || index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   return loada( v_+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE typename CustomVector<Type,AF,padded,TF,Tag,RT>::SIMDType
   CustomVector<Type,AF,padded,TF,Tag,RT>::loadu( size_t index ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );

   return loadu( v_+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense vector. The index
// must be smaller than the number of vector elements and it must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,padded,TF,Tag,RT>::store( size_t index, const SIMDType& value ) noexcept
{
   if( AF )
      storea( index, value );
   else
      storeu( index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense vector. The
// index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,padded,TF,Tag,RT>::storea( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( !AF || index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   storea( v_+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense vector.
// The index must be smaller than the number of vector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly! It
// is used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,padded,TF,Tag,RT>::storeu( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );

   storeu( v_+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the vector.
//
// \param index Access index. The index must be smaller than the number of vector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense vector. The index must be smaller than the number of vector elements and it must be
// a multiple of the number of values inside the SIMD element. This function must \b NOT be
// called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
BLAZE_ALWAYS_INLINE void
   CustomVector<Type,AF,padded,TF,Tag,RT>::stream( size_t index, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= capacity_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( v_+index ), "Invalid vector access index" );

   stream( v_+index, value );
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::assign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] = (*rhs)[i    ];
      v_[i+1UL] = (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] = (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::assign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   if( AF && useStreaming && size_ > ( cacheSize/( sizeof(Type) * 3UL ) ) && !(*rhs).isAliased( this ) )
   {
      size_t i( 0UL );

      for( ; i<ipos; i+=SIMDSIZE ) {
         stream( i, (*rhs).load(i) );
      }
      for( ; remainder && i<size_; ++i ) {
         v_[i] = (*rhs)[i];
      }
   }
   else
   {
      const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
      BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

      size_t i( 0UL );
      ConstIterator_t<VT> it( (*rhs).begin() );

      for( ; i<i4way; i+=SIMDSIZE*4UL ) {
         store( i             , it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE    , it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE*2UL, it.load() ); it += SIMDSIZE;
         store( i+SIMDSIZE*3UL, it.load() ); it += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
         store( i, it.load() );
      }
      for( ; remainder && i<size_; ++i, ++it ) {
         v_[i] = *it;
      }
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] = element->value();
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::addAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] += (*rhs)[i    ];
      v_[i+1UL] += (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] += (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::addAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) + it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) + it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) + it.load() );
   }
   for( ; remainder && i<size_; ++i, ++it ) {
      v_[i] += *it;
   }
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] += element->value();
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::subAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] -= (*rhs)[i    ];
      v_[i+1UL] -= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] -= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::subAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) - it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) - it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) - it.load() );
   }
   for( ; remainder && i<size_; ++i, ++it ) {
      v_[i] -= *it;
   }
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
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::multAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] *= (*rhs)[i    ];
      v_[i+1UL] *= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] *= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::multAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   constexpr bool remainder( !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) * it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) * it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) * it.load() );
   }
   for( ; remainder && i<size_; ++i, ++it ) {
      v_[i] *= *it;
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side sparse vector
inline void CustomVector<Type,AF,padded,TF,Tag,RT>::multAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const ResultType tmp( serial( *this ) );

   reset();

   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      v_[element->index()] = tmp[element->index()] * element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::divAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      v_[i    ] /= (*rhs)[i    ];
      v_[i+1UL] /= (*rhs)[i+1UL];
   }
   if( ipos < (*rhs).size() )
      v_[ipos] /= (*rhs)[ipos];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
template< typename VT >  // Type of the right-hand side dense vector
inline auto CustomVector<Type,AF,padded,TF,Tag,RT>::divAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   const size_t i4way( prevMultiple( size_, SIMDSIZE*4UL ) );
   BLAZE_INTERNAL_ASSERT( i4way <= size_, "Invalid end calculation" );

   size_t i( 0UL );
   ConstIterator_t<VT> it( (*rhs).begin() );

   for( ; i<i4way; i+=SIMDSIZE*4UL ) {
      store( i             , load(i             ) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE    , load(i+SIMDSIZE    ) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*2UL, load(i+SIMDSIZE*2UL) / it.load() ); it += SIMDSIZE;
      store( i+SIMDSIZE*3UL, load(i+SIMDSIZE*3UL) / it.load() ); it += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE, it+=SIMDSIZE ) {
      store( i, load(i) / it.load() );
   }
   for( ; i<size_; ++i, ++it ) {
      v_[i] /= *it;
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CUSTOMVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name CustomVector operators */
//@{
template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
void reset( CustomVector<Type,AF,PF,TF,Tag,RT>& v );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
void clear( CustomVector<Type,AF,PF,TF,Tag,RT>& v );

template< RelaxationFlag RF, typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
bool isDefault( const CustomVector<Type,AF,PF,TF,Tag,RT>& v );

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
bool isIntact( const CustomVector<Type,AF,PF,TF,Tag,RT>& v ) noexcept;

template< typename Type, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
void swap( CustomVector<Type,AF,PF,TF,Tag,RT>& a, CustomVector<Type,AF,PF,TF,Tag,RT>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given custom vector.
// \ingroup custom_vector
//
// \param v The custom vector to be resetted.
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void reset( CustomVector<Type,AF,PF,TF,Tag,RT>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given custom vector.
// \ingroup custom_vector
//
// \param v The custom vector to be cleared.
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void clear( CustomVector<Type,AF,PF,TF,Tag,RT>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given custom vector is in default state.
// \ingroup custom_vector
//
// \param v The custom vector to be tested for its default state.
// \return \a true in case the given vector is component-wise zero, \a false otherwise.
//
// This function checks whether the custom vector is in default state. For instance, in case
// the static vector is instantiated for a built-in integral or floating point data type, the
// function returns \a true in case all vector elements are 0 and \a false in case any vector
// element is not 0. Following example demonstrates the use of the \a isDefault function:

   \code
   using blaze::aligned;
   using blaze::padded;

   blaze::CustomVector<int,aligned,padded> a( ... );
   // ... Resizing and initialization
   if( isDefault( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( a ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the vector
        , AlignmentFlag AF   // Alignment flag
        , PaddingFlag PF     // Padding flag
        , bool TF            // Transpose flag
        , typename Tag       // Type tag
        , typename RT >      // Result type
inline bool isDefault( const CustomVector<Type,AF,PF,TF,Tag,RT>& v )
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given custom vector are intact.
// \ingroup custom_vector
//
// \param v The custom vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the custom vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::aligned;
   using blaze::padded;

   blaze::CustomVector<int,aligned,padded> a( ... );
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline bool isIntact( const CustomVector<Type,AF,PF,TF,Tag,RT>& v ) noexcept
{
   return ( v.size() <= v.capacity() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
// \ingroup custom_vector
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
*/
template< typename Type     // Data type of the vector
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , bool TF           // Transpose flag
        , typename Tag      // Type tag
        , typename RT >     // Result type
inline void swap( CustomVector<Type,AF,PF,TF,Tag,RT>& a, CustomVector<Type,AF,PF,TF,Tag,RT>& b ) noexcept
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
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
struct HasConstDataAccess< CustomVector<T,AF,PF,TF,Tag,RT> >
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
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
struct HasMutableDataAccess< CustomVector<T,AF,PF,TF,Tag,RT> >
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
template< typename T, AlignmentFlag AF, PaddingFlag PF, bool TF, typename Tag, typename RT >
struct IsCustom< CustomVector<T,AF,PF,TF,Tag,RT> >
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
template< typename T, PaddingFlag PF, bool TF, typename Tag, typename RT >
struct IsAligned< CustomVector<T,aligned,PF,TF,Tag,RT> >
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
template< typename T, PaddingFlag PF, bool TF, typename Tag, typename RT >
struct IsContiguous< CustomVector<T,aligned,PF,TF,Tag,RT> >
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
template< typename T, AlignmentFlag AF, bool TF, typename Tag, typename RT >
struct IsPadded< CustomVector<T,AF,padded,TF,Tag,RT> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
