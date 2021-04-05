//=================================================================================================
/*!
//  \file blaze/math/dense/HybridVector.h
//  \brief Header file for the HybridVector class template
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//
//  * The names of its contributors may not be used to endorse or promote products derived
//    from this software without specific prior written permission.
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

#ifndef _BLAZE_MATH_DENSE_HYBRIDVECTOR_H_
#define _BLAZE_MATH_DENSE_HYBRIDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

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
#include <blaze/math/Infinity.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/PaddingFlag.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/NextMultiple.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SolveTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/AlignedArray.h>
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
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blaze/util/typetraits/IsNumeric.h>
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
/*!\defgroup hybrid_vector HybridVector
// \ingroup dense_vector
*/
/*!\brief Efficient implementation of a dynamically sized vector with static memory.
// \ingroup hybrid_vector
//
// The HybridVector class template combines the flexibility of a dynamically sized vector with
// the efficiency and performance of a fixed-size vector. It is implemented as a crossing between
// the blaze::StaticVector and the blaze::DynamicVector class templates: Similar to the static
// vector it uses static stack memory instead of dynamically allocated memory and similar to the
// dynamic vector it can be resized (within the extend of the static memory). The type of the
// elements, the maximum number of elements, the transpose flag, the alignment, the padding, and
// the group tag of the vector can be specified via the six template parameters:

   \code
   namespace blaze {

   template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
   class HybridVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. HybridVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - N   : specifies the maximum number of vector elements, i.e. the maximum size of the vector.
//          It is expected that HybridVector is only used for tiny and small vectors.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - AF  : specifies whether the first element of the vector is properly aligned with respect to
//          the available instruction set (SSE, AVX, ...). Possible values are \c blaze::aligned
//          and \c blaze::unaligned. The default value is \c blaze::defaultAlignmentFlag.
//  - PF  : specifies whether the vector should be padded to maximize the efficiency of vectorized
//          operations. Possible values are \c blaze::padded and \c blaze::unpadded. The default
//          value is \c blaze::defaultPaddingFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// These contiguously stored elements can be directly accessed with the subscript operator. The
// numbering of the vector elements is

                             \f[\left(\begin{array}{*{4}{c}}
                             0 & 1 & \cdots & N-1 \\
                             \end{array}\right)\f]

// The use of HybridVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// vectors with fitting element types. The following example gives an impression of the use of a
// 2-dimensional HybridVector:

   \code
   using blaze::HybridVector;
   using blaze::CompressedVector;
   using blaze::StaticMatrix;

   HybridVector<double,2UL> a( 2 );  // Non-initialized 2D vector of size 2
   a[0] = 1.0;                       // Initialization of the first element
   a[1] = 2.0;                       // Initialization of the second element

   HybridVector<double,2UL> b( 2, 2.0 );  // Directly, homogeneously initialized 2D vector
   CompressedVector<float>  c( 2 );       // Empty sparse single precision vector
   HybridVector<double,2UL> d;            // Default constructed hybrid vector
   StaticMatrix<double,2UL,2UL> A;        // Default constructed static row-major matrix

   d = a + b;  // Vector addition between vectors of equal element type
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
class HybridVector
   : public DenseVector< HybridVector<Type,N,TF,AF,PF,Tag>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This       = HybridVector<Type,N,TF,AF,PF,Tag>;  //!< Type of this HybridVector instance.
   using BaseType   = DenseVector<This,TF>;               //!< Base type of this HybridVector instance.
   using ResultType = This;                               //!< Result type for expression template evaluations.

   //! Transpose type for expression template evaluations.
   using TransposeType = HybridVector<Type,N,!TF,AF,PF,Tag>;

   using ElementType   = Type;                      //!< Type of the vector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations.
   using TagType       = Tag;                       //!< Tag type of this HybridVector instance.
   using CompositeType = const HybridVector&;       //!< Data type for composite expression templates.

   using Reference      = Type&;        //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = Type*;        //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = DenseIterator<Type,AF>;        //!< Iterator over non-constant elements.
   using ConstIterator = DenseIterator<const Type,AF>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a HybridVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = HybridVector<NewType,N,TF,AF,PF,Tag>;  //!< The type of the other HybridVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a HybridVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = HybridVector<Type,NewN,TF,AF,PF,Tag>;  //!< The type of the other HybridVector.
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
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
            constexpr HybridVector();
   explicit constexpr HybridVector( size_t n );
            inline    HybridVector( size_t n, const Type& init );
            constexpr HybridVector( initializer_list<Type> list );

   template< typename Other >
   inline HybridVector( size_t n, const Other* array );

   template< typename Other, size_t Dim >
   constexpr HybridVector( const Other (&array)[Dim] );

   template< typename Other, size_t Dim >
   constexpr HybridVector( const std::array<Other,Dim>& array );

   constexpr HybridVector( const HybridVector& v );

   template< typename VT >
   inline HybridVector( const Vector<VT,TF>& v );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~HybridVector() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr Reference      operator[]( size_t index ) noexcept;
   constexpr ConstReference operator[]( size_t index ) const noexcept;
   inline    Reference      at( size_t index );
   inline    ConstReference at( size_t index ) const;
   constexpr Pointer        data  () noexcept;
   constexpr ConstPointer   data  () const noexcept;
   constexpr Iterator       begin () noexcept;
   constexpr ConstIterator  begin () const noexcept;
   constexpr ConstIterator  cbegin() const noexcept;
   constexpr Iterator       end   () noexcept;
   constexpr ConstIterator  end   () const noexcept;
   constexpr ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   constexpr HybridVector& operator=( const Type& rhs ) &;
   constexpr HybridVector& operator=( initializer_list<Type> list ) &;

   template< typename Other, size_t Dim >
   constexpr HybridVector& operator=( const Other (&array)[Dim] ) &;

   template< typename Other, size_t Dim >
   constexpr HybridVector& operator=( const std::array<Other,Dim>& array ) &;

   constexpr HybridVector& operator=( const HybridVector& rhs ) &;

   template< typename VT > inline HybridVector& operator= ( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline HybridVector& operator+=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline HybridVector& operator-=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline HybridVector& operator*=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline HybridVector& operator/=( const DenseVector<VT,TF>& rhs ) &;
   template< typename VT > inline HybridVector& operator%=( const Vector<VT,TF>& rhs ) &;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
          constexpr size_t size() const noexcept;
   static constexpr size_t spacing() noexcept;
   static constexpr size_t capacity() noexcept;
          inline    size_t nonZeros() const;
          constexpr void   reset();
          constexpr void   clear();
          constexpr void   resize( size_t n, bool preserve=true );
          constexpr void   extend( size_t n, bool preserve=true );
          inline    void   swap( HybridVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline HybridVector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Memory functions****************************************************************************
   /*!\name Memory functions */
   //@{
   static inline void* operator new  ( std::size_t size );
   static inline void* operator new[]( std::size_t size );
   static inline void* operator new  ( std::size_t size, const std::nothrow_t& );
   static inline void* operator new[]( std::size_t size, const std::nothrow_t& );

   static inline void operator delete  ( void* ptr );
   static inline void operator delete[]( void* ptr );
   static inline void operator delete  ( void* ptr, const std::nothrow_t& );
   static inline void operator delete[]( void* ptr, const std::nothrow_t& );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! The number of elements packed within a single SIMD vector.
   static constexpr size_t SIMDSIZE = SIMDTrait<Type>::size;

   //! Alignment adjustment.
   static constexpr size_t NN = ( PF == padded ? nextMultiple( N, SIMDSIZE ) : N );
   //**********************************************************************************************

   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename VT >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        NN >= SIMDSIZE &&
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

 public:
   //**Debugging functions*************************************************************************
   /*!\name Debugging functions */
   //@{
   constexpr bool isIntact() const noexcept;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   static constexpr bool isAligned() noexcept;

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
   //**********************************************************************************************
   //! Alignment of the data elements.
   static constexpr size_t Alignment =
      ( AF == aligned ? AlignmentOf_v<Type> : std::alignment_of<Type>::value );

   //! Type of the aligned storage.
   using AlignedStorage = AlignedArray<Type,NN,Alignment>;
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   AlignedStorage v_;  //!< The statically allocated vector elements.
                       /*!< Access to the vector values is gained via the subscript operator.
                            The order of the elements is
                            \f[\left(\begin{array}{*{4}{c}}
                            0 & 1 & \cdots & N-1 \\
                            \end{array}\right)\f] */
   size_t size_;       //!< The current size/dimension of the vector.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_STATIC_ASSERT( PF == unpadded || NN % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( NN >= N );
   BLAZE_STATIC_ASSERT( IsVectorizable_v<Type> || NN == N );
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

template< typename Type, typename... Ts >
HybridVector( Type, Ts... ) -> HybridVector<Type,1+sizeof...(Ts)>;

template< typename Type, size_t N >
HybridVector( Type (&)[N] ) -> HybridVector< RemoveCV_t<Type>, N >;

template< typename Type, size_t N >
HybridVector( std::array<Type,N> ) -> HybridVector<Type,N>;

#endif
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for HybridVector.
//
// The size of a default constructed HybridVector is initially set to 0.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector()
   : v_   ()       // The statically allocated vector elements
   , size_( 0UL )  // The current size/dimension of the vector
{
   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n.
//
// \param n The size of the vector.
// \exception std::invalid_argument Invalid size for hybrid vector.
//
// This constructor creates a hybrid vector of size \a n and initializes all vector elements to
// the default value (for instance 0 for integral types). In case \a n is larger than the maximum
// allowed number of elements (i.e. \a n > N) a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( size_t n )
   : v_   ()     // The statically allocated vector elements
   , size_( n )  // The current size/dimension of the vector
{
   if( n > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid size for hybrid vector" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogeneous initialization of all \a n vector elements.
//
// \param n The size of the vector.
// \param init The initial value of the vector elements.
// \exception std::invalid_argument Invalid size for hybrid vector.
//
// This constructor creates a hybrid vector of size \a n and initializes all vector elements with
// the specified value. In case \a n is larger than the maximum allowed number of elements (i.e.
// \a n > N) a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( size_t n, const Type& init )
   : size_( n )  // The current size/dimension of the vector
   // v_ is intentionally left uninitialized
{
   using blaze::clear;

   if( n > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid size for hybrid vector" );
   }

   for( size_t i=0UL; i<n; ++i )
      v_[i] = init;

   if( IsNumeric_v<Type> ) {
      for( size_t i=n; i<NN; ++i )
         clear( v_[i] );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all vector elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid setup of hybrid vector.
//
// This constructor provides the option to explicitly initialize the elements of the vector by
// means of an initializer list:

   \code
   blaze::HybridVector<double,6UL> v1{ 4.2, 6.3, -1.2 };
   \endcode

// The vector is sized according to the size of the initializer list and all its elements are
// (copy) assigned the values of the given initializer list. In case the size of the given list
// exceeds the maximum size of the hybrid vector (i.e. is larger than \a N), a
// \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( initializer_list<Type> list )
   : v_   ()               // The statically allocated vector elements
   , size_( list.size() )  // The current size/dimension of the vector
{
   if( size_ > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of hybrid vector" );
   }

   size_t i( 0UL );

   for( const auto& element : list ) {
      v_[i] = element;
      ++i;
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all vector elements.
//
// \param n The size of the vector.
// \param array Dynamic array for the initialization.
// \exception std::invalid_argument Invalid size for hybrid vector.
//
// This assignment operator offers the option to directly initialize the elements of the vector
// with a dynamic array:

   \code
   double* array = new double[6];
   // ... Initialization of the dynamic array
   blaze::HybridVector<double,6> v( 6UL, array );
   delete[] array;
   \endcode

// The vector is sized according to the size of the array and initialized with the values from
// the given array. In case the size of the given array exceeds the maximum size of the hybrid
// vector (i.e. is larger than \a N), a \a std::invalid_argument exception is thrown.\n
// Note that it is expected that the given \a array has at least \a n elements. Providing an
// array with less elements results in undefined behavior!
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the initialization array
inline HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( size_t n, const Other* array )
   : size_( n )  // The current size/dimension of the vector
   // v_ is intentionally left uninitialized
{
   using blaze::clear;

   if( n > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of hybrid vector" );
   }

   for( size_t i=0UL; i<n; ++i )
      v_[i] = array[i];

   if( IsNumeric_v<Type> ) {
      for( size_t i=n; i<NN; ++i )
         clear( v_[i] );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array initialization of all vector elements.
//
// \param array Static array for the initialization.
//
// This assignment operator offers the option to directly initialize the elements of the vector
// with a static array:

   \code
   const double init[2] = { 1.0, 2.0 };
   blaze::HybridVector<double,4> v( init );
   \endcode

// The vector is sized according to the size of the static array and initialized with the values
// from the given static array. This constructor only works for arrays with a size smaller-or-equal
// than the maximum number of elements of the hybrid vector (i.e. Dim <= N). The attempt to use a
// larger static array will result in a compile time error.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other    // Data type of the static array
        , size_t Dim >      // Number of elements of the static array
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( const Other (&array)[Dim] )
   : v_   ( array )  // The statically allocated vector elements
   , size_( Dim )    // The current size/dimension of the vector
{
   BLAZE_STATIC_ASSERT( Dim <= N );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Initialization of all vector elements from the given std::array.
//
// \param array The given std::array for the initialization.
//
// This constructor operator offers the option to directly initialize the elements of the vector
// with a std::array:

   \code
   std::array<double,2> init{ 1.0, 2.0 };
   blaze::HybridVector<double,4> v( init );
   \endcode

// The vector is sized according to the size of the std::array and initialized with the values
// from the given std::array. This constructor only works for arrays with a size smaller-or-equal
// than the maximum number of elements of the hybrid vector (i.e. Dim <= N). The attempt to use a
// larger std::array will result in a compile time error.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other    // Data type of the std::array
        , size_t Dim >      // Dimension of the std::array
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( const std::array<Other,Dim>& array )
   : v_   ( array )  // The statically allocated vector elements
   , size_( Dim )    // The current size/dimension of the vector
{
   BLAZE_STATIC_ASSERT( Dim <= N );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for HybridVector.
//
// \param v Vector to be copied.
//
// The copy constructor is explicitly defined in order to enable/facilitate NRV optimization.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( const HybridVector& v )
   : BaseType()        // Initialization of the base class
   , v_   ( v.v_ )     // The statically allocated vector elements
   , size_( v.size_ )  // The current size/dimension of the vector
{
   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different vectors.
//
// \param v Vector to be copied.
// \exception std::invalid_argument Invalid setup of hybrid vector.
//
// This constructor initializes the hybrid vector from the given vector. In case the size of
// the given vector exceeds the maximum size of the hybrid vector (i.e. is larger than \a N),
// a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >  // Type of the foreign vector
inline HybridVector<Type,N,TF,AF,PF,Tag>::HybridVector( const Vector<VT,TF>& v )
   : size_( (*v).size() )  // The current size/dimension of the vector
   // v_ is intentionally left uninitialized
{
   using blaze::assign;
   using blaze::clear;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*v).size() > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of hybrid vector" );
   }

   for( size_t i=( IsSparseVector_v<VT> ? 0UL : size_ );
               i<( IsNumeric_v<Type>    ? NN  : size_ ); ++i ) {
      clear( v_[i] );
   }

   assign( *this, *v );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::Reference
   HybridVector<Type,N,TF,AF,PF,Tag>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   return v_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference-to-const to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstReference
   HybridVector<Type,N,TF,AF,PF,Tag>::operator[]( size_t index ) const noexcept
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline typename HybridVector<Type,N,TF,AF,PF,Tag>::Reference
   HybridVector<Type,N,TF,AF,PF,Tag>::at( size_t index )
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstReference
   HybridVector<Type,N,TF,AF,PF,Tag>::at( size_t index ) const
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
// This function returns a pointer to the internal storage of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::Pointer
   HybridVector<Type,N,TF,AF,PF,Tag>::data() noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the vector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstPointer
   HybridVector<Type,N,TF,AF,PF,Tag>::data() const noexcept
{
   return v_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the hybrid vector.
//
// \return Iterator to the first element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::Iterator
   HybridVector<Type,N,TF,AF,PF,Tag>::begin() noexcept
{
   return Iterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the hybrid vector.
//
// \return Iterator to the first element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstIterator
   HybridVector<Type,N,TF,AF,PF,Tag>::begin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the hybrid vector.
//
// \return Iterator to the first element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstIterator
   HybridVector<Type,N,TF,AF,PF,Tag>::cbegin() const noexcept
{
   return ConstIterator( v_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the hybrid vector.
//
// \return Iterator just past the last element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::Iterator
   HybridVector<Type,N,TF,AF,PF,Tag>::end() noexcept
{
   BLAZE_INTERNAL_ASSERT( size_ <= N, "Invalid size detected" );
   return Iterator( v_ + size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the hybrid vector.
//
// \return Iterator just past the last element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstIterator
   HybridVector<Type,N,TF,AF,PF,Tag>::end() const noexcept
{
   BLAZE_INTERNAL_ASSERT( size_ <= N, "Invalid size detected" );
   return ConstIterator( v_ + size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the hybrid vector.
//
// \return Iterator just past the last element of the hybrid vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr typename HybridVector<Type,N,TF,AF,PF,Tag>::ConstIterator
   HybridVector<Type,N,TF,AF,PF,Tag>::cend() const noexcept
{
   BLAZE_INTERNAL_ASSERT( size_ <= N, "Invalid size detected" );
   return ConstIterator( v_ + size_ );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( const Type& rhs ) &
{
   BLAZE_INTERNAL_ASSERT( size_ <= N, "Invalid size detected" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] = rhs;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
// \exception Invalid assignment to hybrid vector.
//
// This assignment operator offers the option to directly assign to all elements of the vector
// by means of an initializer list:

   \code
   blaze::HybridVector<double,6UL> v;
   v = { 4.2, 6.3, -1.2 };
   \endcode

// The vector is resized according to the size of the initializer list and all its elements
// are (copy) assigned the values from the given initializer list. In case the size of the
// given list exceeds the maximum size of the hybrid vector (i.e. is larger than \a N), a
// \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( initializer_list<Type> list ) &
{
   if( list.size() > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to hybrid vector" );
   }

   resize( list.size(), false );

   size_t i( 0UL );

   for( const auto& element : list ) {
      v_[i] = element;
      ++i;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Array assignment to all vector elements.
//
// \param array Static array for the assignment.
// \return Reference to the assigned vector.
//
// This assignment operator offers the option to directly set all elements of the vector:

   \code
   const double init[2] = { 1.0, 2.0 };
   blaze::HybridVector<double,4> v;
   v = init;
   \endcode

// The vector is sized according to the size of the static array and assigned the values of the
// given static array. This assignment operator only works for arrays with a size smaller-or-equal
// than the maximum number of elements of the hybrid vector. (i.e. Dim <= N). The attempt to use a
// larger static array will result in a compile time error.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other    // Data type of the static array
        , size_t Dim >      // Number of elements of the static array
constexpr HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( const Other (&array)[Dim] ) &
{
   BLAZE_STATIC_ASSERT( Dim <= N );

   resize( Dim, false );

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
//
// This assignment operator offers the option to directly set all elements of the vector:

   \code
   const std::array<double,2> init{ 1.0, 2.0 };
   blaze::HybridVector<double,4> v;
   v = init;
   \endcode

// The vector is sized according to the size of the std::array and assigned the values of the
// given std::array. This assignment operator only works for arrays with a size smaller-or-equal
// than the maximum number of elements of the hybrid vector. (i.e. Dim <= N). The attempt to use
// a larger std::array will result in a compile time error.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other    // Data type of the std::array
        , size_t Dim >      // Dimension of the std::array
constexpr HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( const std::array<Other,Dim>& array ) &
{
   BLAZE_STATIC_ASSERT( Dim <= N );

   resize( Dim, false );

   for( size_t i=0UL; i<Dim; ++i )
      v_[i] = array[i];

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for HybridVector.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// Explicit definition of a copy assignment operator for performance reasons.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( const HybridVector& rhs ) &
{
   using blaze::assign;

   BLAZE_INTERNAL_ASSERT( size_ <= N, "Invalid size detected" );

   resize( rhs.size() );
   assign( *this, *rhs );

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid assignment to hybrid vector.
//
// This constructor initializes the vector as a copy of the given vector. In case the size
// of the given vector is larger than \a N, a \a std::invalid_argument exception is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >  // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator=( const Vector<VT,TF>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to hybrid vector" );
   }

   if( (*rhs).canAlias( this ) ) {
      HybridVector tmp( *rhs );
      swap( tmp );
   }
   else {
      resize( (*rhs).size(), false );
      if( IsSparseVector_v<VT> )
         reset();
      assign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator+=( const Vector<VT,TF>& rhs ) &
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      addAssign( *this, tmp );
   }
   else {
      addAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator-=( const Vector<VT,TF>& rhs ) &
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const ResultType_t<VT> tmp( *rhs );
      subAssign( *this, tmp );
   }
   else {
      subAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator*=( const Vector<VT,TF>& rhs ) &
{
   using blaze::assign;
   using blaze::multAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( IsSparseVector_v<VT> || (*rhs).canAlias( this ) ) {
      const HybridVector tmp( *this * (*rhs) );
      assign( *this, tmp );
   }
   else {
      multAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator/=( const DenseVector<VT,TF>& rhs ) &
{
   using blaze::assign;
   using blaze::divAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      const HybridVector tmp( *this / (*rhs) );
      assign( *this, tmp );
   }
   else {
      divAssign( *this, *rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact(), "Invariant violation detected" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side vector
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::operator%=( const Vector<VT,TF>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< This, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size_ != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );
   assign( *this, tmp );

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
/*!\brief Returns the current size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr size_t HybridVector<Type,N,TF,AF,PF,Tag>::size() const noexcept
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr size_t HybridVector<Type,N,TF,AF,PF,Tag>::spacing() noexcept
{
   return NN;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
//
// \return The maximum capacity of the vector.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr size_t HybridVector<Type,N,TF,AF,PF,Tag>::capacity() noexcept
{
   return NN;
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline size_t HybridVector<Type,N,TF,AF,PF,Tag>::nonZeros() const
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void HybridVector<Type,N,TF,AF,PF,Tag>::reset()
{
   using blaze::clear;
   for( size_t i=0UL; i<size_; ++i )
      clear( v_[i] );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the vector.
//
// \return void
//
// After the clear() function, the size of the vector is 0.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void HybridVector<Type,N,TF,AF,PF,Tag>::clear()
{
   resize( 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
// \exception std::invalid_argument Invalid size for hybrid vector.
//
// This function resizes the vector to the given size \a n. In case the given size \a n is larger
// than the maximum number of vector elements (i.e. if n > N) a \a std::invalid_argument exception
// is thrown. Note that this function may invalidate all existing views (subvectors, ...) on the
// vector if it used to shrink the vector. Additionally, during this operation all vector elements
// are potentially changed. In order to preserve the old vector values, the \a preserve flag can be
// set to \a true.
//
// Note that in case the size of the vector is increased new vector elements are not initialized!
// This is illustrated by the following example, which demonstrates the resizing of a vector of
// size 2 to a vector of size 4. The new, uninitialized elements are marked with \a x:

                              \f[
                              \left(\begin{array}{*{2}{c}}
                              1 & 2 \\
                              \end{array}\right)

                              \Longrightarrow

                              \left(\begin{array}{*{4}{c}}
                              1 & 2 & x & x \\
                              \end{array}\right)
                              \f]
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void HybridVector<Type,N,TF,AF,PF,Tag>::resize( size_t n, bool preserve )
{
   using blaze::clear;

   MAYBE_UNUSED( preserve );

   if( n > N ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid size for hybrid vector" );
   }

   if( IsNumeric_v<Type> && n < size_ ) {
      for( size_t i=n; i<size_; ++i )
         clear( v_[i] );
   }

   size_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Extending the size of the vector.
//
// \param n Number of additional vector elements.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function increases the vector size by \a n elements. In case the resulting size
// of the vector is larger than the maximum number of vector elements (i.e. if n > N) a
// \a std::invalid_argument exception is thrown. During this operation, all vector elements
// are potentially changed. In order to preserve the old vector values, the \a preserve flag
// can be set to \a true.\n
// Note that new vector elements are not initialized!
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void HybridVector<Type,N,TF,AF,PF,Tag>::extend( size_t n, bool preserve )
{
   MAYBE_UNUSED( preserve );
   resize( size_+n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two hybrid vectors.
//
// \param v The vector to be swapped.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void HybridVector<Type,N,TF,AF,PF,Tag>::swap( HybridVector& v ) noexcept
{
   using std::swap;

   const size_t maxsize( max( size_, v.size_ ) );
   for( size_t i=0UL; i<maxsize; ++i )
      swap( v_[i], v.v_[i] );
   swap( size_, v.size_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scaling of the vector by the scalar value \a scalar (\f$ \vec{a}*=s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the vector.
//
// This function scales the vector by applying the given scalar value \a scalar to each element
// of the vector. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::HybridVector<int,3> a;
   // ... Resizing and initialization
   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline HybridVector<Type,N,TF,AF,PF,Tag>&
   HybridVector<Type,N,TF,AF,PF,Tag>::scale( const Other& scalar )
{
   for( size_t i=0; i<size_; ++i )
      v_[i] *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  MEMORY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Class specific implementation of operator new.
//
// \param size The total number of bytes to be allocated.
// \return Pointer to the newly allocated memory.
// \exception std::bad_alloc Allocation failed.
//
// This class-specific implementation of operator new provides the functionality to allocate
// dynamic memory based on the alignment restrictions of the StaticVector class template.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void* HybridVector<Type,N,TF,AF,PF,Tag>::operator new( std::size_t size )
{
   MAYBE_UNUSED( size );

   BLAZE_INTERNAL_ASSERT( size == sizeof( HybridVector ), "Invalid number of bytes detected" );

   return allocate<HybridVector>( 1UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of operator new[].
//
// \param size The total number of bytes to be allocated.
// \return Pointer to the newly allocated memory.
// \exception std::bad_alloc Allocation failed.
//
// This class-specific implementation of operator new provides the functionality to allocate
// dynamic memory based on the alignment restrictions of the HybridVector class template.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void* HybridVector<Type,N,TF,AF,PF,Tag>::operator new[]( std::size_t size )
{
   BLAZE_INTERNAL_ASSERT( size >= sizeof( HybridVector )       , "Invalid number of bytes detected" );
   BLAZE_INTERNAL_ASSERT( size %  sizeof( HybridVector ) == 0UL, "Invalid number of bytes detected" );

   return allocate<HybridVector>( size/sizeof(HybridVector) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of the no-throw operator new.
//
// \param size The total number of bytes to be allocated.
// \return Pointer to the newly allocated memory.
// \exception std::bad_alloc Allocation failed.
//
// This class-specific implementation of operator new provides the functionality to allocate
// dynamic memory based on the alignment restrictions of the HybridVector class template.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void* HybridVector<Type,N,TF,AF,PF,Tag>::operator new( std::size_t size, const std::nothrow_t& )
{
   MAYBE_UNUSED( size );

   BLAZE_INTERNAL_ASSERT( size == sizeof( HybridVector ), "Invalid number of bytes detected" );

   return allocate<HybridVector>( 1UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of the no-throw operator new[].
//
// \param size The total number of bytes to be allocated.
// \return Pointer to the newly allocated memory.
// \exception std::bad_alloc Allocation failed.
//
// This class-specific implementation of operator new provides the functionality to allocate
// dynamic memory based on the alignment restrictions of the HybridVector class template.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void* HybridVector<Type,N,TF,AF,PF,Tag>::operator new[]( std::size_t size, const std::nothrow_t& )
{
   BLAZE_INTERNAL_ASSERT( size >= sizeof( HybridVector )       , "Invalid number of bytes detected" );
   BLAZE_INTERNAL_ASSERT( size %  sizeof( HybridVector ) == 0UL, "Invalid number of bytes detected" );

   return allocate<HybridVector>( size/sizeof(HybridVector) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of operator delete.
//
// \param ptr The memory to be deallocated.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void HybridVector<Type,N,TF,AF,PF,Tag>::operator delete( void* ptr )
{
   deallocate( static_cast<HybridVector*>( ptr ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of operator delete[].
//
// \param ptr The memory to be deallocated.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void HybridVector<Type,N,TF,AF,PF,Tag>::operator delete[]( void* ptr )
{
   deallocate( static_cast<HybridVector*>( ptr ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of no-throw operator delete.
//
// \param ptr The memory to be deallocated.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void HybridVector<Type,N,TF,AF,PF,Tag>::operator delete( void* ptr, const std::nothrow_t& )
{
   deallocate( static_cast<HybridVector*>( ptr ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Class specific implementation of no-throw operator delete[].
//
// \param ptr The memory to be deallocated.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void HybridVector<Type,N,TF,AF,PF,Tag>::operator delete[]( void* ptr, const std::nothrow_t& )
{
   deallocate( static_cast<HybridVector*>( ptr ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  DEBUGGING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the invariants of the hybrid vector are intact.
//
// \return \a true in case the hybrid vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the hybrid vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false.
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr bool HybridVector<Type,N,TF,AF,PF,Tag>::isIntact() const noexcept
{
   if( size_ > N )
      return false;

   if( IsVectorizable_v<Type> ) {
      for( size_t i=size_; i<NN; ++i ) {
         if( !isDefault<strict>( v_[i] ) )
            return false;
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool HybridVector<Type,N,TF,AF,PF,Tag>::canAlias( const Other* alias ) const noexcept
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool HybridVector<Type,N,TF,AF,PF,Tag>::isAliased( const Other* alias ) const noexcept
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr bool HybridVector<Type,N,TF,AF,PF,Tag>::isAligned() noexcept
{
   return AF == aligned;
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE typename HybridVector<Type,N,TF,AF,PF,Tag>::SIMDType
   HybridVector<Type,N,TF,AF,PF,Tag>::load( size_t index ) const noexcept
{
   if( AF == aligned )
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE typename HybridVector<Type,N,TF,AF,PF,Tag>::SIMDType
   HybridVector<Type,N,TF,AF,PF,Tag>::loada( size_t index ) const noexcept
{
   using blaze::loada;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= NN, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( &v_[index] ), "Invalid alignment detected" );

   return loada( &v_[index] );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE typename HybridVector<Type,N,TF,AF,PF,Tag>::SIMDType
   HybridVector<Type,N,TF,AF,PF,Tag>::loadu( size_t index ) const noexcept
{
   using blaze::loadu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= NN, "Invalid vector access index" );

   return loadu( &v_[index] );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE void
   HybridVector<Type,N,TF,AF,PF,Tag>::store( size_t index, const SIMDType& value ) noexcept
{
   if( AF == aligned )
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE void
   HybridVector<Type,N,TF,AF,PF,Tag>::storea( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storea;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= NN, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( &v_[index] ), "Invalid alignment detected" );

   storea( &v_[index], value );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE void
   HybridVector<Type,N,TF,AF,PF,Tag>::storeu( size_t index, const SIMDType& value ) noexcept
{
   using blaze::storeu;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= NN, "Invalid vector access index" );

   storeu( &v_[index], value );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
BLAZE_ALWAYS_INLINE void
   HybridVector<Type,N,TF,AF,PF,Tag>::stream( size_t index, const SIMDType& value ) noexcept
{
   using blaze::stream;

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= NN, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( checkAlignment( &v_[index] ), "Invalid alignment detected" );

   stream( &v_[index], value );
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::assign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] = (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::assign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   constexpr bool remainder( PF == unpadded || !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; i<ipos; i+=SIMDSIZE ) {
      store( i, (*rhs).load(i) );
   }
   for( ; remainder && i<size_; ++i ) {
      v_[i] = (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side sparse vector
inline void HybridVector<Type,N,TF,AF,PF,Tag>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::addAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] += (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::addAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   constexpr bool remainder( PF == unpadded || !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; i<ipos; i+=SIMDSIZE ) {
      store( i, load(i) + (*rhs).load(i) );
   }
   for( ; remainder && i<size_; ++i ) {
      v_[i] += (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side sparse vector
inline void HybridVector<Type,N,TF,AF,PF,Tag>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::subAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] -= (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::subAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   constexpr bool remainder( PF == unpadded || !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; i<ipos; i+=SIMDSIZE ) {
      store( i, load(i) - (*rhs).load(i) );
   }
   for( ; remainder && i<size_; ++i ) {
      v_[i] -= (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side sparse vector
inline void HybridVector<Type,N,TF,AF,PF,Tag>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::multAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] *= (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::multAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   constexpr bool remainder( PF == unpadded || !IsPadded_v<VT> );

   const size_t ipos( remainder ? prevMultiple( size_, SIMDSIZE ) : size_ );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; i<ipos; i+=SIMDSIZE ) {
      store( i, load(i) * (*rhs).load(i) );
   }
   for( ; remainder && i<size_; ++i ) {
      v_[i] *= (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side sparse vector
inline void HybridVector<Type,N,TF,AF,PF,Tag>::multAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   const HybridVector tmp( serial( *this ) );

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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::divAssign( const DenseVector<VT,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   for( size_t i=0UL; i<size_; ++i )
      v_[i] /= (*rhs)[i];
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
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
template< typename VT >     // Type of the right-hand side dense vector
inline auto HybridVector<Type,N,TF,AF,PF,Tag>::divAssign( const DenseVector<VT,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( (*rhs).size() == size_, "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size_, SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size_, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; i<ipos; i+=SIMDSIZE ) {
      store( i, load(i) / (*rhs).load(i) );
   }
   for( ; i<size_; ++i ) {
      v_[i] /= (*rhs)[i];
   }
}
//*************************************************************************************************








//=================================================================================================
//
//  HYBRIDVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name HybridVector operators */
//@{
template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
constexpr void reset( HybridVector<Type,N,TF,AF,PF,Tag>& v );

template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
constexpr void clear( HybridVector<Type,N,TF,AF,PF,Tag>& v );

template< RelaxationFlag RF, typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
bool isDefault( const HybridVector<Type,N,TF,AF,PF,Tag>& v );

template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
constexpr bool isIntact( const HybridVector<Type,N,TF,AF,PF,Tag>& v );

template< typename Type, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
void swap( HybridVector<Type,N,TF,AF,PF,Tag>& a, HybridVector<Type,N,TF,AF,PF,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given hybrid vector.
// \ingroup hybrid_vector
//
// \param v The vector to be resetted.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void reset( HybridVector<Type,N,TF,AF,PF,Tag>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given hybrid vector.
// \ingroup hybrid_vector
//
// \param v The vector to be cleared.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr void clear( HybridVector<Type,N,TF,AF,PF,Tag>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given hybrid vector is in default state.
// \ingroup hybrid_vector
//
// \param v The hybrid vector to be tested for its default state.
// \return \a true in case the given vector's size is zero, \a false otherwise.
//
// This function checks whether the hybrid vector is in default (constructed) state, i.e. if
// it's size is 0. In case it is in default state, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isDefault() function:

   \code
   blaze::HybridVector<double,3> a;
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
        , size_t N           // Number of elements
        , bool TF            // Transpose flag
        , AlignmentFlag AF   // Alignment flag
        , PaddingFlag PF     // Padding flag
        , typename Tag >     // Type tag
inline bool isDefault( const HybridVector<Type,N,TF,AF,PF,Tag>& v )
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given hybrid vector are intact.
// \ingroup hybrid_vector
//
// \param v The hybrid vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the hybrid vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::HybridVector<double,3> a;
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
constexpr bool isIntact( const HybridVector<Type,N,TF,AF,PF,Tag>& v )
{
   return v.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two hybrid vectors.
// \ingroup hybrid_vector
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
*/
template< typename Type     // Data type of the vector
        , size_t N          // Number of elements
        , bool TF           // Transpose flag
        , AlignmentFlag AF  // Alignment flag
        , PaddingFlag PF    // Padding flag
        , typename Tag >    // Type tag
inline void swap( HybridVector<Type,N,TF,AF,PF,Tag>& a, HybridVector<Type,N,TF,AF,PF,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct MaxSize< HybridVector<T,N,TF,AF,PF,Tag>, 0UL >
   : public Ptrdiff_t<N>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct HasConstDataAccess< HybridVector<T,N,TF,AF,PF,Tag> >
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
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct HasMutableDataAccess< HybridVector<T,N,TF,AF,PF,Tag> >
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
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct IsAligned< HybridVector<T,N,TF,AF,PF,Tag> >
   : public BoolConstant< AF == aligned >
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
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct IsContiguous< HybridVector<T,N,TF,AF,PF,Tag> >
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
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct IsPadded< HybridVector<T,N,TF,AF,PF,Tag> >
   : public BoolConstant< PF == padded >
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
template< typename T, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag >
struct IsResizable< HybridVector<T,N,TF,AF,PF,Tag> >
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                    MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                            , min( size_t( MaxSize_v<T1,0UL> ), size_t( MaxSize_v<T2,0UL> ) )
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                    MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                            , min( size_t( MaxSize_v<T1,0UL> ), size_t( MaxSize_v<T2,0UL> ) )
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , SubTrait_t< TagType_t<T1>, TagType_t<T2> > >;
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
                     , EnableIf_t< IsVector_v<T1> &&
                                   IsScalar_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< MultTrait_t< ElementType_t<T1>, T2 >
                            , MaxSize_v<T1,0UL>
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsVector_v<T2> &&
                                   ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< MultTrait_t< T1, ElementType_t<T2> >
                            , MaxSize_v<T2,0UL>
                            , TransposeFlag_v<T2>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                     ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                   IsDenseVector_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                   ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                   ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                     MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                            , min( size_t( MaxSize_v<T1,0UL> ), size_t( MaxSize_v<T2,0UL> ) )
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsColumnVector_v<T2> &&
                                   ( Size_v<T1,0UL> == DefaultSize_v &&
                                     ( !IsSquare_v<T1> || Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                   ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                     ( IsSquare_v<T1> && MaxSize_v<T2,0UL> != DefaultMaxSize_v ) ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = HybridVector< AddTrait_t<MultType,MultType>
                            , ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ? MaxSize_v<T1,0UL> : MaxSize_v<T2,0UL> )
                            , false
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , AddTrait_t<MultTag,MultTag> >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsRowVector_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( Size_v<T2,1UL> == DefaultSize_v &&
                                     ( !IsSquare_v<T2> || Size_v<T1,0UL> == DefaultSize_v ) ) &&
                                   ( MaxSize_v<T2,1UL> != DefaultMaxSize_v ||
                                     ( IsSquare_v<T2> && MaxSize_v<T1,0UL> != DefaultMaxSize_v ) ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = HybridVector< AddTrait_t<MultType,MultType>
                            , ( MaxSize_v<T2,1UL> != DefaultMaxSize_v ? MaxSize_v<T2,1UL> : MaxSize_v<T1,0UL> )
                            , true
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
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
struct KronTraitEval2< T1, T2
                     , EnableIf_t< IsDenseVector_v<T1> &&
                                   IsDenseVector_v<T2> &&
                                   ( ( Size_v<T1,0UL> == DefaultSize_v ) ||
                                     ( Size_v<T2,0UL> == DefaultSize_v ) ) &&
                                   ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ) &&
                                   ( MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                            , MaxSize_v<T1,0UL> * MaxSize_v<T2,0UL>
                            , TransposeFlag_v<T2>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsScalar_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< DivTrait_t< ElementType_t<T1>, T2 >
                            , MaxSize_v<T1,0UL>
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , DivTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsDenseVector_v<T1> &&
                                  IsDenseVector_v<T2> &&
                                  ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                  ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                  ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                    MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< DivTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                            , min( size_t( MaxSize_v<T1,0UL> ), size_t( MaxSize_v<T2,0UL> ) )
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , DivTrait_t< TagType_t<T1>, TagType_t<T2> > >;
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
                         , EnableIf_t< IsVector_v<T> &&
                                       Size_v<T,0UL> == DefaultSize_v &&
                                       MaxSize_v<T,0UL> != DefaultMaxSize_v > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = HybridVector< EvaluateTrait_t<ElementType>
                            , MaxSize_v<T,0UL>
                            , TransposeFlag_v<T>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval2< T1, T2, OP
                          , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                          ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                        Size_v<T1,0UL> == DefaultSize_v &&
                                        Size_v<T2,0UL> == DefaultSize_v &&
                                        ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ||
                                          MaxSize_v<T2,0UL> != DefaultMaxSize_v ) > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   using Type = HybridVector< EvaluateTrait_t<ElementType>
                            , min( size_t( MaxSize_v<T1,0UL> ), size_t( MaxSize_v<T2,0UL> ) )
                            , TransposeFlag_v<T1>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REDUCETRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP, ReductionFlag RF >
struct PartialReduceTraitEval2< T, OP, RF
                              , EnableIf_t< IsMatrix_v<T> &&
                                            ( Size_v<T,0UL> == DefaultSize_v ||
                                              Size_v<T,1UL> == DefaultSize_v ) &&
                                            MaxSize_v<T,0UL> != DefaultMaxSize_v &&
                                            MaxSize_v<T,1UL> != DefaultMaxSize_v > >
{
   using ET = ElementType_t<T>;

   static constexpr bool TF = ( RF == columnwise );

   using Type = HybridVector< decltype( std::declval<OP>()( std::declval<ET>(), std::declval<ET>() ) )
                            , MaxSize_v< T, TF ? 1UL : 0UL >
                            , TF
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , MapTrait_t< TagType_t<T>, OP > >;
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
template< typename T, size_t R0 >
struct RepeatTraitEval2< T, R0, inf, inf
                       , EnableIf_t< IsDenseVector_v<T> &&
                                     ( R0 != inf ) &&
                                     ( Size_v<T,0UL> == DefaultSize_v ) &&
                                     ( MaxSize_v<T,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< ElementType_t<T>
                            , R0*MaxSize_v<T,0UL>
                            , TransposeFlag_v<T>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
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
                                    IsDenseVector_v<T2> &&
                                    ( Size_v<T1,0UL> == DefaultSize_v ) &&
                                    ( Size_v<T1,1UL> == DefaultSize_v ) &&
                                    ( Size_v<T2,0UL> == DefaultSize_v ) &&
                                    ( ( MaxSize_v<T1,0UL> != DefaultMaxSize_v ) ||
                                      ( MaxSize_v<T1,1UL> != DefaultMaxSize_v ) ||
                                      ( MaxSize_v<T2,0UL> != DefaultMaxSize_v ) ) > >
{
   static constexpr size_t N = ( MaxSize_v<T1,0UL> != DefaultMaxSize_v
                                 ?( MaxSize_v<T1,1UL> != DefaultMaxSize_v
                                    ?( MaxSize_v<T2,0UL> != DefaultMaxSize_v
                                       ? min( MaxSize_v<T1,0UL>, MaxSize_v<T1,1UL>, MaxSize_v<T2,0UL> )
                                       : min( MaxSize_v<T1,0UL>, MaxSize_v<T1,1UL> ) )
                                    :( MaxSize_v<T2,0UL> != DefaultMaxSize_v
                                       ? min( MaxSize_v<T1,0UL>, MaxSize_v<T2,0UL> )
                                       : MaxSize_v<T1,0UL> ) )
                                 :( MaxSize_v<T1,1UL> != DefaultMaxSize_v
                                    ?( MaxSize_v<T2,0UL> != DefaultMaxSize_v
                                       ? min( MaxSize_v<T1,1UL>, MaxSize_v<T2,0UL> )
                                       : MaxSize_v<T1,1UL> )
                                    : MaxSize_v<T2,0UL> ) );

   using Type = HybridVector< ElementType_t<T2>
                            , N
                            , TransposeFlag_v<T2>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
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
template< typename T1, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag, typename T2 >
struct HighType< HybridVector<T1,N,TF,AF,PF,Tag>, HybridVector<T2,N,TF,AF,PF,Tag> >
{
   using Type = StaticVector< typename HighType<T1,T2>::Type, N, TF, AF, PF, Tag >;
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
template< typename T1, size_t N, bool TF, AlignmentFlag AF, PaddingFlag PF, typename Tag, typename T2 >
struct LowType< HybridVector<T1,N,TF,AF,PF,Tag>, HybridVector<T2,N,TF,AF,PF,Tag> >
{
   using Type = StaticVector< typename LowType<T1,T2>::Type, N, TF, AF, PF, Tag >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT >
struct SubvectorTraitEval2< VT, inf, inf
                          , EnableIf_t< IsDenseVector_v<VT> &&
                                        ( Size_v<VT,0UL> != DefaultSize_v ||
                                          MaxSize_v<VT,0UL> != DefaultMaxSize_v ) > >
{
   using Type = HybridVector< RemoveConst_t< ElementType_t<VT> >
                            , max( Size_v<VT,0UL>, MaxSize_v<VT,0UL> )
                            , TransposeFlag_v<VT>
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<VT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I >
struct RowTraitEval2< MT, I
                    , EnableIf_t< IsDenseMatrix_v<MT> &&
                                  Size_v<MT,1UL> == DefaultSize_v &&
                                  MaxSize_v<MT,1UL> != DefaultMaxSize_v > >
{
   using Type = HybridVector< RemoveConst_t< ElementType_t<MT> >
                            , MaxSize_v<MT,1UL>
                            , true
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I >
struct ColumnTraitEval2< MT, I
                       , EnableIf_t< IsDenseMatrix_v<MT> &&
                                     Size_v<MT,0UL> == DefaultSize_v &&
                                     MaxSize_v<MT,0UL> != DefaultMaxSize_v > >
{
   using Type = HybridVector< RemoveConst_t< ElementType_t<MT> >
                            , MaxSize_v<MT,0UL>
                            , false
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, ptrdiff_t I >
struct BandTraitEval2< MT, I
                     , EnableIf_t< IsDenseMatrix_v<MT> &&
                                   ( Size_v<MT,0UL> == DefaultSize_v ||
                                     Size_v<MT,1UL> == DefaultSize_v ) &&
                                   MaxSize_v<MT,0UL> != DefaultMaxSize_v &&
                                   MaxSize_v<MT,1UL> != DefaultMaxSize_v > >
{
   static constexpr size_t M   = MaxSize_v<MT,0UL>;
   static constexpr size_t N   = MaxSize_v<MT,1UL>;
   static constexpr size_t Min = min( M - ( I >= 0L ? 0UL : -I ), N - ( I >= 0L ? I : 0UL ) );

   using Type = HybridVector< RemoveConst_t< ElementType_t<MT> >
                            , Min
                            , defaultTransposeFlag
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<MT> >;
};

template< typename MT >
struct BandTraitEval2< MT, inf
                     , EnableIf_t< IsDenseMatrix_v<MT> &&
                                   Size_v<MT,0UL> != DefaultSize_v &&
                                   Size_v<MT,1UL> != DefaultSize_v > >
{
   using Type = HybridVector< RemoveConst_t< ElementType_t<MT> >
                            , min( Size_v<MT,0UL>, Size_v<MT,1UL> )
                            , defaultTransposeFlag
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<MT> >;
};

template< typename MT >
struct BandTraitEval2< MT, inf
                     , EnableIf_t< IsDenseMatrix_v<MT> &&
                                   ( Size_v<MT,0UL> == DefaultSize_v ||
                                     Size_v<MT,1UL> == DefaultSize_v ) &&
                                   MaxSize_v<MT,0UL> != DefaultMaxSize_v &&
                                   MaxSize_v<MT,1UL> != DefaultMaxSize_v > >
{
   using Type = HybridVector< RemoveConst_t< ElementType_t<MT> >
                            , min( MaxSize_v<MT,0UL>, MaxSize_v<MT,1UL> )
                            , defaultTransposeFlag
                            , defaultAlignmentFlag
                            , defaultPaddingFlag
                            , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
