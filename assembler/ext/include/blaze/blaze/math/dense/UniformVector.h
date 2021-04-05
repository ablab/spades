//=================================================================================================
/*!
//  \file blaze/math/dense/UniformVector.h
//  \brief Header file for the implementation of a uniform vector
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

#ifndef _BLAZE_MATH_DENSE_UNIFORMVECTOR_H_
#define _BLAZE_MATH_DENSE_UNIFORMVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/dense/UniformIterator.h>
#include <blaze/math/dense/Forward.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/ReduceTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/Thresholds.h>
#include <blaze/system/TransposeFlag.h>
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
/*!\defgroup uniform_vector UniformVector
// \ingroup dense_vector
*/
/*!\brief Efficient implementation of a uniform vector.
// \ingroup uniform_vector
//
// The UniformVector class template is the representation of an arbitrary sized uniform vector
// with elements of arbitrary type. The type of the elements, the transpose flag, and the group
// tag of the vector can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class UniformVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. UniformVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// These uniform elements can be directly accessed with the subscript operator. The numbering
// of the vector elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right)\f]

// The use of UniformVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// vectors with fitting element types. The following example gives an impression of the use of
// UniformVector:

   \code
   using blaze::UniformVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;

   UniformVector<double> a( 2 );  // Default initialized 2D vector of size 2
   a = 1.0;                       // Assignment to all elements of the uniform vector

   UniformVector<double>   b( 2, 2.0  );  // Directly, uniformly initialized 2D vector
   CompressedVector<float> c( 2 );        // Empty sparse single precision vector
   UniformVector<double>   d;             // Default constructed uniform vector
   DynamicMatrix<double>   A;             // Default constructed row-major matrix

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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
class UniformVector
   : public Expression< DenseVector< UniformVector<Type,TF,Tag>, TF > >
{
 public:
   //**Type definitions****************************************************************************
   using This     = UniformVector<Type,TF,Tag>;          //!< Type of this UniformVector instance.
   using BaseType = Expression< DenseVector<This,TF> >;  //!< Base type of this UniformVector instance.

   //! Result type for expression template evaluations.
   using ResultType = This;

   //! Transpose type for expression template evaluations.
   using TransposeType = UniformVector<Type,!TF,Tag>;

   using ElementType   = Type;                      //!< Type of the vector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;  //!< SIMD type of the vector elements.
   using TagType       = Tag;                       //!< Tag type of this UniformVector instance.
   using ReturnType    = const Type&;               //!< Return type for expression template evaluations
   using CompositeType = const UniformVector&;      //!< Data type for composite expression templates.

   using Reference      = const Type&;  //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = const Type*;  //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = UniformIterator<const Type,aligned>;  //!< Iterator over non-constant elements.
   using ConstIterator = UniformIterator<const Type,aligned>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a UniformVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = UniformVector<NewType,TF,Tag>;  //!< The type of the other UniformVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a UniformVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = UniformVector<Type,TF,Tag>;  //!< The type of the other UniformVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operationss. In case the element type of the vector is a
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
            constexpr UniformVector() noexcept;
   explicit constexpr UniformVector( size_t n );
            constexpr UniformVector( size_t n, const Type& init );

   template< typename VT >
   inline UniformVector( const Vector<VT,TF>& v );

   UniformVector( const UniformVector& ) = default;
   UniformVector( UniformVector&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~UniformVector() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr ConstReference operator[]( size_t index ) const noexcept;
   inline    ConstReference at( size_t index ) const;
   constexpr ConstPointer   data  () const noexcept;
   constexpr ConstIterator  begin () const noexcept;
   constexpr ConstIterator  cbegin() const noexcept;
   constexpr ConstIterator  end   () const noexcept;
   constexpr ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   constexpr UniformVector& operator=( const Type& rhs ) &;

   UniformVector& operator=( const UniformVector& ) & = default;
   UniformVector& operator=( UniformVector&& ) & = default;

   template< typename VT > inline UniformVector& operator= ( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline UniformVector& operator+=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline UniformVector& operator-=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline UniformVector& operator*=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline UniformVector& operator/=( const DenseVector<VT,TF>& rhs ) &;

   template< typename ST >
   inline auto operator*=( ST rhs ) & -> EnableIf_t< IsScalar_v<ST>, UniformVector& >;

   template< typename ST >
   inline auto operator/=( ST rhs ) & -> EnableIf_t< IsScalar_v<ST>, UniformVector& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   constexpr size_t size() const noexcept;
   constexpr size_t spacing() const noexcept;
   constexpr size_t capacity() const noexcept;
   inline    size_t nonZeros() const;
   constexpr void   reset();
   constexpr void   clear();
   constexpr void   resize( size_t n, bool preserve=true );
   constexpr void   extend( size_t n, bool preserve=true );
   constexpr void   swap( UniformVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline UniformVector& scale( const Other& scalar );
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

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;   //!< The current size/dimension of the uniform vector.
   Type   value_;  //!< The value of all elements of the uniform vector.
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
/*!\brief The default constructor for UniformVector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr UniformVector<Type,TF,Tag>::UniformVector() noexcept
   : size_ ( 0UL )  // The current size/dimension of the uniform vector
   , value_()       // The value of all elements of the uniform vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a vector of size \a n.
//
// \param n The size of the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr UniformVector<Type,TF,Tag>::UniformVector( size_t n )
   : size_ ( n )  // The current size/dimension of the uniform vector
   , value_()     // The value of all elements of the uniform vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a homogeneous initialization of all \a n vector elements.
//
// \param n The size of the vector.
// \param init The initial value of the vector elements.
//
// All vector elements are initialized with the specified value.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr UniformVector<Type,TF,Tag>::UniformVector( size_t n, const Type& init )
   : size_ ( n    )  // The current size/dimension of the uniform vector
   , value_( init )  // The value of all elements of the uniform vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor for different uniform vectors.
//
// \param v Uniform vector to be copied.
// \exception std::invalid_argument Invalid setup of uniform vector.
//
// The vector is sized according to the given \f$ N \f$-dimensional uniform vector and
// initialized as a copy of this vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the foreign vector
inline UniformVector<Type,TF,Tag>::UniformVector( const Vector<VT,TF>& v )
   : size_ ( (*v).size() )  // The current size/dimension of the uniform vector
   , value_()               // The value of all elements of the uniform vector
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( !IsUniform_v<VT> && !isUniform( *v ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of uniform vector" );
   }

   if( size_ > 0UL ) {
      value_ = (*v)[0UL];
   }
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstReference
   UniformVector<Type,TF,Tag>::operator[]( size_t index ) const noexcept
{
   MAYBE_UNUSED( index );

   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );

   return value_;
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename UniformVector<Type,TF,Tag>::ConstReference
   UniformVector<Type,TF,Tag>::at( size_t index ) const
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
// This function returns a pointer to the internal storage of the uniform vector. Note that you
// can NOT assume that elements of the uniform vector lie adjacent to each other!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstPointer
   UniformVector<Type,TF,Tag>::data() const noexcept
{
   return &value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the uniform vector.
//
// \return Iterator to the first element of the uniform vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstIterator
   UniformVector<Type,TF,Tag>::begin() const noexcept
{
   return ConstIterator( &value_, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the uniform vector.
//
// \return Iterator to the first element of the uniform vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstIterator
   UniformVector<Type,TF,Tag>::cbegin() const noexcept
{
   return ConstIterator( &value_, 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the uniform vector.
//
// \return Iterator just past the last element of the uniform vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstIterator
   UniformVector<Type,TF,Tag>::end() const noexcept
{
   return ConstIterator( &value_, size_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the uniform vector.
//
// \return Iterator just past the last element of the uniform vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename UniformVector<Type,TF,Tag>::ConstIterator
   UniformVector<Type,TF,Tag>::cend() const noexcept
{
   return ConstIterator( &value_, size_ );
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator=( const Type& rhs ) &
{
   value_ = rhs;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be copied.
// \return Reference to the assigned vector.
//
// The vector is resized according to the given uniform vector and initialized as a copy of this
// vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side vector
inline UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator=( const Vector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( !IsUniform_v<VT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment of uniform vector" );
   }

   size_ = (*rhs).size();

   if( (*rhs).size() > 0UL ) {
      value_ = (*rhs)[0UL];
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
// \exception std::invalid_argument Invalid addition assignment to uniform vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side vector
inline UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator+=( const Vector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsUniform_v<VT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid addition assignment to uniform vector" );
   }

   if( (*rhs).size() > 0UL ) {
      value_ += (*rhs)[0UL];
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid subtraction assignment to uniform vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side vector
inline UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator-=( const Vector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsUniform_v<VT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subtraction assignment to uniform vector" );
   }

   if( (*rhs).size() > 0UL ) {
      value_ -= (*rhs)[0UL];
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the subtraction of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the vector.
// \return Reference to the vector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid multiplication assignment to uniform vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side vector
inline UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator*=( const Vector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsUniform_v<VT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid multiplication assignment to uniform vector" );
   }

   if( (*rhs).size() > 0UL ) {
      value_ *= (*rhs)[0UL];
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
// \exception std::invalid_argument Invalid multiplication assignment to uniform vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side vector
inline UniformVector<Type,TF,Tag>&
   UniformVector<Type,TF,Tag>::operator/=( const DenseVector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsUniform_v<VT> && !isUniform( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid division assignment to uniform vector" );
   }

   if( (*rhs).size() > 0UL ) {
      value_ /= (*rhs)[0UL];
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication between a vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
//
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename ST >   // Data type of the right-hand side scalar
inline auto UniformVector<Type,TF,Tag>::operator*=( ST scalar ) &
   -> EnableIf_t< IsScalar_v<ST>, UniformVector& >
{
   if( size() > 0UL ) {
      value_ *= scalar;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division between a vector and a scalar value
//        (\f$ \vec{a}/=s \f$).
//
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename ST >   // Data type of the right-hand side scalar
inline auto UniformVector<Type,TF,Tag>::operator/=( ST scalar ) &
   -> EnableIf_t< IsScalar_v<ST>, UniformVector& >
{
   if( size() > 0UL ) {
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
/*!\brief Returns the current size/dimension of the vector.
//
// \return The size of the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t UniformVector<Type,TF,Tag>::size() const noexcept
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t UniformVector<Type,TF,Tag>::spacing() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the vector.
//
// \return The maximum capacity of the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t UniformVector<Type,TF,Tag>::capacity() const noexcept
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline size_t UniformVector<Type,TF,Tag>::nonZeros() const
{
   if( size_ == 0UL || isDefault<strict>( value_ ) )
      return 0UL;
   else
      return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void UniformVector<Type,TF,Tag>::reset()
{
   using blaze::clear;

   clear( value_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the vector.
//
// \return void
//
// After the clear() function, the size of the vector is 0.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void UniformVector<Type,TF,Tag>::clear()
{
   size_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the vector.
//
// \param n The new size of the vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function resizes the vector using the given size to \a n. Note that this function may
// invalidate all existing views (subvectors, ...) on the vector if it is used to shrink the
// vector. Additionally, the resize operation potentially changes all vector elements. In order
// to preserve the old vector values, the \a preserve flag can be set to \a true.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void UniformVector<Type,TF,Tag>::resize( size_t n, bool preserve )
{
   MAYBE_UNUSED( preserve );

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
// This function increases the vector size by \a n elements. Note that this function potentially
// changes all vector elements. In order to preserve the old vector values, the \a preserve flag
// can be set to \a true.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void UniformVector<Type,TF,Tag>::extend( size_t n, bool preserve )
{
   resize( size_+n, preserve );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
//
// \param v The vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void UniformVector<Type,TF,Tag>::swap( UniformVector& v ) noexcept
{
   using std::swap;

   swap( size_, v.size_ );
   swap( value_, v.value_ );
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
   blaze::UniformVector<int> a;
   // ... Initialization
   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline UniformVector<Type,TF,Tag>& UniformVector<Type,TF,Tag>::scale( const Other& scalar )
{
   if( size_ > 0UL ) {
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
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool UniformVector<Type,TF,Tag>::canAlias( const Other* alias ) const noexcept
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
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool UniformVector<Type,TF,Tag>::isAliased( const Other* alias ) const noexcept
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline bool UniformVector<Type,TF,Tag>::isAligned() const noexcept
{
   return true;
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline bool UniformVector<Type,TF,Tag>::canSMPAssign() const noexcept
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformVector<Type,TF,Tag>::SIMDType
   UniformVector<Type,TF,Tag>::load( size_t index ) const noexcept
{
   return loada( index );
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformVector<Type,TF,Tag>::SIMDType
   UniformVector<Type,TF,Tag>::loada( size_t index ) const noexcept
{
   MAYBE_UNUSED( index );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL, "Invalid vector access index" );

   return set( value_ );
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
BLAZE_ALWAYS_INLINE typename UniformVector<Type,TF,Tag>::SIMDType
   UniformVector<Type,TF,Tag>::loadu( size_t index ) const noexcept
{
   MAYBE_UNUSED( index );

   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( Type );

   BLAZE_INTERNAL_ASSERT( index < size_, "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size_, "Invalid vector access index" );

   return set( value_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  UNIFORMVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name UniformVector operators */
//@{
template< typename Type, bool TF, typename Tag >
constexpr void reset( UniformVector<Type,TF,Tag>& v );

template< typename Type, bool TF, typename Tag >
constexpr void clear( UniformVector<Type,TF,Tag>& v );

template< RelaxationFlag RF, typename Type, bool TF, typename Tag >
constexpr bool isDefault( const UniformVector<Type,TF,Tag>& v );

template< typename Type, bool TF, typename Tag >
constexpr bool isIntact( const UniformVector<Type,TF,Tag>& v ) noexcept;

template< typename Type, bool TF, typename Tag >
constexpr void swap( UniformVector<Type,TF,Tag>& a, UniformVector<Type,TF,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given uniform vector.
// \ingroup uniform_vector
//
// \param v The uniform vector to be resetted.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void reset( UniformVector<Type,TF,Tag>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given uniform vector.
// \ingroup uniform_vector
//
// \param v The uniform vector to be cleared.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void clear( UniformVector<Type,TF,Tag>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given uniform vector is in default state.
// \ingroup uniform_vector
//
// \param v The uniform vector to be tested for its default state.
// \return \a true in case the given vector's size is zero, \a false otherwise.
//
// This function checks whether the uniform vector is in default (constructed) state, i.e. if
// it's size is 0. In case it is in default state, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isDefault() function:

   \code
   blaze::UniformVector<int> a;
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
        , bool TF            // Transpose flag
        , typename Tag >     // Type tag
constexpr bool isDefault( const UniformVector<Type,TF,Tag>& v )
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given uniform vector are intact.
// \ingroup uniform_vector
//
// \param v The uniform vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the uniform vector are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it
// will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::UniformVector<int> a;
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr bool isIntact( const UniformVector<Type,TF,Tag>& v ) noexcept
{
   MAYBE_UNUSED( v );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
// \ingroup uniform_vector
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void swap( UniformVector<Type,TF,Tag>& a, UniformVector<Type,TF,Tag>& b ) noexcept
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
/*!\brief Creating a uniform vector.
// \ingroup uniform_vector
//
// \param n The size of the vector.
// \param init The initial value of the vector elements.
// \return A uniform vector of the given size.
//
// This function creates a uniform vector of the given size. By default, the resulting uniform
// vector is a column vector, but it is possible to specify the transpose flag explicitly:

   \code
   using blaze::uniform;
   using blaze::columnVector;
   using blaze::rowVector;

   // Creates the uniform column vector ( 1, 1, 1, 1, 1 )
   auto u1 = uniform( 5UL, 1 );

   // Creates the uniform column vector ( 1.2, 1.2, 1.2 )
   auto u2 = uniform<columnVector>( 3UL, 1.2 );

   // Creates the uniform row vector ( 5U, 5U, 5U, 5U )
   auto u3 = uniform<rowVector>( 4UL, 5U );
   \endcode
*/
template< bool TF = defaultTransposeFlag, typename T >
constexpr decltype(auto) uniform( size_t n, T&& init )
{
   return UniformVector< RemoveCVRef_t<T>, TF >( n, std::forward<T>( init ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename Type, bool TF, typename Tag >
struct IsUniform< UniformVector<Type,TF,Tag> >
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
template< typename Type, bool TF, typename Tag >
struct IsAligned< UniformVector<Type,TF,Tag> >
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
template< typename Type, bool TF, typename Tag >
struct IsResizable< UniformVector<Type,TF,Tag> >
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformVector< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , TransposeFlag_v<T1>
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformVector< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , TransposeFlag_v<T1>
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
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsVector_v<T1> &&
                                   IsScalar_v<T2> &&
                                   IsUniform_v<T1> && !IsZero_v<T1> > >
{
   using Type = UniformVector< MultTrait_t< ElementType_t<T1>, T2 >
                             , TransposeFlag_v<T1>
                             , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsVector_v<T2> &&
                                   IsUniform_v<T2> && !IsZero_v<T2> > >
{
   using Type = UniformVector< MultTrait_t< T1, ElementType_t<T2> >
                             , TransposeFlag_v<T2>
                             , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                     ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                   ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , TransposeFlag_v<T1>
                             , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsColumnVector_v<T2> &&
                                   IsUniform_v<T1> &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = UniformVector< AddTrait_t<MultType,MultType>
                             , false
                             , AddTrait_t<MultTag,MultTag> >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsRowVector_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   IsUniform_v<T2> &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = UniformVector< AddTrait_t<MultType,MultType>
                             , true
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
                     , EnableIf_t< IsVector_v<T1> &&
                                   IsVector_v<T2> &&
                                   ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , TransposeFlag_v<T2>
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
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsScalar_v<T2> &&
                                  IsUniform_v<T1> && !IsZero_v<T1> > >
{
   using Type = UniformVector< DivTrait_t< ElementType_t<T1>, T2 >
                             , TransposeFlag_v<T1>
                             , DivTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct DivTraitEval1< T1, T2
                    , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                    ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                  ( IsUniform_v<T1> && IsUniform_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = UniformVector< DivTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                             , TransposeFlag_v<T1>
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
struct UnaryMapTraitEval1< T, OP
                         , EnableIf_t< IsVector_v<T> &&
                                       YieldsUniform_v<OP,T> &&
                                       !YieldsZero_v<OP,T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = UniformVector< EvaluateTrait_t<ElementType>
                             , TransposeFlag_v<T>
                             , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                          ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                        YieldsUniform_v<OP,T1,T2> &&
                                        !YieldsZero_v<OP,T1,T2> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   using Type = UniformVector< EvaluateTrait_t<ElementType>
                             , TransposeFlag_v<T1>
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
struct PartialReduceTraitEval1< T, OP, RF
                              , EnableIf_t< IsMatrix_v<T> &&
                                            ( IsUniform_v<T> || IsIdentity_v<T> ) > >
{
   using ET = ElementType_t<T>;

   using Type = UniformVector< decltype( std::declval<OP>()( std::declval<ET>(), std::declval<ET>() ) )
                             , ( RF == columnwise )
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
struct RepeatTraitEval1< T, R0, inf, inf
                       , EnableIf_t< IsVector_v<T> &&
                                     IsUniform_v<T> && !IsZero_v<T> > >
{
   using Type = UniformVector< ElementType_t<T>
                             , TransposeFlag_v<T>
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
template< typename T1, bool TF, typename Tag, typename T2 >
struct HighType< UniformVector<T1,TF,Tag>, UniformVector<T2,TF,Tag> >
{
   using Type = UniformVector< typename HighType<T1,T2>::Type, TF, Tag >;
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
template< typename T1, bool TF, typename Tag, typename T2 >
struct LowType< UniformVector<T1,TF,Tag>, UniformVector<T2,TF,Tag> >
{
   using Type = UniformVector< typename LowType<T1,T2>::Type, TF, Tag >;
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
template< typename VT, size_t I, size_t N >
struct SubvectorTraitEval1< VT, I, N
                          , EnableIf_t< IsUniform_v<VT> && !IsZero_v<VT> > >
{
   using Type = UniformVector< RemoveConst_t< ElementType_t<VT> >
                             , TransposeFlag_v<VT>
                             , TagType_t<VT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENTSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, size_t N >
struct ElementsTraitEval1< VT, N
                         , EnableIf_t< IsUniform_v<VT> && !IsZero_v<VT> > >
{
   using Type = UniformVector< RemoveConst_t< ElementType_t<VT> >
                             , TransposeFlag_v<VT>
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
struct RowTraitEval1< MT, I
                    , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformVector< RemoveConst_t< ElementType_t<MT> >
                             , true
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
struct ColumnTraitEval1< MT, I
                       , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformVector< RemoveConst_t< ElementType_t<MT> >
                             , false
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
struct BandTraitEval1< MT, I
                     , EnableIf_t< IsUniform_v<MT> && !IsZero_v<MT> > >
{
   using Type = UniformVector< RemoveConst_t< ElementType_t<MT> >
                             , defaultTransposeFlag
                             , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
