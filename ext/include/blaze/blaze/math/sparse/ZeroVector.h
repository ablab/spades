//=================================================================================================
/*!
//  \file blaze/math/sparse/ZeroVector.h
//  \brief Implementation of a zero vector
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

#ifndef _BLAZE_MATH_SPARSE_ZEROVECTOR_H_
#define _BLAZE_MATH_SPARSE_ZEROVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup zero_vector ZeroVector
// \ingroup sparse_vector
*/
/*!\brief Efficient implementation of an arbitrary sized zero vector.
// \ingroup zero_vector
//
// The ZeroVector class template is the representation of an immutable, arbitrary sized zero
// vector with N elements of arbitrary type. The type of the elements, the transpose flag, and
// the group tag of the vector can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class ZeroVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. ZeroVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//          vector (\a blaze::columnVector). The default value is \a blaze::defaultTransposeFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \a blaze::Group0.
//          See \ref grouping_tagging for details.
//
// It is not possible to insert, erase or modify the elements of a zero vector. It is only
// possible to read from the elements:

   \code
   using blaze::ZeroVector;
   using blaze::columnVector;

   // Creating a 4D zero column vector
   ZeroVector<double,columnVector> a( 4 );

   // The subscript operator provides access to all possible elements of the zero vector,
   // including the zero elements.
   a[1] = 2.0;       // Compilation error: It is not possible to write to a zero vector
   double d = a[2];  // Access to the element at index 2

   // In order to traverse all non-zero elements currently stored in the vector, the begin()
   // and end() functions can be used.
   for( ZeroVector<double,columnVector>::Iterator i=a.begin(); i!=a.end(); ++i ) {
      ... = i->value();  // Access to the value of the non-zero element
      ... = i->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of ZeroVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, ...) can be performed on all possible combinations of dense and sparse vectors
// with fitting element types. The following example gives an impression of the use of ZeroVector:

   \code
   using blaze::ZeroVector;
   using blaze::CompressedVector;
   using blaze::DynamicVector;

   ZeroVector<double> z( 3 );  // 3-dimensional zero vector

   DynamicVector<double,columnMajor> a( 3 );  // 3-dimensional dynamic dense vector
   CompressedVector<double,rowMajor> b( 3 );  // 3-dimensional double precision sparse vector
   CompressedVector<float,rowMajor>  c( 3 );  // 3-dimensional single precision sparse vector
   // ... Initialization of a, b, and c

   DynamicVector<double>    d( z );  // Creation of a new dense vector as a copy of z
   CompressedVector<double> e;       // Creation of a default sparse vector

   d = z + a;    // Addition of a zero vector and a dense vector
   d = b - z;    // Subtraction of a sparse vector and a zero vector
   e = z * c;    // Vector multiplication between two vectors of different element types

   d = 2.0 * z;  // Scaling of a zero vector
   e = z * 2.0;  // Scaling of a zero vector
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
class ZeroVector
   : public Expression< SparseVector< ZeroVector<Type,TF,Tag>, TF > >
{
 private:
   //**Type definitions****************************************************************************
   using Element = ValueIndexPair<Type>;  //!< Value-index-pair for the ZeroVector class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This     = ZeroVector<Type,TF,Tag>;              //!< Type of this ZeroVector instance.
   using BaseType = Expression< SparseVector<This,TF> >;  //!< Base type of this ZeroVector instance.

   //! Result type for expression template evaluations.
   using ResultType = This;

   //! Transpose type for expression template evaluations.
   using TransposeType = ZeroVector<Type,!TF,Tag>;

   using ElementType = Type;         //!< Type of the zero vector elements.
   using TagType     = Tag;          //!< Tag type of this ZeroVector instance.
   using ReturnType  = const Type&;  //!< Return type for expression template evaluations.

   using CompositeType  = const This&;     //!< Data type for composite expression templates.
   using Reference      = const Type&;     //!< Reference to a zero vector element.
   using ConstReference = const Type&;     //!< Reference to a constant zero vector element.
   using Iterator       = Element*;        //!< Iterator over non-constant elements.
   using ConstIterator  = const Element*;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a ZeroVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = ZeroVector<NewType,TF,Tag>;  //!< The type of the other ZeroVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a ZeroVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = ZeroVector<Type,TF,Tag>;  //!< The type of the other ZeroVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
            constexpr ZeroVector() noexcept;
   explicit constexpr ZeroVector( size_t size ) noexcept;

   template< typename VT > inline ZeroVector( const Vector<VT,TF>& v );

   ZeroVector( const ZeroVector& ) = default;
   ZeroVector( ZeroVector&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ZeroVector() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr ConstReference operator[]( size_t index ) const noexcept;
   inline    ConstReference at( size_t index ) const;
   constexpr ConstIterator  begin () const noexcept;
   constexpr ConstIterator  cbegin() const noexcept;
   constexpr ConstIterator  end   () const noexcept;
   constexpr ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   template< typename VT >
   inline ZeroVector& operator=( const Vector<VT,TF>& rhs ) &;

   ZeroVector& operator=( const ZeroVector& ) & = default;
   ZeroVector& operator=( ZeroVector&& ) & = default;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   constexpr size_t size() const noexcept;
   constexpr size_t capacity() const noexcept;
   constexpr size_t nonZeros() const noexcept;
   constexpr void   clear() noexcept;
   constexpr void   resize( size_t n ) noexcept;
   constexpr void   swap( ZeroVector& v ) noexcept;
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline ConstIterator find      ( size_t index ) const;
   inline ConstIterator lowerBound( size_t index ) const;
   inline ConstIterator upperBound( size_t index ) const;
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
   size_t size_;  //!< The current size/dimension of the zero vector.

   static const Type zero_;  //!< The zero element.
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

template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
const Type ZeroVector<Type,TF,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for ZeroVector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr ZeroVector<Type,TF,Tag>::ZeroVector() noexcept
   : size_( 0UL )  // The current size/dimension of the zero vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a zero vector of size \a n.
//
// \param n The size of the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr ZeroVector<Type,TF,Tag>::ZeroVector( size_t n ) noexcept
   : size_( n )  // The current size/dimension of the zero vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor for different zero vectors.
//
// \param v Zero vector to be copied.
// \exception std::invalid_argument Invalid setup of zero vector.
//
// The vector is sized according to the given N-dimensional zero vector and initialized as a
// copy of this vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the foreign zero vector
inline ZeroVector<Type,TF,Tag>::ZeroVector( const Vector<VT,TF>& v )
   : size_( (*v).size() )  // The current size/dimension of the zero vector
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( !IsZero_v<VT> && !isZero( *v ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid setup of zero vector" );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the zero vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename ZeroVector<Type,TF,Tag>::ConstReference
   ZeroVector<Type,TF,Tag>::operator[]( size_t index ) const noexcept
{
   MAYBE_UNUSED( index );

   BLAZE_USER_ASSERT( index < size(), "Invalid zero vector access index" );

   return zero_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the zero vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid zero vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename ZeroVector<Type,TF,Tag>::ConstReference
   ZeroVector<Type,TF,Tag>::at( size_t index ) const
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid zero vector access index" );
   }

   return zero_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the zero vector.
//
// \return Iterator to the first non-zero element of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::begin() const noexcept
{
   return nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the zero vector.
//
// \return Iterator to the first non-zero element of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::cbegin() const noexcept
{
   return nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the zero vector.
//
// \return Iterator just past the last non-zero element of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::end() const noexcept
{
   return nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the zero vector.
//
// \return Iterator just past the last non-zero element of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::cend() const noexcept
{
   return nullptr;
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different zero vectors.
//
// \param rhs Zero vector to be copied.
// \return Reference to the assigned vector.
// \exception std::invalid_argument Invalid assignment to zero vector.
//
// The vector is resized according to the given zero vector and initialized as a copy of this
// vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >   // Type of the right-hand side dense vector
inline ZeroVector<Type,TF,Tag>&
   ZeroVector<Type,TF,Tag>::operator=( const Vector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( !IsZero_v<VT> && !isZero( *rhs ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment of zero vector" );
   }

   size_ = (*rhs).size();

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the zero vector.
//
// \return The size of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t ZeroVector<Type,TF,Tag>::size() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the zero vector.
//
// \return The capacity of the zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t ZeroVector<Type,TF,Tag>::capacity() const noexcept
{
   return 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the zero vector.
//
// \return The number of non-zero elements in the zero vector.
//
// Note that the number of non-zero elements is always smaller than the current size of the
// zero vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr size_t ZeroVector<Type,TF,Tag>::nonZeros() const noexcept
{
   return 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the zero vector.
//
// \return void
//
// After the clear() function, the size of the zero vector is 0.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void ZeroVector<Type,TF,Tag>::clear() noexcept
{
   size_ = 0UL;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the zero vector.
//
// \param n The new size of the zero vector.
// \return void
//
// This function resizes the zero vector using the given size to \a n. Note that this function
// may invalidate all existing views (subvectors, ...) on the vector if it is used to shrink the
// vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void ZeroVector<Type,TF,Tag>::resize( size_t n ) noexcept
{
   size_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two zero vectors.
//
// \param v The zero vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void ZeroVector<Type,TF,Tag>::swap( ZeroVector& v ) noexcept
{
   const size_t tmp( size_ );
   size_ = v.size_;
   v.size_ = tmp;
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the zero vector (the end() iterator) is returned.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::find( size_t index ) const
{
   MAYBE_UNUSED( index );

   BLAZE_USER_ASSERT( index < size(), "Invalid zero vector access index" );

   return nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::lowerBound( size_t index ) const
{
   MAYBE_UNUSED( index );

   BLAZE_USER_ASSERT( index < size(), "Invalid zero vector access index" );

   return nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename ZeroVector<Type,TF,Tag>::ConstIterator
   ZeroVector<Type,TF,Tag>::upperBound( size_t index ) const
{
   MAYBE_UNUSED( index );

   BLAZE_USER_ASSERT( index < size(), "Invalid zero vector access index" );

   return nullptr;
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
inline bool ZeroVector<Type,TF,Tag>::canAlias( const Other* alias ) const noexcept
{
   MAYBE_UNUSED( alias );

   return false;
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
inline bool ZeroVector<Type,TF,Tag>::isAliased( const Other* alias ) const noexcept
{
   MAYBE_UNUSED( alias );

   return false;
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
inline bool ZeroVector<Type,TF,Tag>::canSMPAssign() const noexcept
{
   return false;
}
//*************************************************************************************************








//=================================================================================================
//
//  ZEROVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name ZeroVector operators */
//@{
template< typename Type, bool TF, typename Tag >
constexpr void reset( ZeroVector<Type,TF,Tag>& m ) noexcept;

template< typename Type, bool TF, typename Tag >
constexpr void clear( ZeroVector<Type,TF,Tag>& m ) noexcept;

template< RelaxationFlag RF, typename Type, bool TF, typename Tag >
constexpr bool isDefault( const ZeroVector<Type,TF,Tag>& m ) noexcept;

template< typename Type, bool TF, typename Tag >
constexpr bool isIntact( const ZeroVector<Type,TF,Tag>& m ) noexcept;

template< typename Type, bool TF, typename Tag >
constexpr void swap( ZeroVector<Type,TF,Tag>& a, ZeroVector<Type,TF,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given zero vector.
// \ingroup zero_vector
//
// \param v The zero vector to be resetted.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void reset( ZeroVector<Type,TF,Tag>& v ) noexcept
{
   MAYBE_UNUSED( v );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given zero vector.
// \ingroup zero_vector
//
// \param v The zero vector to be cleared.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void clear( ZeroVector<Type,TF,Tag>& v ) noexcept
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given zero vector is in default state.
// \ingroup zero_vector
//
// \param v The zero vector to be tested for its default state.
// \return \a true in case the given vector's size is zero, \a false otherwise.
//
// This function checks whether the zero vector is in default (constructed) state, i.e. if
// it's size is 0. In case it is in default state, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isDefault() function:

   \code
   blaze::ZeroVector<double> z;
   // ... Resizing and initialization
   if( isDefault( z ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( z ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the vector
        , bool TF            // Transpose flag
        , typename Tag >     // Type tag
constexpr bool isDefault( const ZeroVector<Type,TF,Tag>& v ) noexcept
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given zero vector are intact.
// \ingroup zero_vector
//
// \param v The zero vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the zero vector are intact, i.e. if its state
// is valid. In case the invariants are intact, the function returns \a true, else it will return
// \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::ZeroVector<double> z;
   // ... Resizing and initialization
   if( isIntact( z ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr bool isIntact( const ZeroVector<Type,TF,Tag>& v ) noexcept
{
   MAYBE_UNUSED( v );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two zero vectors.
// \ingroup zero_vector
//
// \param a The first zero vector to be swapped.
// \param b The second zero vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
constexpr void swap( ZeroVector<Type,TF,Tag>& a, ZeroVector<Type,TF,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the given zero vector.
// \ingroup zero_vector
//
// \param v The given zero vector.
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void erase( ZeroVector<Type,TF,Tag>& v, size_t index )
{
   MAYBE_UNUSED( v, index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the given zero vector.
// \ingroup zero_vector
//
// \param v The given zero vector.
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the zero vector.
*/
template< typename Type        // Data type of the vector
        , bool TF              // Transpose flag
        , typename Tag         // Type tag
        , typename Iterator >  // Type of the vector iterator
inline Iterator erase( ZeroVector<Type,TF,Tag>& v, Iterator pos )
{
   MAYBE_UNUSED( v, pos );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the given zero vector.
// \ingroup zero_vector
//
// \param v The given zero vector.
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the zero vector.
*/
template< typename Type        // Data type of the vector
        , bool TF              // Transpose flag
        , typename Tag         // Type tag
        , typename Iterator >  // Type of the vector iterator
inline Iterator erase( ZeroVector<Type,TF,Tag>& m, Iterator first, Iterator last )
{
   MAYBE_UNUSED( m, first, last );

   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the given zero vector.
// \ingroup zero_vector
//
// \param v The given zero vector.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the given zero vector. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::ZeroVector<double,blaze::columnVector> z;
   // ... Resizing and initialization

   erase( z, []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the vector
        , bool TF          // Transpose flag
        , typename Tag     // Type tag
        , typename Pred >  // Type of the unary predicate
inline void erase( ZeroVector<Type,TF,Tag>& v, Pred predicate )
{
   MAYBE_UNUSED( v, predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the given zero vector.
// \ingroup zero_vector
//
// \param v The given zero vector.
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void
//
// This function erases specific elements from a range of elements of the zero vector. The
// elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements and to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::ZeroVector<double,blaze::columnVector> z;
   // ... Resizing and initialization

   erase( z, z.begin(), z.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type      // Data type of the vector
        , bool TF            // Transpose flag
        , typename Tag       // Type tag
        , typename Iterator  // Type of the vector iterator
        , typename Pred >    // Type of the unary predicate
inline void erase( ZeroVector<Type,TF,Tag>& m, Iterator first, Iterator last, Pred predicate )
{
   MAYBE_UNUSED( m, first, last, predicate );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a zero vector.
// \ingroup zero_vector
//
// \param n The size of the zero vector.
// \return A zero vector of the given size.
//
// This function creates a zero vector of the given element type and size. By default, the
// resulting zero vector is a column vector, but it is possible to specify the transpose flag
// explicitly:

   \code
   using blaze::zero;
   using blaze::columnVector;
   using blaze::rowVector;

   // Creates the zero column vector ( 0, 0, 0, 0, 0 )
   auto z1 = zero<int>( 5UL );

   // Creates the zero column vector ( 0.0, 0.0, 0.0 )
   auto z2 = zero<double,columnVector>( 3UL );

   // Creates the zero row vector ( 0U, 0U, 0U, 0U )
   auto z3 = zero<unsigned int,rowVector>( 4UL );
   \endcode
*/
template< typename T, bool TF = defaultTransposeFlag >
constexpr decltype(auto) zero( size_t n ) noexcept
{
   return ZeroVector<T,TF>( n );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Declares the given vector expression \a v as zero vector.
// \ingroup zero_vector
//
// \param v The input vector.
// \return The redeclared vector.
//
// The \a declzero function declares the given dense or sparse vector expression \a m as zero
// vector. The following example demonstrates the use of the \a declzero function:

   \code
   blaze::ZeroVector<double> a, b;
   // ... Resizing and initialization
   b = declzero( a );
   \endcode
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline ZeroVector<ElementType_t<VT>,TF>
   declzero( const Vector<VT,TF>& v )
{
   BLAZE_FUNCTION_TRACE;

   return ZeroVector<ElementType_t<VT>,TF>( (*v).size() );
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
struct IsUniform< ZeroVector<Type,TF,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISZERO SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename Type, bool TF, typename Tag >
struct IsZero< ZeroVector<Type,TF,Tag> >
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
struct IsResizable< ZeroVector<Type,TF,Tag> >
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
                                  !IsZero_v<T1> && IsZero_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T1>, AddTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct AddTraitEval1< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  IsZero_v<T1> && !IsZero_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T2>, AddTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct AddTraitEval1< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  IsZero_v<T1> && IsZero_v<T2> > >
{
   using Type = ZeroVector< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
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
                                  !IsZero_v<T1> && IsZero_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T1>, SubTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct SubTraitEval1< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  IsZero_v<T1> && !IsZero_v<T2> > >
{
   using Type = Rebind_t< ResultType_t<T2>, SubTrait_t< ElementType_t<T1>, ElementType_t<T2> > >;
};

template< typename T1, typename T2 >
struct SubTraitEval1< T1, T2
                    , EnableIf_t< IsVector_v<T1> &&
                                  IsVector_v<T2> &&
                                  IsZero_v<T1> && IsZero_v<T2> > >
{
   using Type = ZeroVector< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
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
                                   IsZero_v<T1> > >
{
   using Type = ZeroVector< MultTrait_t< ElementType_t<T1>, T2 >
                          , TransposeFlag_v<T1>
                          , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsVector_v<T2> &&
                                   IsZero_v<T2> > >
{
   using Type = ZeroVector< MultTrait_t< T1, ElementType_t<T2> >
                          , TransposeFlag_v<T2>
                          , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                     ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                   ( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = ZeroVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                          , TransposeFlag_v<T1>
                          , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsColumnVector_v<T2> &&
                                   ( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = ZeroVector< AddTrait_t<MultType,MultType>
                          , false
                          , AddTrait_t<MultTag,MultTag> >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsRowVector_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = ZeroVector< AddTrait_t<MultType,MultType>
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
                                   ( IsZero_v<T1> ||
                                     IsZero_v<T2> ) > >
{
   using Type = ZeroVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                          , ( IsZero_v<T2> ? TransposeFlag_v<T2> : TransposeFlag_v<T1> )
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
                                  IsZero_v<T1> > >
{
   using Type = ZeroVector< DivTrait_t< ElementType_t<T1>, T2 >
                          , TransposeFlag_v<T1>
                          , DivTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct DivTraitEval1< T1, T2
                    , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                    ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                  ( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = ZeroVector< DivTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                          , TransposeFlag_v<T1>
                          , DivTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CROSSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct CrossTraitEval1< T1, T2
                      , EnableIf_t< IsVector_v<T1> &&
                                    IsVector_v<T2> &&
                                    ( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = ZeroVector< SubTrait_t<MultType,MultType>
                          , TransposeFlag_v<T1>
                          , SubTrait_t<MultTag,MultTag> >;
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
                                       YieldsZero_v<OP,T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = ZeroVector< EvaluateTrait_t<ElementType>
                          , TransposeFlag_v<T>
                          , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< IsVector_v<T1> &&
                                        IsVector_v<T2> &&
                                        YieldsZero_v<OP,T1,T2> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T1> >()
                                                   , std::declval< ElementType_t<T2> >() ) );

   using Type = ZeroVector< EvaluateTrait_t<ElementType>
                          , TransposeFlag_v<T1>
                          , MapTrait_t< TagType_t<T1>, TagType_t<T2>, OP > >;
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
                                     IsZero_v<T> > >
{
   using Type = ZeroVector< ElementType_t<T>
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
struct HighType< ZeroVector<T1,TF,Tag>, ZeroVector<T2,TF,Tag> >
{
   using Type = ZeroVector< typename HighType<T1,T2>::Type, TF, Tag >;
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
struct LowType< ZeroVector<T1,TF,Tag>, ZeroVector<T2,TF,Tag> >
{
   using Type = ZeroVector< typename LowType<T1,T2>::Type, TF, Tag >;
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
                          , EnableIf_t< IsZero_v<VT> > >
{
   using Type = ZeroVector< RemoveConst_t< ElementType_t<VT> >
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
                         , EnableIf_t< IsZero_v<VT> > >
{
   using Type = ZeroVector< RemoveConst_t< ElementType_t<VT> >
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
                    , EnableIf_t< IsZero_v<MT> > >
{
   using Type = ZeroVector< RemoveConst_t< ElementType_t<MT> >
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
                       , EnableIf_t< IsZero_v<MT> > >
{
   using Type = ZeroVector< RemoveConst_t< ElementType_t<MT> >
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
                     , EnableIf_t< IsZero_v<MT> > >
{
   using Type = ZeroVector< RemoveConst_t< ElementType_t<MT> >
                          , defaultTransposeFlag
                          , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
