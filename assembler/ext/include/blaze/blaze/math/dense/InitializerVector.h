//=================================================================================================
/*!
//  \file blaze/math/dense/InitializerVector.h
//  \brief Header file for the implementation of a vector representation of an initializer list
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

#ifndef _BLAZE_MATH_DENSE_INITIALIZERVECTOR_H_
#define _BLAZE_MATH_DENSE_INITIALIZERVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/dense/InitializerIterator.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/IsDefault.h>
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
/*!\defgroup initializer_vector InitializerVector
// \ingroup dense_vector
*/
/*!\brief Dense vector representation of an initializer list.
// \ingroup initializer_vector
//
// The InitializerVector class template is a dense vector representation of an (extended)
// initializer list of arbitrary type. The type of the elements, the transpose flag, and the
// group tag of the vector can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class InitializerVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. InitializerVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\c blaze::rowVector) or a column
//          vector (\c blaze::columnVector). The default value is \c blaze::defaultTransposeFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \c blaze::Group0.
//          See \ref grouping_tagging for details.
//
// On construction, an InitializerVector is immediately bound to an initializer list:

   \code
   const auto list = { 2, 6, -1, 3, 5 };

   blaze::InitializerVector<int> a( list );  // Representation of the initializer list as dense column vector
   \endcode

// It is possible to only represent an extended initializer list by providing an additional
// size argument:

   \code
   const auto list = { 2, 6, -1, 3, 5 };

   blaze::InitializerVector<int> b( list, 5UL );  // Representation of the original initializer list
   blaze::InitializerVector<int> c( list, 8UL );  // Representing the initializer list { 2, 6, -1, 3, 5, 0, 0, 0 }
   \endcode

// Since an InitializerVector represents a specific initializer list, its lifetime is bound to
// the lifetime of the according initializer list. When the initializer list goes out of scope
// access to the initializer list via an InitializerVector results in undefined behavior:

   \code
   blaze::InitializerVector<int> d{ 1, 2, 3, 4, 5 };     // Undefined behavior!
   blaze::InitializerVector<int> e( { 0, 3, 2 }, 3UL );  // Undefined behavior!
   \endcode

// Also, an InitializerVector can only be used on the right of an assignment as its elements
// are considered to be immutable. The following example gives an impression on the usage of an
// InitializerVector:

   \code
   const auto list = { 2, 6, -1, 3, 5 };

   blaze::InitializerVector<int> f( list );  // Representation of the initializer list as dense vector
   blaze::DynamicVector<int> g;

   g = f;  // Initialize vector g via vector f
   f = g;  // Compilation error! Cannot assign to an initializer vector
   \endcode

// An initializer vector can be used as operand in arithmetic operations. All operations (addition,
// subtraction, multiplication, scaling, ...) can be performed on all possible combinations of
// dense and sparse vectors with fitting element types. The following example gives an impression
// of the use of InitializerVector:

   \code
   using blaze::InitializerVector;
   using blaze::DynamicVector;
   using blaze::CompressedVector;
   using blaze::DynamicMatrix;

   const auto list = { 1.0, 2.0 };

   InitializerVector<double> a( list );

   DynamicVector<double>   b( 2, 2.0 );  // Directly, homogeneously initialized 2D vector
   CompressedVector<float> c( 2 );       // Empty sparse single precision vector
   DynamicVector<double>   d;            // Default constructed dynamic vector
   DynamicMatrix<double>   A;            // Default constructed row-major matrix

   d = a + b;    // Vector addition between vectors of equal element type
   d = a - c;    // Vector subtraction between a dense and sparse vector with different element types
   d = a * b;    // Component-wise vector multiplication

   d = a * 2.0;  // Scaling of vector a
   d = 2.0 * a;  // Scaling of vector a

   d += a - b;   // Addition assignment
   d -= a + c;   // Subtraction assignment
   d *= a * b;   // Multiplication assignment

   double scalar = trans( a ) * b;  // Scalar/dot/inner product between two vectors

   A = a * trans( b );  // Outer product between two vectors
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
class InitializerVector
   : public DenseVector< InitializerVector<Type,TF,Tag>, TF >
{
 public:
   //**Type definitions****************************************************************************
   using This     = InitializerVector<Type,TF,Tag>;  //!< Type of this InitializerVector instance.
   using BaseType = DenseVector<This,TF>;            //!< Base type of this InitializerVector instance.

   //! Result type for expression template evaluations.
   using ResultType = DynamicVector<Type,TF>;

   //! Transpose type for expression template evaluations.
   using TransposeType = DynamicVector<Type,!TF>;

   using ElementType   = Type;                     //!< Type of the vector elements.
   using TagType       = Tag;                      //!< Tag type of this InitializerVector instance.
   using ReturnType    = const Type&;              //!< Return type for expression template evaluations
   using CompositeType = const InitializerVector&; //!< Data type for composite expression templates.

   using Reference      = const Type&;  //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;  //!< Reference to a constant vector value.
   using Pointer        = const Type*;  //!< Pointer to a non-constant vector value.
   using ConstPointer   = const Type*;  //!< Pointer to a constant vector value.

   using Iterator      = InitializerIterator<Type>;  //!< Iterator over non-constant elements.
   using ConstIterator = InitializerIterator<Type>;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a InitializerVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = InitializerVector<NewType,TF,Tag>;  //!< The type of the other InitializerVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a InitializerVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = InitializerVector<Type,TF,Tag>;  //!< The type of the other InitializerVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SIMD optimization.
   /*! The \a simdEnabled compilation flag indicates whether expressions the vector is involved
       in can be optimized via SIMD operationss. In case the element type of the vector is a
       vectorizable data type, the \a simdEnabled compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool simdEnabled = false;

   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline InitializerVector( initializer_list<Type> list ) noexcept;
   inline InitializerVector( initializer_list<Type> list, size_t n );

   InitializerVector( const InitializerVector& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~InitializerVector() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline ConstReference at( size_t index ) const;
   inline ConstPointer   data  () const noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   InitializerVector& operator=( const InitializerVector& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size() const noexcept;
   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   swap( InitializerVector& v ) noexcept;
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
   using ListType = initializer_list<Type>;  //!< Type of the represented initializer list.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;    //!< The current size/dimension of the vector.
   ListType list_;  //!< The initializer list represented by the vector.
                    /*!< Access to the vector elements is gained via the subscript operator.
                         The order of the elements is
                         \f[\left(\begin{array}{*{5}{c}}
                         0 & 1 & 2 & \cdots & N-1 \\
                         \end{array}\right)\f] */

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

template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
const Type InitializerVector<Type,TF,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constructor for InitializerVector.
//
// \param list The initializer list represented by the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline InitializerVector<Type,TF,Tag>::InitializerVector( initializer_list<Type> list ) noexcept
   : size_( list.size() )  // The current size/dimension of the vector
   , list_( list )         // The initializer list represented by the vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for InitializerVector.
//
// \param list The initializer list represented by the vector.
// \param n The size of the vector.
// \exception std::invalid_argument Invalid initializer list dimension.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline InitializerVector<Type,TF,Tag>::InitializerVector( initializer_list<Type> list, size_t n )
   : size_( n    )  // The current size/dimension of the vector
   , list_( list )  // The initializer list represented by the vector
{
   if( n < list.size() ) {
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
inline typename InitializerVector<Type,TF,Tag>::ConstReference
   InitializerVector<Type,TF,Tag>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid vector access index" );
   if( index < list_.size() )
      return list_.begin()[index];
   else
      return zero_;
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
inline typename InitializerVector<Type,TF,Tag>::ConstReference
   InitializerVector<Type,TF,Tag>::at( size_t index ) const
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
// This function returns a pointer to the internal storage of the initializer vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename InitializerVector<Type,TF,Tag>::ConstPointer
   InitializerVector<Type,TF,Tag>::data() const noexcept
{
   return list_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the initializer vector.
//
// \return Iterator to the first element of the initializer vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename InitializerVector<Type,TF,Tag>::ConstIterator
   InitializerVector<Type,TF,Tag>::begin() const noexcept
{
   return ConstIterator( 0UL, list_ );;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the initializer vector.
//
// \return Iterator to the first element of the initializer vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename InitializerVector<Type,TF,Tag>::ConstIterator
   InitializerVector<Type,TF,Tag>::cbegin() const noexcept
{
   return ConstIterator( 0UL, list_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the initializer vector.
//
// \return Iterator just past the last element of the initializer vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename InitializerVector<Type,TF,Tag>::ConstIterator
   InitializerVector<Type,TF,Tag>::end() const noexcept
{
   return ConstIterator( size_, list_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the initializer vector.
//
// \return Iterator just past the last element of the initializer vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename InitializerVector<Type,TF,Tag>::ConstIterator
   InitializerVector<Type,TF,Tag>::cend() const noexcept
{
   return ConstIterator( size_, list_ );
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
inline size_t InitializerVector<Type,TF,Tag>::size() const noexcept
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
inline size_t InitializerVector<Type,TF,Tag>::spacing() const noexcept
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
inline size_t InitializerVector<Type,TF,Tag>::capacity() const noexcept
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
inline size_t InitializerVector<Type,TF,Tag>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( size_t i=0UL; i<list_.size(); ++i ) {
      if( !isDefault<strict>( list_.begin()[i] ) )
         ++nonzeros;
   }

   return nonzeros;
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
inline void InitializerVector<Type,TF,Tag>::swap( InitializerVector& v ) noexcept
{
   using std::swap;

   swap( size_, v.size_ );
   swap( list_, v.list_ );
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
inline bool InitializerVector<Type,TF,Tag>::canAlias( const Other* alias ) const noexcept
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
inline bool InitializerVector<Type,TF,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************




//=================================================================================================
//
//  INITIALIZERVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name InitializerVector operators */
//@{
template< typename Type, bool TF, typename Tag >
bool isIntact( const InitializerVector<Type,TF,Tag>& v ) noexcept;

template< typename Type, bool TF, typename Tag >
void swap( InitializerVector<Type,TF,Tag>& a, InitializerVector<Type,TF,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given initializer vector are intact.
// \ingroup initializer_vector
//
// \param v The initializer vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the initializer vector are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::InitializerVector<int> a{ 1, 2, 0, 3 };
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline bool isIntact( const InitializerVector<Type,TF,Tag>& v ) noexcept
{
   MAYBE_UNUSED( v );

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two vectors.
// \ingroup initializer_vector
//
// \param a The first vector to be swapped.
// \param b The second vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void swap( InitializerVector<Type,TF,Tag>& a, InitializerVector<Type,TF,Tag>& b ) noexcept
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
template< typename T, bool TF, typename Tag >
struct HasConstDataAccess< InitializerVector<T,TF,Tag> >
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
template< typename T, bool TF, typename Tag >
struct IsInitializer< InitializerVector<T,TF,Tag> >
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
template< typename T1, bool TF, typename Tag, typename T2 >
struct HighType< InitializerVector<T1,TF,Tag>, InitializerVector<T2,TF,Tag> >
{
   using Type = InitializerVector< typename HighType<T1,T2>::Type, TF, Tag >;
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
struct LowType< InitializerVector<T1,TF,Tag>, InitializerVector<T2,TF,Tag> >
{
   using Type = InitializerVector< typename LowType<T1,T2>::Type, TF, Tag >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
