//=================================================================================================
/*!
//  \file blaze/util/SmallArray.h
//  \brief Header file for the SmallArray implementation
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

#ifndef _BLAZE_UTIL_SMALLARRAY_H_
#define _BLAZE_UTIL_SMALLARRAY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <memory>
#include <blaze/util/algorithms/Destroy.h>
#include <blaze/util/algorithms/DestroyAt.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/UninitializedDefaultConstruct.h>
#include <blaze/util/algorithms/UninitializedMove.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/Exception.h>
#include <blaze/util/InitializerList.h>
#include <blaze/util/smallarray/SmallArrayData.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConstructible.h>
#include <blaze/util/typetraits/IsAssignable.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of a dynamic array with small array optimization.
// \ingroup util
//
// The SmallArray class template is a hybrid data structure between a static array and a dynamic
// array. It provides static, in-place memory for up to \a N elements of type \a T, but can grow
// beyond this size by allocating dynamic memory via its allocator of type \a A.
*/
template< typename T                        // Data type of the elements
        , size_t N                          // Number of preallocated elements
        , typename A = std::allocator<T> >  // Type of the allocator
class SmallArray
   : private SmallArrayData<T,N>
   , private A
{
 public:
   //**Type definitions****************************************************************************
   using ElementType    = T;         //!< Type of the array elements.
   using Pointer        = T*;        //!< Pointer to a non-constant array element.
   using ConstPointer   = const T*;  //!< Pointer to a constant array element.
   using Reference      = T&;        //!< Reference to a non-constant array element.
   using ConstReference = const T&;  //!< Reference to a constant array element.
   using Iterator       = T*;        //!< Iterator over non-constant elements.
   using ConstIterator  = const T*;  //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline SmallArray( const A& alloc = A() );
   explicit inline SmallArray( size_t n, const A& alloc = A() );
            inline SmallArray( size_t n, const T& init, const A& alloc = A() );

   template< typename InputIt >
   inline SmallArray( InputIt first, InputIt last, const A& alloc = A() );

   template< typename U >
   inline SmallArray( initializer_list<U> list, const A& alloc = A() );

   inline SmallArray( const SmallArray& sa );
   inline SmallArray( SmallArray&& sa );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~SmallArray();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data() noexcept;
   inline ConstPointer   data() const noexcept;
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
   template< typename U >
   inline SmallArray& operator=( initializer_list<U> list );

   inline SmallArray& operator=( const SmallArray& rhs );
   inline SmallArray& operator=( SmallArray&& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool   empty() const noexcept;
   inline size_t size() const noexcept;
   inline size_t capacity() const noexcept;

   inline void     clear();
          void     resize( size_t n );
          void     resize( size_t n, const T& value );
          void     reserve( size_t n );
          void     shrinkToFit();
          void     pushBack( const T& value );
          void     pushBack( T&& value );
          Iterator insert( Iterator pos, const T& value );
          Iterator insert( Iterator pos, T&& value );
          Iterator erase( Iterator pos );
          Iterator erase( Iterator first, Iterator last );
          void     swap( SmallArray& sa ) noexcept( IsNothrowMoveConstructible_v<T> );
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
   inline SmallArray( size_t n, const A& alloc, Uninitialized );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using A::allocate;
   using A::deallocate;

   inline bool isDynamic() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   T* begin_;  //!< Pointer to the beginning of the currently used storage.
   T* end_;    //!< Pointer to the end of the currently used storage.
   T* final_;  //!< Pointer to the very end of the currently used storage.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST   ( T );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE( T );
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
/*!\brief The (default) constructor for SmallArray.
//
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( const A& alloc )
   : SmallArray( 0UL, alloc, Uninitialized() )
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an array of size \a n. No element initialization is performed!
//
// \param n The initial size of the array.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( size_t n, const A& alloc )
   : SmallArray( n, alloc, Uninitialized() )
{
   std::uninitialized_fill( begin_, end_, T() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for an array of size \a n.
//
// \param n The initial size of the array.
// \param init The initial value of the array elements.
// \param alloc Allocator for all memory allocations of this container.
//
// All array elements are initialized with the specified value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( size_t n, const T& init, const A& alloc )
   : SmallArray( n, alloc, Uninitialized() )
{
   std::uninitialized_fill( begin_, end_, init );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a range of elements.
//
// \param first Iterator to the be first element of the input range.
// \param last Iterator to the element one past the last element of the input range.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T          // Data type of the elements
        , size_t N            // Number of preallocated elements
        , typename A >        // Type of the allocator
template< typename InputIt >  // Type of the iterators
inline SmallArray<T,N,A>::SmallArray( InputIt first, InputIt last, const A& alloc )
   : SmallArray( std::distance( first, last ), alloc, Uninitialized() )
{
   std::uninitialized_copy( first, last, begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all array elements.
//
// \param list The initializer list.
// \param alloc Allocator for all memory allocations of this container.
//
// This constructor provides the option to explicitly initialize the elements of the small array
// within a constructor call:

   \code
   blaze::SmallArray<double,8UL> v1{ 4.2, 6.3, -1.2 };
   \endcode

// The array is sized according to the size of the initializer list and all its elements are
// copy initialized by the elements of the given initializer list.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
template< typename U >  // Type of the initializer list elements
inline SmallArray<T,N,A>::SmallArray( initializer_list<U> list, const A& alloc )
   : SmallArray( list.size(), alloc, Uninitialized() )
{
   std::uninitialized_copy( list.begin(), list.end(), begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for SmallArray.
//
// \param sa The small array to be copied.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( const SmallArray& sa )
   : SmallArray( sa.size(), A(), Uninitialized() )
{
   std::uninitialized_copy( sa.begin(), sa.end(), begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for SmallArray.
//
// \param sa The small array to be moved into this instance.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( SmallArray&& sa )
   // Base class initialization is intentionally omitted
   : begin_( sa.begin_ )  // Pointer to the beginning of the currently used storage
   , end_  ( sa.end_   )  // Pointer to the end of the currently used storage
   , final_( sa.final_ )  // Pointer to the very end of the currently used storage
{
   if( !sa.isDynamic() ) {
      begin_ = this->array();
      end_   = begin_ + sa.size();
      final_ = begin_ + N;
      blaze::uninitialized_move( sa.begin_, sa.end_, begin_ );
      blaze::destroy( sa.begin_, sa.end_ );
   }

   sa.begin_ = sa.array();
   sa.end_   = sa.begin_;
   sa.final_ = sa.begin_ + N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary constructor for SmallArray.
//
// \param n The initial size of the array.
// \param alloc Allocator for all memory allocations of this container.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::SmallArray( size_t n, const A& alloc, Uninitialized )
   // SmallArrayData is initialized
   // begin_ is intentionally not initialized
   // end_ is intentionally not initialized
   // final_ is intentionally not initialized
   : A( alloc )  // Base class initialization
{
   if( n <= N ) {
      begin_ = this->array();
      end_   = begin_ + n;
      final_ = begin_ + N;
   }
   else {
      begin_ = allocate( n );
      end_   = begin_ + n;
      final_ = begin_ + n;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for SmallArray.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>::~SmallArray()
{
   blaze::destroy( begin_, end_ );

   if( isDynamic() ) {
      deallocate( begin_, capacity() );
   }
}
//*************************************************************************************************





//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::Reference
   SmallArray<T,N,A>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size(), "Invalid small array access index" );
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference-to-const to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstReference
   SmallArray<T,N,A>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size(), "Invalid small array access index" );
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid small array access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::Reference
   SmallArray<T,N,A>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid small array access index" );
   }
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the array elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid small array access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstReference
   SmallArray<T,N,A>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid small array access index" );
   }
   return begin_[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the array elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::Pointer
   SmallArray<T,N,A>::data() noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level data access to the array elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstPointer
   SmallArray<T,N,A>::data() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small array.
//
// \return Iterator to the first element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::begin() noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small array.
//
// \return Iterator to the first element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstIterator
   SmallArray<T,N,A>::begin() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the small array.
//
// \return Iterator to the first element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstIterator
   SmallArray<T,N,A>::cbegin() const noexcept
{
   return begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small array.
//
// \return Iterator just past the last element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::end() noexcept
{
   return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small array.
//
// \return Iterator just past the last element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstIterator
   SmallArray<T,N,A>::end() const noexcept
{
   return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the small array.
//
// \return Iterator just past the last element of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline typename SmallArray<T,N,A>::ConstIterator
   SmallArray<T,N,A>::cend() const noexcept
{
   return end_;
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief List assignment to all array elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the small
// array by means of an initializer list:

   \code
   blaze::SmallArray<double,8UL> v;
   v = { 4.2, 6.3, -1.2 };
   \endcode

// The array is resized according to the size of the initializer list and all its elements are
// assigned the values from the given initializer list.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
template< typename U >  // Type of the initializer list elements
inline SmallArray<T,N,A>& SmallArray<T,N,A>::operator=( initializer_list<U> list )
{
   resize( list.size() );
   std::copy( list.begin(), list.end(), begin_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for SmallArray.
//
// \param rhs Array to be copied.
// \return Reference to the assigned array.
//
// The array is resized according to the given small array and initialized as a copy of this
// array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>& SmallArray<T,N,A>::operator=( const SmallArray& rhs )
{
   resize( rhs.size() );
   std::copy( rhs.begin(), rhs.end(), begin_ );

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for SmallArray.
//
// \param rhs The array to be moved into this instance.
// \return Reference to the assigned array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline SmallArray<T,N,A>& SmallArray<T,N,A>::operator=( SmallArray&& rhs )
{
   resize( rhs.size() );
   std::move( rhs.begin_, rhs.end_, begin_ );
   blaze::destroy( rhs.begin_, rhs.end_ );

   rhs.begin_ = rhs.array();
   rhs.end_   = rhs.begin_;
   rhs.final_ = rhs.begin_ + N;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the array is empty.
//
// \return \a true in case the array is empty, \a false if not.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline bool SmallArray<T,N,A>::empty() const noexcept
{
   return begin_ == end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the current size/dimension of the small array.
//
// \return The current size of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline size_t SmallArray<T,N,A>::size() const noexcept
{
   return end_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the small array.
//
// \return The capacity of the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline size_t SmallArray<T,N,A>::capacity() const noexcept
{
   return final_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the array.
//
// \return void
//
// After the clear() function, the size of the array is 0.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void SmallArray<T,N,A>::clear()
{
   blaze::destroy( begin_, end_ );

   if( isDynamic() ) {
      deallocate( begin_, capacity() );
   }

   begin_ = this->array();
   end_   = begin_;
   final_ = begin_ + N;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the array.
//
// \param n The new size of the array.
// \return void
//
// This function resizes the array using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the array is too small. New array
// elements are not initialized!
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::resize( size_t n )
{
   if( n > size() )
   {
      reserve( n );
      blaze::uninitialized_default_construct( begin_+size(), begin_+n );
      end_ = begin_ + n;
   }
   else if( n < size() )
   {
      blaze::destroy( begin_+n, end_ );
      end_ = begin_ + n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the array.
//
// \param n The new size of the array.
// \param value The initial value of new array elements.
// \return void
//
// This function resizes the array using the given size to \a n. During this operation, new
// dynamic memory may be allocated in case the capacity of the array is too small. New array
// elements are initialized to \a value!
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::resize( size_t n, const T& value )
{
   if( n > size() )
   {
      reserve( n );
      std::uninitialized_fill( begin_+size(), begin_+n, value );
      end_ = begin_ + n;
   }
   else if( n < size() )
   {
      blaze::destroy( begin_+n, end_ );
      end_ = begin_ + n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the array.
//
// \param n The new minimum capacity of the array.
// \return void
//
// This function increases the capacity of the array to at least \a n elements. The current
// values of the array elements are preserved.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::reserve( size_t n )
{
   const size_t oldCapacity( capacity() );

   if( n > oldCapacity )
   {
      const size_t oldSize( size() );
      T* tmp( allocate( n ) );

      if( IsNothrowMoveConstructible_v<T> ) {
         blaze::uninitialized_move( begin_, end_, tmp );
      }
      else {
         std::uninitialized_copy( begin_, end_, tmp );
      }

      blaze::destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, oldCapacity );
      }

      final_ = tmp + n;
      end_   = tmp + oldSize;
      begin_ = tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the array by removing unused capacity. Please note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this array are invalidated.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::shrinkToFit()
{
   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( isDynamic() && oldCapacity > oldSize )
   {
      T* tmp( allocate( oldSize ) );

      if( IsNothrowMoveConstructible_v<T> ) {
         blaze::uninitialized_move( begin_, end_, tmp );
      }
      else {
         std::uninitialized_copy( begin_, end_, tmp );
      }

      blaze::destroy( begin_, end_ );
      deallocate( begin_, oldCapacity );

      final_ = tmp + oldSize;
      end_   = final_;
      begin_ = tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adding an element to the end of the small array.
//
// \param value The element to be added to the end of the small array.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::pushBack( const T& value )
{
   using blaze::max;

   const size_t oldCapacity( capacity() );

   if( size() == oldCapacity ) {
      reserve( max( 2UL*oldCapacity, 7UL ) );
   }

   ::new ( end_ ) T( value );
   ++end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Adding an element to the end of the small array.
//
// \param value The element to be added to the end of the small array.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::pushBack( T&& value )
{
   using blaze::max;

   const size_t oldCapacity( capacity() );

   if( size() == oldCapacity ) {
      reserve( max( 2UL*oldCapacity, 7UL ) );
   }

   ::new ( end_ ) T( std::move( value ) );
   ++end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element at the specified position into the small array.
//
// \param pos The position of the new element.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::insert( Iterator pos, const T& value )
{
   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( oldSize == oldCapacity )
   {
      const size_t newCapacity( max( 1UL, 2UL*oldCapacity ) );
      const size_t index( pos - begin_ );

      T* tmp   ( allocate( newCapacity ) );
      T* newpos( tmp + index );

      blaze::uninitialized_move( begin_, pos, tmp );
      ::new ( newpos ) T( value );
      blaze::uninitialized_move( pos, end_, tmp+index+1UL );
      blaze::destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, capacity() );
      }

      final_ = tmp + newCapacity;
      end_   = tmp + oldSize + 1UL;
      begin_ = tmp;

      return newpos;
   }
   else if( pos == end_ )
   {
      ::new( pos ) T( value );
      ++end_;
      return pos;
   }
   else
   {
      const auto tmp( end_ - 1UL );
      ::new ( end_ ) T( std::move( *tmp ) );

      try {
         std::move_backward( pos, tmp, end_ );
         blaze::destroy_at( pos );
         ::new ( pos ) T( value );
         ++end_;
      }
      catch( ... ) {
         blaze::destroy_at( end_ );
         throw;
      }

      return pos;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element at the specified position into the small array.
//
// \param pos The position of the new element.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::insert( Iterator pos, T&& value )
{
   const size_t oldCapacity( capacity() );
   const size_t oldSize    ( size() );

   if( oldSize == oldCapacity )
   {
      const size_t newCapacity( max( 1UL, 2UL*oldCapacity ) );
      const size_t index( pos - begin_ );

      T* tmp   ( allocate( newCapacity ) );
      T* newpos( tmp + index );

      blaze::uninitialized_move( begin_, pos, tmp );
      ::new ( newpos ) T( std::move( value ) );
      blaze::uninitialized_move( pos, end_, tmp+index+1UL );
      blaze::destroy( begin_, end_ );

      if( isDynamic() ) {
         deallocate( begin_, capacity() );
      }

      final_ = tmp + newCapacity;
      end_   = tmp + oldSize + 1UL;
      begin_ = tmp;

      return newpos;
   }
   else if( pos == end_ )
   {
      ::new( pos ) T( std::move( value ) );
      ++end_;
      return pos;
   }
   else
   {
      const auto tmp( end_ - 1UL );
      ::new ( end_ ) T( std::move( *tmp ) );

      try {
         std::move_backward( pos, tmp, end_ );
         blaze::destroy_at( pos );
         ::new ( pos ) T( std::move( value ) );
         ++end_;
      }
      catch( ... ) {
         blaze::destroy_at( end_ );
         throw;
      }

      return pos;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the small array.
//
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::erase( Iterator pos )
{
   std::move( pos+1UL, end_, pos );
   --end_;
   blaze::destroy_at( end_ );

   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the small array.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the small array.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
typename SmallArray<T,N,A>::Iterator
   SmallArray<T,N,A>::erase( Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( first <= last, "Invalid range detected" );

   const size_t n( last - first );

   std::move( last, end_, first );
   end_ -= n;
   blaze::destroy( end_, end_+n );

   return first;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two small arrays.
//
// \param sa The small array to be swapped.
// \return void
//
// This function swaps the contents of two small arrays. Please note that this function is only
// guaranteed to not throw an exception if the move constructor of the underlying data type \a T
// is guaranteed to be noexcept.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
void SmallArray<T,N,A>::swap( SmallArray& sa ) noexcept( IsNothrowMoveConstructible_v<T> )
{
   using std::swap;

   if( isDynamic() && sa.isDynamic() )
   {
      swap( begin_, sa.begin_ );
      swap( end_  , sa.end_   );
      swap( final_, sa.final_ );
   }
   else if( isDynamic() )
   {
      const size_t n( sa.size() );

      blaze::uninitialized_move( sa.begin_, sa.end_, this->array() );
      blaze::destroy( sa.begin_, sa.end_ );

      sa.begin_ = begin_;
      sa.end_   = end_;
      sa.final_ = final_;

      begin_ = this->array();
      end_   = begin_ + n;
      final_ = begin_ + N;
   }
   else if( sa.isDynamic() )
   {
      const size_t n( size() );

      blaze::uninitialized_move( begin_, end_, sa.array() );
      blaze::destroy( begin_, end_ );

      begin_ = sa.begin_;
      end_   = sa.end_;
      final_ = sa.final_;

      sa.begin_ = sa.array();
      sa.end_   = sa.begin_ + n;
      sa.final_ = sa.begin_ + N;
   }
   else if( size() > sa.size() )
   {
      const size_t n( size() - sa.size() );
      const auto pos = std::swap_ranges( sa.begin_, sa.end_, begin_ );

      blaze::uninitialized_move( pos, end_, sa.end_ );
      blaze::destroy( pos, end_ );
      end_    -= n;
      sa.end_ += n;
   }
   else
   {
      const size_t n( sa.size() - size() );
      const auto pos = std::swap_ranges( begin_, end_, sa.begin_ );

      blaze::uninitialized_move( pos, sa.end_, end_ );
      blaze::destroy( pos, sa.end_ );
      end_    += n;
      sa.end_ -= n;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the small array uses its dynamic storage.
//
// \return \a true in case the dynamic storage is in use, \a false if not.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline bool SmallArray<T,N,A>::isDynamic() const noexcept
{
   return ( N == 0UL || begin_ != this->array() );
}
//*************************************************************************************************




//=================================================================================================
//
//  SMALLARRAY OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SmallArray operators */
//@{
template< typename T1, size_t N1, typename A1, typename T2, size_t N2, typename A2 >
inline bool operator==( const SmallArray<T1,N1,A1>& lhs, const SmallArray<T2,N2,A2>& rhs );

template< typename T1, size_t N1, typename A1, typename T2, size_t N2, typename A2 >
inline bool operator!=( const SmallArray<T1,N1,A1>& lhs, const SmallArray<T2,N2,A2>& rhs );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::Iterator begin( SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator begin( const SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator cbegin( const SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::Iterator end( SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator end( const SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator cend( const SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline void clear( SmallArray<T,N,A>& sa );

template< typename T, size_t N, typename A >
inline void swap( SmallArray<T,N,A>& a, SmallArray<T,N,A>& b )
   noexcept( IsNothrowMoveConstructible_v<T> );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense arrays.
// \ingroup util
//
// \param lhs The left-hand side small array for the comparison.
// \param rhs The right-hand side small array for the comparison.
// \return \a true if the two arrays are equal, \a false if not.
*/
template< typename T1    // Data type of the elements of the left-hand side array
        , size_t N1      // Number of elements of the left-hand side array
        , typename A1    // Type of the allocator of the left-hand side array
        , typename T2    // Data type of the elements of the right-hand side array
        , size_t N2      // Number of elements of the right-hand side array
        , typename A2 >  // Type of the allocator of the right-hand side array
inline bool operator==( const SmallArray<T1,N1,A1>& lhs, const SmallArray<T2,N2,A2>& rhs )
{
   if( lhs.size() != rhs.size() ) return false;

   return std::equal( lhs.begin(), lhs.end(), rhs.begin() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense arrays.
// \ingroup util
//
// \param lhs The left-hand side small array for the comparison.
// \param rhs The right-hand side small array for the comparison.
// \return \a true if the two arrays are not equal, \a false if they are equal.
*/
template< typename T1    // Data type of the elements of the left-hand side array
        , size_t N1      // Number of elements of the left-hand side array
        , typename A1    // Type of the allocator of the left-hand side array
        , typename T2    // Data type of the elements of the right-hand side array
        , size_t N2      // Number of elements of the right-hand side array
        , typename A2 >  // Type of the allocator of the right-hand side array
inline bool operator!=( const SmallArray<T1,N1,A1>& lhs, const SmallArray<T2,N2,A2>& rhs )
{
   return !( lhs == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator to the first element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::Iterator begin( SmallArray<T,N,A>& sa )
{
   return sa.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator to the first element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator begin( const SmallArray<T,N,A>& sa )
{
   return sa.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator to the first element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator cbegin( const SmallArray<T,N,A>& sa )
{
   return sa.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator just past the last element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::Iterator end( SmallArray<T,N,A>& sa )
{
   return sa.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator just past the last element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator end( const SmallArray<T,N,A>& sa )
{
   return sa.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the given small array.
// \ingroup util
//
// \param sa The given small array.
// \return Iterator just past the last element of the given small array.
*/
template< typename T, size_t N, typename A >
inline typename SmallArray<T,N,A>::ConstIterator cend( const SmallArray<T,N,A>& sa )
{
   return sa.end();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given small array.
// \ingroup util
//
// \param sa The small array to be cleared.
// \return void
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void clear( SmallArray<T,N,A>& sa )
{
   sa.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two small arrays.
// \ingroup util
//
// \param a The first array to be swapped.
// \param b The second array to be swapped.
// \return void
//
// This function swaps the contents of two small arrays. Please note that this function is only
// guaranteed to not throw an exception if the move constructor of the underlying data type \a T
// is guaranteed to be noexcept.
*/
template< typename T    // Data type of the elements
        , size_t N      // Number of preallocated elements
        , typename A >  // Type of the allocator
inline void swap( SmallArray<T,N,A>& a, SmallArray<T,N,A>& b )
   noexcept( IsNothrowMoveConstructible_v<T> )
{
   a.swap( b );
}
//*************************************************************************************************

} // namespace blaze

#endif
