//=================================================================================================
/*!
//  \file blaze/util/AlignedAllocator.h
//  \brief Header file for the AlignedAllocator implementation
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

#ifndef _BLAZE_UTIL_ALIGNEDALLOCATOR_H_
#define _BLAZE_UTIL_ALIGNEDALLOCATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/AlignmentOf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Allocator for type-specific aligned memory.
// \ingroup util
//
// The AlignedAllocator class template represents an implementation of the allocator concept of
// the standard library for the allocation of type-specific, aligned, uninitialized memory. The
// allocator performs its allocation via the blaze::alignedAllocate() and blaze::alignedDeallocate()
// functions to guarantee properly aligned memory based on the alignment restrictions of the
// specified type \a T. For instance, in case the given type is a fundamental, built-in data type
// and in case SSE vectorization is possible, the returned memory is guaranteed to be at least
// 16-byte aligned. In case AVX is active, the memory is even guaranteed to be at least 32-byte
// aligned.
*/
template< typename T >
class AlignedAllocator
{
 public:
   //**Type definitions****************************************************************************
   using ValueType      = T;          //!< Type of the allocated values.
   using SizeType       = size_t;     //!< Size type of the aligned allocator.
   using DifferenceType = ptrdiff_t;  //!< Difference type of the aligned allocator.

   // STL allocator requirements
   using value_type      = ValueType;       //!< Type of the allocated values.
   using size_type       = SizeType;        //!< Size type of the aligned allocator.
   using difference_type = DifferenceType;  //!< Difference type of the aligned allocator.
   //**********************************************************************************************

   //**rebind class definition*********************************************************************
   /*!\brief Implementation of the AlignedAllocator rebind mechanism.
   */
   template< typename U >
   struct rebind
   {
      using other = AlignedAllocator<U>;  //!< Type of the other allocator.
   };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   AlignedAllocator() = default;

   template< typename U >
   inline AlignedAllocator( const AlignedAllocator<U>& );
   //@}
   //**********************************************************************************************

   //**Allocation functions************************************************************************
   /*!\name Allocation functions */
   //@{
   inline T*   allocate  ( size_t numObjects );
   inline void deallocate( T* ptr, size_t numObjects ) noexcept;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion constructor from different AlignedAllocator instances.
//
// \param allocator The foreign aligned allocator to be copied.
*/
template< typename T >
template< typename U >
inline AlignedAllocator<T>::AlignedAllocator( const AlignedAllocator<U>& allocator )
{
   MAYBE_UNUSED( allocator );
}
//*************************************************************************************************




//=================================================================================================
//
//  ALLOCATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Allocates aligned memory for the specified number of objects.
//
// \param numObjects The number of objects to be allocated.
// \return Pointer to the newly allocated memory.
//
// This function allocates a junk of memory for the specified number of objects of type \a T.
// The returned pointer is guaranteed to be aligned according to the alignment restrictions of
// the data type \a T. For instance, in case the type is a fundamental, built-in data type and
// in case SSE vectorization is possible, the returned memory is guaranteed to be at least
// 16-byte aligned. In case AVX is active, the memory is even guaranteed to be 32-byte aligned.
*/
template< typename T >
inline T* AlignedAllocator<T>::allocate( size_t numObjects )
{
   const size_t alignment( AlignmentOf_v<T> );

   if( alignment >= 8UL ) {
      return reinterpret_cast<T*>( alignedAllocate( numObjects*sizeof(T), alignment ) );
   }
   else {
      return static_cast<T*>( operator new[]( numObjects * sizeof( T ) ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Deallocation of memory.
//
// \param ptr The address of the first element of the array to be deallocated.
// \param numObjects The number of objects to be deallocated.
// \return void
//
// This function deallocates a junk of memory that was previously allocated via the allocate()
// function. Note that the argument \a numObjects must be equal ot the first argument of the call
// to allocate() that origianlly produced \a ptr.
*/
template< typename T >
inline void AlignedAllocator<T>::deallocate( T* ptr, size_t numObjects ) noexcept
{
   MAYBE_UNUSED( numObjects );

   if( ptr == nullptr )
      return;

   const size_t alignment( AlignmentOf_v<T> );

   if( alignment >= 8UL ) {
      alignedDeallocate( ptr );
   }
   else {
      operator delete[]( ptr );
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name AlignedAllocator operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const AlignedAllocator<T1>& lhs, const AlignedAllocator<T2>& rhs ) noexcept;

template< typename T1, typename T2 >
inline bool operator!=( const AlignedAllocator<T1>& lhs, const AlignedAllocator<T2>& rhs ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two AlignedAllocator objects.
//
// \param lhs The left-hand side aligned allocator.
// \param rhs The right-hand side aligned allocator.
// \return \a true.
*/
template< typename T1    // Type of the left-hand side aligned allocator
        , typename T2 >  // Type of the right-hand side aligned allocator
inline bool operator==( const AlignedAllocator<T1>& lhs, const AlignedAllocator<T2>& rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two AlignedAllocator objects.
//
// \param lhs The left-hand side aligned allocator.
// \param rhs The right-hand side aligned allocator.
// \return \a false.
*/
template< typename T1    // Type of the left-hand side aligned allocator
        , typename T2 >  // Type of the right-hand side aligned allocator
inline bool operator!=( const AlignedAllocator<T1>& lhs, const AlignedAllocator<T2>& rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );
   return false;
}
//*************************************************************************************************

} // namespace blaze

#endif
