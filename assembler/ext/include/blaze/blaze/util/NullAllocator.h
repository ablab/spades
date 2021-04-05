//=================================================================================================
/*!
//  \file blaze/util/NullAllocator.h
//  \brief Header file for the NullAllocator implementation
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

#ifndef _BLAZE_UTIL_NULLALLOCATOR_H_
#define _BLAZE_UTIL_NULLALLOCATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Allocator returning nullptr.
// \ingroup util
//
// The NullAllocator class template represents an implementation of the allocator concept of
// the standard library and acts as a stand-in for situation where no memory allocation is
// required. The NullAllocator will never allocate any memory and will always return nullptr.
*/
template< typename T >
class NullAllocator
{
 public:
   //**Type definitions****************************************************************************
   using ValueType      = T;          //!< Type of the allocated values.
   using SizeType       = size_t;     //!< Size type of the null allocator.
   using DifferenceType = ptrdiff_t;  //!< Difference type of the null allocator.

   // STL allocator requirements
   using value_type      = ValueType;       //!< Type of the allocated values.
   using size_type       = SizeType;        //!< Size type of the null allocator.
   using difference_type = DifferenceType;  //!< Difference type of the null allocator.
   //**********************************************************************************************

   //**rebind class definition*********************************************************************
   /*!\brief Implementation of the NullAllocator rebind mechanism.
   */
   template< typename U >
   struct rebind
   {
      using other = NullAllocator<U>;  //!< Type of the other allocator.
   };
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   NullAllocator() = default;

   template< typename U >
   inline NullAllocator( const NullAllocator<U>& );
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
/*!\brief Conversion constructor from different NullAllocator instances.
//
// \param allocator The foreign null allocator to be copied.
*/
template< typename T >
template< typename U >
inline NullAllocator<T>::NullAllocator( const NullAllocator<U>& allocator )
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
/*!\brief Performs no memory allocation and returns nullptr.
//
// \param numObjects The number of objects to be allocated.
// \return Pointer to the newly allocated memory.
//
// This function does not perform any memory allocation and always returns nullptr.
*/
template< typename T >
inline T* NullAllocator<T>::allocate( size_t numObjects )
{
   MAYBE_UNUSED( numObjects );

   return static_cast<T*>( nullptr );
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
inline void NullAllocator<T>::deallocate( T* ptr, size_t numObjects ) noexcept
{
   MAYBE_UNUSED( ptr, numObjects );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name NullAllocator operators */
//@{
template< typename T1, typename T2 >
inline bool operator==( const NullAllocator<T1>& lhs, const NullAllocator<T2>& rhs ) noexcept;

template< typename T1, typename T2 >
inline bool operator!=( const NullAllocator<T1>& lhs, const NullAllocator<T2>& rhs ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two NullAllocator objects.
//
// \param lhs The left-hand side null allocator.
// \param rhs The right-hand side null allocator.
// \return \a true.
*/
template< typename T1    // Type of the left-hand side null allocator
        , typename T2 >  // Type of the right-hand side null allocator
inline bool operator==( const NullAllocator<T1>& lhs, const NullAllocator<T2>& rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two NullAllocator objects.
//
// \param lhs The left-hand side null allocator.
// \param rhs The right-hand side null allocator.
// \return \a false.
*/
template< typename T1    // Type of the left-hand side null allocator
        , typename T2 >  // Type of the right-hand side null allocator
inline bool operator!=( const NullAllocator<T1>& lhs, const NullAllocator<T2>& rhs ) noexcept
{
   MAYBE_UNUSED( lhs, rhs );
   return false;
}
//*************************************************************************************************

} // namespace blaze

#endif
