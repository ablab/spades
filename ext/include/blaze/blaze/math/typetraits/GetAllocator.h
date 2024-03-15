//=================================================================================================
/*!
//  \file blaze/math/typetraits/GetAllocator.h
//  \brief Header file for the GetAllocator type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_GETALLOCATOR_H_
#define _BLAZE_MATH_TYPETRAITS_GETALLOCATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/typetraits/UnderlyingElement.h>
#include <blaze/util/NullAllocator.h>
#include <blaze/util/typetraits/IsDetected.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Determines the type of allocator of the given type.
// \ingroup math_type_traits
//
// This type trait determines the type of allocator of the given type \a T. In case the type
// \a T exposes this allocator via a nested type \a AllocatorType, the nested type \a Type is
// set to the type of allocator. Otherwise the nested type \a Type is set to NullAllocator.

   \code
   using blaze::DynamicVector;
   using blaze::DynamicMatrix;
   using blaze::StaticVector;
   using blaze::CompressedMatrix;

   blaze::GetAllocator< DynamicVector<int> >::Type              // Results in 'AlignedAllocator<int>'
   blaze::GetAllocator< const DynamicVector<double> >::Type     // Results in 'AlignedAllocator<double>'
   blaze::GetAllocator< volatile DynamicMatrix<int> >::Type     // Results in 'AlignedAllocator<int>'
   blaze::GetAllocator< int >::Type                             // Results in 'NullAllocator<int>'
   blaze::GetAllocator< const StaticVector<float,3UL> >::Type   // Results in 'NullAllocator<float>'
   blaze::GetAllocator< volatile CompressedMatrix<int> >::Type  // Results in 'NullAllocator<int>'
   \endcode
*/
template< typename T >
using GetAllocator = DetectedOr< NullAllocator< UnderlyingElement_t<T> >, AllocatorType_t, T >;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the GetAllocator type trait.
// \ingroup math_type_traits
//
// The GetAllocator_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the GetAllocator class template. For instance, given the type \a T the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::GetAllocator<T>::Type;
   using Type2 = blaze::GetAllocator_t<T>;
   \endcode
*/
template< typename T >
using GetAllocator_t = typename GetAllocator<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
