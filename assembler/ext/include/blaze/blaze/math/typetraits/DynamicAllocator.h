//=================================================================================================
/*!
//  \file blaze/math/typetraits/DynamicAllocator.h
//  \brief Header file for the DynamicAllocator type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_DYNAMICALLOCATOR_H_
#define _BLAZE_MATH_TYPETRAITS_DYNAMICALLOCATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/util/AlignedAllocator.h>
#include <blaze/util/NullAllocator.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Deduction of an allocator type for dynamic vectors and matrices.
// \ingroup math_type_traits
//
// \section dynamicallocator_general General
//
// The DynamicAllocator type trait deduces the allocator type for dynamic vectors and matrices.
// Given one or two allocators, it provides a nested \a Type alias template, which results in
// the according allocator type:

   \code
   using A1 = AlignedAllocator<int>;
   using A2 = AlignedAllocator<double>;

   blaze::DynamicAllocator<A1,A2>::Type<double>  // Results in 'AlignedAllocator<double>'
   \endcode

// In case no resulting allocator type can be determined, the nested \a Type template will result
// in \a blaze::NullAllocator for all possible types.
//
//
// \n \section dynamicallocator_specializations Creating custom specializations
//
// DynamicAllocator is guaranteed to work only for all \b Blaze allocator types. In order to add
// support for user-defined allocator types it is possible to specialize the DynamicAllocator
// template. The following example demonstrates the according specialization for the
// \a blaze::AlignedAllocator class template:

   \code
   template< typename T1, typename T2 >
   struct DynamicAllocator< AlignedAllocator<T1>, AlignedAllocator<T2> >
   {
      template< typename U >
      using Type = AlignedAllocator<U>;
   };
   \endcode
*/
template< typename A1, typename... As >
struct DynamicAllocator
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename U >
   using Type = NullAllocator<U>;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for 'AlignedAllocator'.
// \ingroup math_type_traits
*/
template< typename T >
struct DynamicAllocator< AlignedAllocator<T> >
{
   template< typename U >
   using Type = AlignedAllocator<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for 'NullAllocator'.
// \ingroup math_type_traits
*/
template< typename T >
struct DynamicAllocator< NullAllocator<T> >
{
   template< typename U >
   using Type = AlignedAllocator<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for two 'AlignedAllocator'.
// \ingroup math_type_traits
*/
template< typename T1, typename T2 >
struct DynamicAllocator< AlignedAllocator<T1>, AlignedAllocator<T2> >
{
   template< typename U >
   using Type = AlignedAllocator<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for any allocator and 'NullAllocator'.
// \ingroup math_type_traits
*/
template< typename A1, typename T2 >
struct DynamicAllocator< A1, NullAllocator<T2> >
{
   template< typename U >
   using Type = typename std::allocator_traits<A1>::template rebind_alloc<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for 'NullAllocator' and any allocator.
// \ingroup math_type_traits
*/
template< typename A2, typename T1 >
struct DynamicAllocator< NullAllocator<T1>, A2 >
{
   template< typename U >
   using Type = typename std::allocator_traits<A2>::template rebind_alloc<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DynamicAllocator type trait for two 'NullAllocator'.
// \ingroup math_type_traits
*/
template< typename T1, typename T2 >
struct DynamicAllocator< NullAllocator<T1>, NullAllocator<T2> >
{
   template< typename U >
   using Type = AlignedAllocator<U>;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DynamicAllocator type trait.
// \ingroup math_type_traits
//
// The DynamicAllocator_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the DynamicAllocator class template. For instance, given the types \a A1 and \a A2
// the following  two type definitions are identical:

   \code
   using Type1 = typename blaze::DynamicAllocator<A1,A2>::template Type<T>;
   using Type2 = blaze::DynamicAllocator_t<T,A1,A2>;
   \endcode
*/
template< typename T, typename... As >
using DynamicAllocator_t = typename DynamicAllocator<As...>::template Type<T>;
//*************************************************************************************************

} // namespace blaze

#endif
