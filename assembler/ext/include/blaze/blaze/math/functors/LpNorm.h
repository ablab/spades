//=================================================================================================
/*!
//  \file blaze/math/functors/LpNorm.h
//  \brief Header file for the LpNorm functor
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

#ifndef _BLAZE_MATH_FUNCTORS_LPNORM_H_
#define _BLAZE_MATH_FUNCTORS_LPNORM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/util/StaticAssert.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for the lpNorm() function.
// \ingroup functors
*/
template< size_t... P >  // Compile time norm parameter
struct LpNorm
{
   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \return The Lp norm of the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( T&& a ) const
   {
      BLAZE_STATIC_ASSERT_MSG( sizeof...( P ) == 1UL, "Missing norm parameter detected" );
      return lpNorm( std::forward<T>( a ), P... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \param p The runtime norm parameter.
   // \return The Lp norm of the given object/value.
   */
   template< typename T, typename ST >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( T&& a, ST p ) const
   {
      BLAZE_STATIC_ASSERT_MSG( sizeof...( P ) == 0UL, "Over-specified norm parameter detected" );
      return lpNorm( std::forward<T>( a ), p );
   }
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CLASS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of LpNorm class template for the L1 norm.
// \ingroup functors
*/
template<>
struct LpNorm<1UL>
{
   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \return The Lp norm of the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( const T& a ) const
   {
      return l1Norm( a );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of LpNorm class template for the L2 norm.
// \ingroup functors
*/
template<>
struct LpNorm<2UL>
{
   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \return The Lp norm of the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( const T& a ) const
   {
      return l2Norm( a );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of LpNorm class template for the L3 norm.
// \ingroup functors
*/
template<>
struct LpNorm<3UL>
{
   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \return The Lp norm of the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( const T& a ) const
   {
      return l3Norm( a );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of LpNorm class template for the L4 norm.
// \ingroup functors
*/
template<>
struct LpNorm<4UL>
{
   //**********************************************************************************************
   /*!\brief Calls the lpNorm() function with the given object/value.
   //
   // \param a The given object/value.
   // \return The Lp norm of the given object/value.
   */
   template< typename T >
   BLAZE_ALWAYS_INLINE BLAZE_DEVICE_CALLABLE decltype(auto) operator()( const T& a ) const
   {
      return l4Norm( a );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
