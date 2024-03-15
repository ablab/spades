//=================================================================================================
/*!
//  \file blaze/math/simd/SIMDPack.h
//  \brief Header file for the SIMDPack base class
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

#ifndef _BLAZE_MATH_SIMD_SIMDPACK_H_
#define _BLAZE_MATH_SIMD_SIMDPACK_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all SIMD data types.
// \ingroup simd
//
// The SIMDPack class template is a base class for all SIMD data types within the Blaze library.
// It provides an abstraction from the actual type of the SIMD pack, but enables a conversion
// back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
*/
template< typename T >  // Type of the SIMD pack
class SIMDPack
{
 public:
   //**Conversion operators************************************************************************
   /*!\name Conversion operators */
   //@{
   [[deprecated]] BLAZE_ALWAYS_INLINE constexpr T&       operator~()       noexcept;
   [[deprecated]] BLAZE_ALWAYS_INLINE constexpr const T& operator~() const noexcept;

   constexpr T&       operator*()       noexcept;
   constexpr const T& operator*() const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   SIMDPack() = default;
   SIMDPack( const SIMDPack& ) = default;
   SIMDPack( SIMDPack&& ) = default;
   ~SIMDPack() = default;
   SIMDPack& operator=( const SIMDPack& ) = default;
   SIMDPack& operator=( SIMDPack&& ) = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion operator for non-constant SIMD packs.
//
// \return Mutable reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a T of the
// SIMD pack. It will return a mutable reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
[[deprecated]] BLAZE_ALWAYS_INLINE constexpr T& SIMDPack<T>::operator~() noexcept
{
   return static_cast<T&>( *this );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Conversion operator for constant SIMD packs.
//
// \return Constant reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a T of the
// SIMD pack. It will return a constant reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
[[deprecated]] BLAZE_ALWAYS_INLINE constexpr const T& SIMDPack<T>::operator~() const noexcept
{
   return static_cast<const T&>( *this );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Conversion operator for non-constant SIMD packs.
//
// \return Mutable reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based down-cast to the actual type \a T of the SIMD pack.
// It will return a mutable reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
BLAZE_ALWAYS_INLINE constexpr T& SIMDPack<T>::operator*() noexcept
{
   return static_cast<T&>( *this );
}
//*************************************************************************************************

//*************************************************************************************************
/*!\brief Conversion operator for constant SIMD packs.
//
// \return Constant reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based down-cast to the actual type \a T of the SIMD pack.
// It will return a constant reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
BLAZE_ALWAYS_INLINE constexpr const T& SIMDPack<T>::operator*() const noexcept
{
   return static_cast<const T&>( *this );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SIMDPack global functions */
//@{
template< typename T >
T& crtp_cast( SIMDPack<T>& pack );

template< typename T >
const T& crtp_cast( const SIMDPack<T>& pack );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for non-constant SIMD packs.
//
// \param pack The SIMD pack to be downcast.
// \return Mutable reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a T of the
// SIMD pack. It will return a mutable reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
BLAZE_ALWAYS_INLINE T& crtp_cast( SIMDPack<T>& pack )
{
   return *pack;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief CRTP-based conversion operation for constant SIMD packs.
//
// \param pack The SIMD pack to be downcast.
// \return Const reference of the actual type of the SIMD pack.
//
// This operator performs the CRTP-based type-safe downcast to the actual type \a T of the
// SIMD pack. It will return a constant reference to the actual type \a T.
*/
template< typename T >  // Type of the SIMD pack
BLAZE_ALWAYS_INLINE const T& crtp_cast( const SIMDPack<T>& pack )
{
   return *pack;
}
//*************************************************************************************************

} // namespace blaze

#endif
