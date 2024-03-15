//=================================================================================================
/*!
//  \file blaze/util/smallarray/SmallArrayData.h
//  \brief Header file for the SmallArrayData class template
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

#ifndef _BLAZE_UTIL_SMALLARRAY_SMALLARRAYDATA_H_
#define _BLAZE_UTIL_SMALLARRAY_SMALLARRAYDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>
#include <blaze/util/typetraits/AlignmentOf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class for the SmallArray class template.
// \ingroup util
//
// This auxiliary class template is used as backend by the SmallArray class template. It is
// responsible to provide a static storage of size \a N.
*/
template< typename T  // Data type of the elements
        , size_t N >  // Number of preallocated elements
struct SmallArrayData
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr SmallArrayData() noexcept;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr T*       array()       noexcept;
   constexpr const T* array() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   alignas( AlignmentOf_v<T> ) byte_t v_[N*sizeof(T)];  //!< The static storage.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default constructor for SmallArrayData. No element initialization is performed!
//
// \note This constructor does not perform any kind of element initialization!
*/
template< typename T  // Data type of the elements
        , size_t N >  // Number of preallocated elements
constexpr SmallArrayData<T,N>::SmallArrayData() noexcept
   // v_ is intentionally not initialized
{}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a pointer to the first element of the static array.
//
// \return A pointer to the first element of the static array.
*/
template< typename T  // Data type of the elements
        , size_t N >  // Number of preallocated elements
constexpr T* SmallArrayData<T,N>::array() noexcept
{
   return reinterpret_cast<T*>( v_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a pointer to the first element of the static array.
//
// \return A pointer to the first element of the static array.
*/
template< typename T  // Data type of the elements
        , size_t N >  // Number of preallocated elements
constexpr const T* SmallArrayData<T,N>::array() const noexcept
{
   return reinterpret_cast<const T*>( v_ );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR N = 0
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the SmallArrayData class template for N = 0.
// \ingroup util
//
// This specialization of SmallArrayData handles the request for zero static elements of type T.
*/
template< typename T >  // Data type of the elements
struct SmallArrayData<T,0UL>
{
 public:
   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   constexpr T*       array()       noexcept;
   constexpr const T* array() const noexcept;
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a pointer to the first element of the static array.
//
// \return A pointer to the first element of the static array.
*/
template< typename T >  // Data type of the elements
constexpr T* SmallArrayData<T,0UL>::array() noexcept
{
   return nullptr;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a pointer to the first element of the static array.
//
// \return A pointer to the first element of the static array.
*/
template< typename T >  // Data type of the elements
constexpr const T* SmallArrayData<T,0UL>::array() const noexcept
{
   return nullptr;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
