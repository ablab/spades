//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsShrinkable.h
//  \brief Header file for the IsShrinkable type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSHRINKABLE_H_
#define _BLAZE_MATH_TYPETRAITS_ISSHRINKABLE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for shrinkable data types.
// \ingroup math_type_traits
//
// This type trait tests whether the given data type is a shrinkable data type. In case the
// data type can be shrunk (via the shrinkToFit() function), the \a value member constant is
// set to \a true, the nested type definition \a Type is \a TrueType, and the class derives from
// \a TrueType. Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class
// derives from \a FalseType. Examples:

   \code
   blaze::IsShrinkable< DynamicVector<double,false> >::value       // Evaluates to 1
   blaze::IsShrinkable< const DynamicMatrix<double,false> >::Type  // Results in TrueType
   blaze::IsShrinkable< volatile CompressedMatrix<int,true> >      // Is derived from TrueType
   blaze::IsShrinkable< int >::value                               // Evaluates to 0
   blaze::IsShrinkable< const StaticVector<float,3U,false >::Type  // Results in FalseType
   blaze::IsShrinkable< volatile HybridMatrix<int,3U,4U,false> >   // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsShrinkable
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsShrinkable type trait for const types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsShrinkable< const T >
   : public IsShrinkable<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsShrinkable type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsShrinkable< volatile T >
   : public IsShrinkable<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsShrinkable type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename T >
struct IsShrinkable< const volatile T >
   : public IsShrinkable<T>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsShrinkable type trait.
// \ingroup math_type_traits
//
// The IsShrinkable_v variable template provides a convenient shortcut to access the nested
// \a value of the IsShrinkable class template. For instance, given the type \a T the following
// two statements are identical:

   \code
   constexpr bool value1 = blaze::IsShrinkable<T>::value;
   constexpr bool value2 = blaze::IsShrinkable_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsShrinkable_v = IsShrinkable<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
