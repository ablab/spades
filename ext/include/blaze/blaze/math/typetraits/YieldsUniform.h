//=================================================================================================
/*!
//  \file blaze/math/typetraits/YieldsUniform.h
//  \brief Header file for the YieldsUniform type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_YIELDSUNIFORM_H_
#define _BLAZE_MATH_TYPETRAITS_YIELDSUNIFORM_H_


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
/*!\brief Compile time check for operations on vectors and matrices.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given operation \a OP yields a uniform vector or
// matrix when applied to several vectors or matrices of types \a T and \a Ts. In case the
// operation yields a uniform vector or matrix, the \a value member constant is set to \a true,
// the nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class derives from
// \a FalseType.
*/
template< typename OP, typename T, typename... Ts >
struct YieldsUniform
   : public FalseType
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsUniform type trait for const types.
// \ingroup math_type_traits
*/
template< typename OP, typename T, typename... Ts >
struct YieldsUniform< const OP, T, Ts... >
   : public YieldsUniform<OP,T,Ts...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsUniform type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename OP, typename T, typename... Ts >
struct YieldsUniform< volatile OP, T, Ts... >
   : public YieldsUniform<OP,T,Ts...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsUniform type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename OP, typename T, typename... Ts >
struct YieldsUniform< const volatile OP, T, Ts... >
   : public YieldsUniform<OP,T,Ts...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the YieldsUniform type trait.
// \ingroup math_type_traits
//
// The YieldsUniform_v variable template provides a convenient shortcut to access the nested
// \a value of the YieldsUniform class template. For instance, given the operation \a OP and
// the vector or matrix type \a T the following two statements are identical:

   \code
   constexpr bool value1 = blaze::YieldsUniform<OP,T>::value;
   constexpr bool value2 = blaze::YieldsUniform_v<OP,T>;
   \endcode
*/
template< typename OP, typename T, typename... Ts >
constexpr bool YieldsUniform_v = YieldsUniform<OP,T,Ts...>::value;
//*************************************************************************************************

} // namespace blaze

#endif
