//=================================================================================================
/*!
//  \file blaze/math/typetraits/YieldsTriangular.h
//  \brief Header file for the YieldsTriangular type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_YIELDSTRIANGULAR_H_
#define _BLAZE_MATH_TYPETRAITS_YIELDSTRIANGULAR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsUpper.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Compile time check for operations on matrices.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given operation \a OP yields a lower or upper
// triangular matrix when applied to several matrices of types \a MT and \a MTs. In case the
// operation yields a triangular matrix, the \a value member constant is set to \a true, the
// nested type definition \a Type is \a TrueType, and the class derives from \a TrueType.
// Otherwise \a value is set to \a false, \a Type is \a FalseType, and the class derives from
// \a FalseType.
*/
template< typename OP, typename MT, typename... MTs >
struct YieldsTriangular
   : public BoolConstant< YieldsLower_v<OP,MT,MTs...> || YieldsUpper_v<OP,MT,MTs...> >
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsTriangular type trait for const types.
// \ingroup math_type_traits
*/
template< typename OP, typename MT, typename... MTs >
struct YieldsTriangular< const OP, MT, MTs... >
   : public YieldsTriangular<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsTriangular type trait for volatile types.
// \ingroup math_type_traits
*/
template< typename OP, typename MT, typename... MTs >
struct YieldsTriangular< volatile OP, MT, MTs... >
   : public YieldsTriangular<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the YieldsTriangular type trait for cv qualified types.
// \ingroup math_type_traits
*/
template< typename OP, typename MT, typename... MTs >
struct YieldsTriangular< const volatile OP, MT, MTs... >
   : public YieldsTriangular<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the YieldsTriangular type trait.
// \ingroup math_type_traits
//
// The YieldsTriangular_v variable template provides a convenient shortcut to access the nested
// \a value of the YieldsTriangular class template. For instance, given the operation \a OP and
// the matrix type \a MT the following two statements are identical:

   \code
   constexpr bool value1 = blaze::YieldsTriangular<OP,MT>::value;
   constexpr bool value2 = blaze::YieldsTriangular_v<OP,MT>;
   \endcode
*/
template< typename OP, typename MT, typename... MTs >
constexpr bool YieldsTriangular_v = YieldsTriangular<OP,MT,MTs...>::value;
//*************************************************************************************************

} // namespace blaze

#endif
