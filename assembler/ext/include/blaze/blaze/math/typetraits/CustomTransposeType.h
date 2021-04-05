//=================================================================================================
/*!
//  \file blaze/math/typetraits/CustomTransposeType.h
//  \brief Header file for the CustomTransposeType type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_CUSTOMTRANSPOSETYPE_H_
#define _BLAZE_MATH_TYPETRAITS_CUSTOMTRANSPOSETYPE_H_


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the custom transpose type.
// \ingroup math_type_traits
//
// The CustomTransposeType type trait provides a compile time indirection for the evaluation of
// the transpose type of a custom vector or custom matrix. By default, the nested type alias
// \c Type is set to the \a TransposeType of the given vector or matrix \a T, but it is possible
// to specialize CustomTransposeType. Note that cv-qualifiers are not explicitly handled.
*/
template< typename T >
struct CustomTransposeType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename T::TransposeType;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the CustomTransposeType type trait.
// \ingroup math_type_traits
//
// The CustomTransposeType_t alias declaration provides a convenient shortcut to access the
// nested \a Type of the CustomTransposeType class template. For instance, given the type \a T
// the following two type definitions are identical:

   \code
   using Type1 = typename blaze::CustomTransposeType<T>::Type;
   using Type2 = blaze::CustomTransposeType_t<T>;
   \endcode
*/
template< typename T >
using CustomTransposeType_t = typename CustomTransposeType<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
