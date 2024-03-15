//=================================================================================================
/*!
//  \file blaze/math/traits/EvaluateTrait.h
//  \brief Header file for the EvaluateTrait class template
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

#ifndef _BLAZE_MATH_TRAITS_EVALUATETRAIT_H_
#define _BLAZE_MATH_TRAITS_EVALUATETRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/typetraits/RemoveCVRef.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the return type of the evaluate function.
// \ingroup math_traits
//
// Via this type trait it is possible to evaluate the resulting return type of the evaluate
// function. Given the data type \a T, the nested type \a Type corresponds to the resulting
// return type. In case \a T is a dense or sparse vector or matrix type, the resulting \a Type
// is set to the according result type. Otherwise, \a Type is set to the unqualified and
// unmodified \a T.
*/
template< typename T         // Type of the operand
        , typename = void >  // Restricting condition
struct EvaluateTrait
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = RemoveCVRef_t<T>;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of EvaluateTrait for vector and matrix types.
// \ingroup math_traits
*/
template< typename T >  // Type of the operand
struct EvaluateTrait< T, EnableIf_t< IsVector_v<T> || IsMatrix_v<T> > >
{
   //**********************************************************************************************
   using Type = typename T::ResultType;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the EvaluateTrait type trait.
// \ingroup math_traits
//
// The EvaluateTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the EvaluateTrait class template. For instance, given the type \a T the following
// two type definitions are identical:

   \code
   using Type1 = typename blaze::EvaluateTrait<T>::Type;
   using Type2 = blaze::EvaluateTrait_t<T>;
   \endcode
*/
template< typename T >
using EvaluateTrait_t = typename EvaluateTrait<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
