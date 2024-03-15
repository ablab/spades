//=================================================================================================
/*!
//  \file blaze/math/constraints/VecVecAddExpr.h
//  \brief Constraint on the data type
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

#ifndef _BLAZE_MATH_CONSTRAINTS_VECVECADDEXPR_H_
#define _BLAZE_MATH_CONSTRAINTS_VECVECADDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsVecVecAddExpr.h>
#include <blaze/math/typetraits/Size.h>


namespace blaze {

//=================================================================================================
//
//  MUST_BE_VECVECADDEXPR_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is not a vector/vector addition expression (i.e. a type
// derived from the VecVecAddExpr base class), a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_BE_VECVECADDEXPR_TYPE(T) \
   static_assert( ::blaze::IsVecVecAddExpr_v<T>, "Non-vector/vector addition expression type detected" )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_NOT_BE_VECVECADDEXPR_TYPE CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data type \a T is a vector/vector addition expression (i.e. a type derived
// from the VecVecAddExpr base class), a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_NOT_BE_VECVECADDEXPR_TYPE(T) \
   static_assert( !::blaze::IsVecVecAddExpr_v<T>, "Vector/vector addition expression type detected" )
//*************************************************************************************************




//=================================================================================================
//
//  MUST_FORM_VALID_VECVECADDEXPR CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Constraint on the data type.
// \ingroup math_constraints
//
// In case the given data types \a T1 and \a T2 do not form a valid vector/vector addition,
// a compilation error is created.
*/
#define BLAZE_CONSTRAINT_MUST_FORM_VALID_VECVECADDEXPR(T1,T2) \
   static_assert( ( ( ::blaze::IsRowVector_v<T1> && ::blaze::IsRowVector_v<T2> ) || \
                    ( ::blaze::IsColumnVector_v<T1> && ::blaze::IsColumnVector_v<T2> ) ) && \
                  ( ( ::blaze::Size_v<T1,0UL> == -1L ) || \
                    ( ::blaze::Size_v<T2,0UL> == -1L ) || \
                    ( ::blaze::Size_v<T1,0UL> == ::blaze::Size_v<T2,0UL> ) ) \
                , "Invalid vector/vector addition expression detected" )
//*************************************************************************************************

} // namespace blaze

#endif
