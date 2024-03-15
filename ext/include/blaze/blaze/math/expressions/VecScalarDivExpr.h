//=================================================================================================
/*!
//  \file blaze/math/expressions/VecScalarDivExpr.h
//  \brief Header file for the VecScalarDivExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECSCALARDIVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECSCALARDIVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************


#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DivExpr.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all vector/scalar division expression templates.
// \ingroup math
//
// The VecScalarDivExpr class serves as a tag for all expression templates that implement a
// vector/scalar division. All classes, that represent a vector/scalar division and that are
// used within the expression template environment of the Blaze library have to derive publicly
// from this class in order to qualify as vector/scalar division expression template. Only in
// case a class is derived publicly from the VecScalarDivExpr base class, the IsVecScalarDivExpr
// type trait recognizes the class as valid vector/scalar division expression template.
*/
template< typename VT >  // Vector base type of the expression
struct VecScalarDivExpr
   : public DivExpr<VT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar division
//        expression and a scalar value (\f$ \vec{a}=(\vec{b}/s1)*s2 \f$).
// \ingroup math
//
// \param vec The left-hand side vector-scalar division.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// vector-scalar division expression and a scalar value.
*/
template< typename VT  // Vector base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> &&
                      ( IsInvertible_v<ST> ||
                        IsInvertible_v< RightOperand_t< VectorType_t<VT> > > ) >* = nullptr >
inline decltype(auto) operator*( const VecScalarDivExpr<VT>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return (*vec).leftOperand() * ( scalar / (*vec).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a scalar value and a vector-scalar
//        division expression (\f$ \vec{a}=s2*(\vec{b}/s1) \f$).
// \ingroup math
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side vector-scalar division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a vector-scalar division expression.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Vector base type of the expression
        , EnableIf_t< IsScalar_v<ST> &&
                      ( IsInvertible_v<ST> ||
                        IsInvertible_v< RightOperand_t< VectorType_t<VT> > > ) >* = nullptr >
inline decltype(auto) operator*( ST scalar, const VecScalarDivExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (*vec).leftOperand() * ( scalar / (*vec).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a vector-scalar division expression
//        and a scalar value (\f$ \vec{a}=(\vec{b}/s1)/s2 \f$).
// \ingroup math
//
// \param vec The left-hand side vector-scalar division.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the division of a vector-scalar
// division expression and a scalar value.
*/
template< typename VT  // Vector base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator/( const VecScalarDivExpr<VT>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != ST(0), "Division by zero detected" );

   return (*vec).leftOperand() / ( (*vec).rightOperand() * scalar );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculation of the transpose of the given vector-scalar division.
// \ingroup math
//
// \param vector The vector-scalar division expression to be transposed.
// \return The transpose of the expression.
//
// This operator implements the performance optimized treatment of the transpose of a
// vector-scalar division. It restructures the expression \f$ a=trans(b/s1) \f$ to
// the expression \f$ a=trans(b)/s1 \f$.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) trans( const VecScalarDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return trans( (*vector).leftOperand() ) / (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculation of the complex conjugate of the given vector-scalar division.
// \ingroup math
//
// \param vector The vector-scalar division expression to be conjugated.
// \return The complex conjugate of the expression.
//
// This operator implements the performance optimized treatment of the complex conjugate
// of a vector-scalar division. It restructures the expression \f$ a=conj(b/s1) \f$ to the
// expression \f$ a=conj(b)/s1 \f$.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) conj( const VecScalarDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return conj( (*vector).leftOperand() ) / (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
