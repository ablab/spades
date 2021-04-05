//=================================================================================================
/*!
//  \file blaze/math/expressions/MatMatMultExpr.h
//  \brief Header file for the MatMatMultExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATMATMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATMATMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MultExpr.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all matrix/matrix multiplication expression templates.
// \ingroup math
//
// The MatMatMultExpr class serves as a tag for all expression templates that implement a
// matrix/matrix multiplication. All classes, that represent a matrix multiplication and
// that are used within the expression template environment of the Blaze library have to
// derive publicly from this class in order to qualify as matrix multiplication expression
// template. Only in case a class is derived publicly from the MatMatMultExpr base class,
// the IsMatMatMultExpr type trait recognizes the class as valid matrix multiplication
// expression template.
*/
template< typename MT >  // Matrix base type of the expression
struct MatMatMultExpr
   : public MultExpr<MT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-matrix multiplication
//        expression and a dense vector (\f$ \vec{y}=(A*B)*\vec{x} \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-matrix multiplication.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// matrix-matrix multiplication expression and a dense vector. It restructures the expression
// \f$ \vec{x}=(A*B)*\vec{x} \f$ to the expression \f$ \vec{y}=A*(B*\vec{x}) \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator*( const MatMatMultExpr<MT>& mat, const DenseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (*mat).leftOperand() * ( (*mat).rightOperand() * vec );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-matrix multiplication
//        expression and a sparse vector (\f$ \vec{y}=(A*B)*\vec{x} \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-matrix multiplication.
// \param vec The right-hand side sparse vector for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// matrix-matrix multiplication expression and a sparse vector. It restructures the expression
// \f$ \vec{y}=(A*B)*\vec{x} \f$ to the expression \f$ \vec{y}=A*(B*\vec{x}) \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Type of the right-hand side sparse vector
inline decltype(auto)
   operator*( const MatMatMultExpr<MT>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (*mat).leftOperand() * ( (*mat).rightOperand() * vec );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose dense vector and a
//        matrix-matrix multiplication expression (\f$ \vec{y}^T=\vec{x}^T*(A*B) \f$).
// \ingroup math
//
// \param vec The left-hand side dense vector for the multiplication.
// \param mat The right-hand side matrix-matrix multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a dense
// vector and a matrix-matrix multiplication expression. It restructures the expression
// \f$ \vec{y}^T=\vec{x}^T*(A*B) \f$ to the expression \f$ \vec{y}^T=(\vec{x}^T*A)*B \f$.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseVector<VT,true>& vec, const MatMatMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( vec * (*mat).leftOperand() ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose sparse vector and a
//        matrix-matrix multiplication expression (\f$ \vec{y}^T=\vec{x}^T*(A*B) \f$).
// \ingroup math
//
// \param vec The left-hand side sparse vector for the multiplication.
// \param mat The right-hand side matrix-matrix multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a sparse
// vector and a matrix-matrix multiplication expression. It restructures the expression
// \f$ \vec{y}^T=\vec{x}^T*(A*B) \f$ to the expression \f$ \vec{y}^T=(\vec{x}^T*A)*B \f$.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const SparseVector<VT,true>& vec, const MatMatMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( vec * (*mat).leftOperand() ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
