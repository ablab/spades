//=================================================================================================
/*!
//  \file blaze/math/expressions/VecNoSIMDExpr.h
//  \brief Header file for the VecNoSIMDExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECNOSIMDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECNOSIMDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/CrossExpr.h>
#include <blaze/math/expressions/MatReduceExpr.h>
#include <blaze/math/expressions/MatVecMultExpr.h>
#include <blaze/math/expressions/NoSIMDExpr.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/expressions/VecEvalExpr.h>
#include <blaze/math/expressions/VecMapExpr.h>
#include <blaze/math/expressions/VecNoAliasExpr.h>
#include <blaze/math/expressions/VecRepeatExpr.h>
#include <blaze/math/expressions/VecScalarMultExpr.h>
#include <blaze/math/expressions/VecScalarDivExpr.h>
#include <blaze/math/expressions/VecSerialExpr.h>
#include <blaze/math/expressions/VecTransExpr.h>
#include <blaze/math/expressions/VecVecAddExpr.h>
#include <blaze/math/expressions/VecVecDivExpr.h>
#include <blaze/math/expressions/VecVecKronExpr.h>
#include <blaze/math/expressions/VecVecMapExpr.h>
#include <blaze/math/expressions/VecVecMultExpr.h>
#include <blaze/math/expressions/VecVecSubExpr.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all vector no-SIMD expression templates.
// \ingroup math
//
// The VecNoSIMDExpr class serves as a tag for all expression templates that implement a vector
// no-SIMD operation. All classes, that represent a vector no-SIMD operation and that are used
// within the expression template environment of the Blaze library have to derive publicly from
// this class in order to qualify as vector no-SIMD expression template. Only in case a class is
// derived publicly from the VecNoSIMDExpr base class, the IsVecNoSIMDExpr type trait recognizes
// the class as valid vector no-SIMD expression template.
*/
template< typename VT >  // Vector base type of the expression
struct VecNoSIMDExpr
   : public NoSIMDExpr<VT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector addition.
// \ingroup math
//
// \param vector The constant vector/vector addition.
// \return The SIMD-disabled addition.
//
// This function returns an expression representing the SIMD-disabled vector/vector addition.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecAddExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) + nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector subtraction.
// \ingroup math
//
// \param vector The constant vector/vector subtraction.
// \return The SIMD-disabled subtraction.
//
// This function returns an expression representing the SIMD-disabled vector/vector subtraction.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecSubExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) - nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector multiplication.
// \ingroup math
//
// \param vector The constant vector/vector multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled vector/vector multiplication.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) * nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector Kronecker product.
// \ingroup math
//
// \param vector The constant vector/vector Kronecker product.
// \return The SIMD-disabled Kronecker product.
//
// This function returns an expression representing the SIMD-disabled vector/vector Kronecker
// product.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecKronExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return kron( nosimd( (*vector).leftOperand() ), nosimd( (*vector).rightOperand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector division.
// \ingroup math
//
// \param vector The constant vector/vector division.
// \return The SIMD-disabled division.
//
// This function returns an expression representing the SIMD-disabled vector/vector division.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) / nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/vector cross product.
// \ingroup math
//
// \param vector The constant vector/vector cross product.
// \return The SIMD-disabled cross product.
//
// This function returns an expression representing the SIMD-disabled vector/vector cross product.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const CrossExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) % nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/scalar multiplication.
// \ingroup math
//
// \param vector The constant vector/scalar multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled vector/scalar multiplication.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecScalarMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) * (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/scalar division.
// \ingroup math
//
// \param vector The constant vector/scalar division.
// \return The SIMD-disabled division.
//
// This function returns an expression representing the SIMD-disabled vector/scalar division.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecScalarDivExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) / (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given unary vector map operation.
// \ingroup math
//
// \param vector The constant unary vector map operation.
// \return The SIMD-disabled unary map operation.
//
// This function returns an expression representing the SIMD-disabled unary vector map operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecMapExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return map( nosimd( (*vector).operand() ), (*vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given binary vector map operation.
// \ingroup math
//
// \param vector The constant binary vector map operation.
// \return The SIMD-disabled binary map operation.
//
// This function returns an expression representing the SIMD-disabled binary vector map operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecVecMapExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return map( nosimd( (*vector).leftOperand() ), nosimd( (*vector).rightOperand() ),
               (*vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector evaluation operation.
// \ingroup math
//
// \param vector The constant vector evaluation operation.
// \return The SIMD-disabled evaluation operation.
//
// This function returns an expression representing the SIMD-disabled vector evaluation operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecEvalExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return eval( nosimd( (*vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector serialization operation.
// \ingroup math
//
// \param vector The constant vector serialization operation.
// \return The SIMD-disabled serialization operation.
//
// This function returns an expression representing the SIMD-disabled vector serialization
// operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecSerialExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return serial( nosimd( (*vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector no-alias operation.
// \ingroup math
//
// \param vector The constant vector no-alias operation.
// \return The SIMD-disabled no-alias operation.
//
// This function returns an expression representing the SIMD-disabled vector no-alias operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecNoAliasExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( nosimd( (*vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector no-SIMD operation.
// \ingroup math
//
// \param vector The constant vector no-SIMD operation.
// \return The SIMD-disabled no-SIMD operation.
//
// This function returns an expression representing the SIMD-disabled vector no-SIMD operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecNoSIMDExpr<VT>& vector )
{
   return *vector;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector transpose operation.
// \ingroup math
//
// \param vector The constant vector transpose operation.
// \return The SIMD-disabled transpose operation.
//
// This function returns an expression representing the SIMD-disabled vector transpose
// operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecTransExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return trans( nosimd( (*vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/vector multiplication.
// \ingroup math
//
// \param vector The constant matrix/vector multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled matrix/vector
// multiplication.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const MatVecMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) * nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector/matrix multiplication.
// \ingroup math
//
// \param vector The constant vector/matrix multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled vector/matrix
// multiplication.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const TVecMatMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*vector).leftOperand() ) * nosimd( (*vector).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix reduction operation.
// \ingroup math
//
// \param vector The constant matrix reduction operation.
// \return The SIMD-disabled reduction operation.
//
// This function returns an expression representing the SIMD-disabled matrix reduction
// operation.
*/
template< typename VT         // Vector base type of the expression
        , ReductionFlag RF >  // Reduction flag
inline decltype(auto) nosimd( const MatReduceExpr<VT,RF>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return reduce<RF>( nosimd( (*vector).operand() ), (*vector).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector repeat operation.
// ingroup math
//
// \param vector The constant vector repeat operation.
// \return The SIMD-disabled repeat operation.
//
// This function returns an expression representing the SIMD-disabled vector repeat operation.
*/
template< typename VT  // Vector base type of the expression
        , size_t R0 >  // Compile time repetitions
inline decltype(auto) nosimd( const VecRepeatExpr<VT,R0>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return repeat<R0>( nosimd( (*vector).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector repeat operation.
// ingroup math
//
// \param vector The constant vector repeat operation.
// \return The SIMD-disabled repeat operation.
//
// This function returns an expression representing the SIMD-disabled vector repeat operation.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) nosimd( const VecRepeatExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return repeat( nosimd( (*vector).operand() ), (*vector).template repetitions<0UL>() );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
