//=================================================================================================
/*!
//  \file blaze/math/expressions/MatNoSIMDExpr.h
//  \brief Header file for the MatNoSIMDExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATNOSIMDEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATNOSIMDEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/DeclDiagExpr.h>
#include <blaze/math/expressions/DeclHermExpr.h>
#include <blaze/math/expressions/DeclLowExpr.h>
#include <blaze/math/expressions/DeclStrLowExpr.h>
#include <blaze/math/expressions/DeclStrUppExpr.h>
#include <blaze/math/expressions/DeclSymExpr.h>
#include <blaze/math/expressions/DeclUniLowExpr.h>
#include <blaze/math/expressions/DeclUniUppExpr.h>
#include <blaze/math/expressions/DeclUppExpr.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatKronExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/MatNoAliasExpr.h>
#include <blaze/math/expressions/MatRepeatExpr.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/NoSIMDExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/VecExpandExpr.h>
#include <blaze/math/expressions/VecTVecMapExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all matrix no-SIMD expression templates.
// \ingroup math
//
// The MatNoSIMDExpr class serves as a tag for all expression templates that implement a matrix
// no-SIMD operation. All classes, that represent a matrix no-SIMD operation and that are used
// within the expression template environment of the Blaze library have to derive publicly from
// this class in order to qualify as matrix no-SIMD expression template. Only in case a class is
// derived publicly from the MatNoSIMDExpr base class, the IsMatNoSIMDExpr type trait recognizes
// the class as valid matrix no-SIMD expression template.
*/
template< typename MT >  // Matrix base type of the expression
struct MatNoSIMDExpr
   : public NoSIMDExpr<MT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/matrix addition.
// ingroup math
//
// \param matrix The constant matrix/matrix addition.
// \return The SIMD-disabled addition.
//
// This function returns an expression representing the SIMD-disabled matrix/matrix addition.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMatAddExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) + nosimd( (*matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/matrix subtraction.
// ingroup math
//
// \param matrix The constant matrix/matrix subtraction.
// \return The SIMD-disabled subtraction.
//
// This function returns an expression representing the SIMD-disabled matrix/matrix subtraction.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMatSubExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) - nosimd( (*matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given Schur product.
// ingroup math
//
// \param matrix The constant Schur product.
// \return The SIMD-disabled Schur product.
//
// This function returns an expression representing the SIMD-disabled Schur product.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const SchurExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) % nosimd( (*matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/matrix multiplication.
// ingroup math
//
// \param matrix The constant matrix/matrix multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled matrix/matrix multiplication.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMatMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) * nosimd( (*matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/matrix Kronecker product.
// ingroup math
//
// \param matrix The constant matrix/matrix Kronecker product.
// \return The SIMD-disabled matrix/matrix Kronecker product.
//
// This function returns an expression representing the SIMD-disabled matrix/matrix Kronecker
// product.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMatKronExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return kron( nosimd( (*matrix).leftOperand() ), nosimd( (*matrix).rightOperand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given outer product.
// ingroup math
//
// \param matrix The constant outer product.
// \return The SIMD-disabled outer product.
//
// This function returns an expression representing the SIMD-disabled outer product.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const VecTVecMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) * nosimd( (*matrix).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/scalar multiplication.
// ingroup math
//
// \param matrix The constant matrix/scalar multiplication.
// \return The SIMD-disabled multiplication.
//
// This function returns an expression representing the SIMD-disabled matrix/scalar multiplication.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix/scalar division.
// ingroup math
//
// \param matrix The constant matrix/scalar division.
// \return The SIMD-disabled division.
//
// This function returns an expression representing the SIMD-disabled matrix/scalar division.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatScalarDivExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( (*matrix).leftOperand() ) / (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given unary matrix map operation.
// ingroup math
//
// \param matrix The constant unary matrix map operation.
// \return The SIMD-disabled unary map operation.
//
// This function returns an expression representing the SIMD-disabled unary matrix map operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( nosimd( (*matrix).operand() ), (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given binary matrix map operation.
// ingroup math
//
// \param matrix The constant binary matrix map operation.
// \return The SIMD-disabled binary map operation.
//
// This function returns an expression representing the SIMD-disabled binary matrix map operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatMatMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( nosimd( (*matrix).leftOperand() ), nosimd( (*matrix).rightOperand() ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given outer map operation.
// ingroup math
//
// \param matrix The constant outer map operation.
// \return The SIMD-disabled outer map operation.
//
// This function returns an expression representing the SIMD-disabled outer map operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const VecTVecMapExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return map( nosimd( (*matrix).leftOperand() ), nosimd( (*matrix).rightOperand() ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix evaluation operation.
// ingroup math
//
// \param matrix The constant matrix evaluation operation.
// \return The SIMD-disabled evaluation operation.
//
// This function returns an expression representing the SIMD-disabled matrix evaluation operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatEvalExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return eval( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix serialization operation.
// ingroup math
//
// \param matrix The constant matrix serialization operation.
// \return The SIMD-disabled serialization operation.
//
// This function returns an expression representing the SIMD-disabled matrix serialization
// operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatSerialExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return serial( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix no-alias operation.
// ingroup math
//
// \param matrix The constant matrix no-alias operation.
// \return The SIMD-disabled no-alias operation.
//
// This function returns an expression representing the SIMD-disabled matrix no-alias operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatNoAliasExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix no-SIMD operation.
// ingroup math
//
// \param matrix The constant matrix no-SIMD operation.
// \return The SIMD-disabled no-SIMD operation.
//
// This function returns an expression representing the SIMD-disabled matrix no-SIMD operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatNoSIMDExpr<MT>& matrix )
{
   return *matrix;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declsym operation.
// ingroup math
//
// \param matrix The constant matrix declsym operation.
// \return The SIMD-disabled declsym operation.
//
// This function returns an expression representing the SIMD-disabled matrix declsym operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclSymExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declsym( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declherm operation.
// ingroup math
//
// \param matrix The constant matrix declherm operation.
// \return The SIMD-disabled declherm operation.
//
// This function returns an expression representing the SIMD-disabled matrix declherm operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclHermExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declherm( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix decllow operation.
// ingroup math
//
// \param matrix The constant matrix decllow operation.
// \return The SIMD-disabled decllow operation.
//
// This function returns an expression representing the SIMD-disabled matrix decllow operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclLowExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return decllow( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declunilow operation.
// ingroup math
//
// \param matrix The constant matrix declunilow operation.
// \return The SIMD-disabled declunilow operation.
//
// This function returns an expression representing the SIMD-disabled matrix declunilow operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclUniLowExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declunilow( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declstrlow operation.
// ingroup math
//
// \param matrix The constant matrix declstrlow operation.
// \return The SIMD-disabled declstrlow operation.
//
// This function returns an expression representing the SIMD-disabled matrix declstrlow operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclStrLowExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declstrlow( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declupp operation.
// ingroup math
//
// \param matrix The constant matrix declupp operation.
// \return The SIMD-disabled declupp operation.
//
// This function returns an expression representing the SIMD-disabled matrix declupp operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclUppExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declupp( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix decluniupp operation.
// ingroup math
//
// \param matrix The constant matrix decluniupp operation.
// \return The SIMD-disabled decluniupp operation.
//
// This function returns an expression representing the SIMD-disabled matrix decluniupp operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclUniUppExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return decluniupp( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix declstrupp operation.
// ingroup math
//
// \param matrix The constant matrix declstrupp operation.
// \return The SIMD-disabled declstrupp operation.
//
// This function returns an expression representing the SIMD-disabled matrix declstrupp operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclStrUppExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return declstrupp( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix decldiag operation.
// ingroup math
//
// \param matrix The constant matrix decldiag operation.
// \return The SIMD-disabled decldiag operation.
//
// This function returns an expression representing the SIMD-disabled matrix decldiag operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const DeclDiagExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return decldiag( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix transpose operation.
// ingroup math
//
// \param matrix The constant matrix transpose operation.
// \return The SIMD-disabled transpose operation.
//
// This function returns an expression representing the SIMD-disabled matrix transpose operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatTransExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return trans( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector expansion operation.
// ingroup math
//
// \param matrix The constant vector expansion operation.
// \return The SIMD-disabled expansion operation.
//
// This function returns an expression representing the SIMD-disabled vector expansion operation.
*/
template< typename MT  // Matrix base type of the expression
        , size_t E >   // Compile time expansion argument
inline decltype(auto) nosimd( const VecExpandExpr<MT,E>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return expand<E>( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given vector expansion operation.
// ingroup math
//
// \param matrix The constant vector expansion operation.
// \return The SIMD-disabled expansion operation.
//
// This function returns an expression representing the SIMD-disabled vector expansion operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const VecExpandExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return expand( nosimd( (*matrix).operand() ), (*matrix).expansion() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix repeat operation.
// ingroup math
//
// \param matrix The constant matrix repeat operation.
// \return The SIMD-disabled repeat operation.
//
// This function returns an expression representing the SIMD-disabled matrix repeat operation.
*/
template< typename MT  // Matrix base type of the expression
        , size_t R0    // Compile time row-wise repetitions
        , size_t R1 >  // Compile time column-wise repetitions
inline decltype(auto) nosimd( const MatRepeatExpr<MT,R0,R1>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return repeat<R0,R1>( nosimd( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Disable the SIMD evaluation of the given matrix repeat operation.
// ingroup math
//
// \param matrix The constant matrix repeat operation.
// \return The SIMD-disabled repeat operation.
//
// This function returns an expression representing the SIMD-disabled matrix repeat operation.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) nosimd( const MatRepeatExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return repeat( nosimd( (*matrix).operand() )
                , (*matrix).template repetitions<0UL>()
                , (*matrix).template repetitions<1UL>() );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
