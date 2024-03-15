//=================================================================================================
/*!
//  \file blaze/math/expressions/MatInvExpr.h
//  \brief Header file for the MatInvExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATINVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATINVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/Expression.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/typetraits/IsMatScalarMultExpr.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Base class for all matrix inversion expression templates.
// \ingroup math
//
// The MatInvExpr class serves as a tag for all expression templates that implement a matrix
// inversion operation. All classes, that represent a matrix inversion operation and that are
// used within the expression template environment of the Blaze library have to derive publicly
// from this class in order to qualify as matrix inversion expression template. Only in case
// a class is derived publicly from the MatInvExpr base class, the IsMatInvExpr type trait
// recognizes the class as valid matrix inversion expression template.
*/
template< typename MT >  // Matrix base type of the expression
struct MatInvExpr
   : public Expression<MT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix inversion expression
//        and a dense vector (\f$ \vec{y}=inv(A)*\vec{x} \f$).
// \ingroup math
//
// \param mat The left-hand side matrix inversion.
// \param vec The right-hand side dense vector for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// matrix inversion expression and a dense vector. It restructures the expression
// \f$ \vec{x}=inv(A)*\vec{x} \f$ to the expression \f$ \vec{y}=solve(A,\vec{x}) \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator*( const MatInvExpr<MT>& mat, const DenseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return solve( (*mat).operand(), *vec );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose dense vector and a
//        matrix inversion expression (\f$ \vec{y}=\vec{x}*inv(A) \f$).
// \ingroup math
//
// \param vec The left-hand side transpose dense vector for the multiplication.
// \param mat The right-hand side matrix inversion for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// transpose dense vector and a matrix inversion expression. It restructures the expression
// \f$ \vec{x}=\vec{x}*inv(A) \f$ to the expression \f$ \vec{y}=solve(A^T,\vec{x}) \f$.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseVector<VT,true>& vec, const MatInvExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return solve( trans( (*mat).operand() ), *vec );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix inversion expression
//        and a dense matrix (\f$ C=inv(A)*B \f$).
// \ingroup math
//
// \param lhs The left-hand side matrix inversion.
// \param rhs The right-hand side dense matrix for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// matrix inversion expression and a dense matrix. It restructures the expression
// \f$ C=inv(A)*B \f$ to the expression \f$ C=solve(A,B) \f$.
*/
template< typename MT1  // Matrix base type of the left-hand side expression
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO       // Storage order of the right-hand side dense matrix
        , DisableIf_t< IsMatScalarMultExpr_v<MT2> >* = nullptr >
inline decltype(auto)
   operator*( const MatInvExpr<MT1>& lhs, const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return solve( (*lhs).operand(), *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and a matrix inversion
//        expression (\f$ C=B*inv(A) \f$).
// \ingroup math
//
// \param lhs The left-hand side dense matrix for the multiplication.
// \param rhs The right-hand side matrix inversion for the multiplication.
// \return The resulting vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// dense matrix and a matrix inversion expression. It restructures the expression
// \f$ C=B*inv(A) \f$ to the expression \f$ C=solve(A^T,B^T)^T \f$.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , bool SO         // Storage order of the left-hand side dense matrix
        , typename MT2    // Matrix base type of the right-hand side expression
        , DisableIf_t< IsMatScalarMultExpr_v<MT1> >* = nullptr >
inline decltype(auto)
   operator*( const DenseMatrix<MT1,SO>& lhs, const MatInvExpr<MT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return trans( solve( trans( (*rhs).operand() ), trans( *lhs ) ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculating the inverse of a matrix inversion.
// \ingroup math
//
// \param matrix The matrix to be (re-)inverted.
// \return The inverse of the inverted matrix.
//
// This function implements a performance optimized treatment of the inversion operation on a
// matrix inversion expression. It returns an expression representing the inverse of a matrix
// inversion:

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> A, B;
   // ... Resizing and initialization
   B = inv( inv( A ) );
   \endcode
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) inv( const MatInvExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return matrix.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Computation of the determinant of the given matrix inversion.
// \ingroup math
//
// \param matrix The given matrix inversion.
// \return The determinant of the given matrix inversion.
//
// This function computes the determinant of the given matrix inversion.
//
// \note The computation of the determinant is numerically unreliable since especially for large
// matrices the value can overflow during the computation. Please note that this function does
// not guarantee that it is possible to compute the determinant with the given matrix!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) det( const MatInvExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return inv( det( (*matrix).operand() ) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
