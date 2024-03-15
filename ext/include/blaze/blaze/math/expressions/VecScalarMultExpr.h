//=================================================================================================
/*!
//  \file blaze/math/expressions/VecScalarMultExpr.h
//  \brief Header file for the VecScalarMultExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_VECSCALARMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_VECSCALARMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MultExpr.h>
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
/*!\brief Base class for all vector/scalar multiplication expression templates.
// \ingroup math
//
// The VecScalarMultExpr class serves as a tag for all expression templates that implement a
// vector/scalar multiplication. All classes, that represent a vector/scalar multiplication
// and that are used within the expression template environment of the Blaze library have
// to derive publicly from this class in order to qualify as vector/scalar multiplication
// expression template. Only in case a class is derived publicly from the VecScalarMultExpr
// base class, the IsVecScalarMultExpr type trait recognizes the class as valid vector/scalar
// multiplication expression template.
*/
template< typename VT >  // Vector base type of the expression
struct VecScalarMultExpr
   : public MultExpr<VT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unary minus operator for the negation of a vector-scalar multiplication
//        (\f$ \vec{a} = -(\vec{b} * s) \f$).
// \ingroup math
//
// \param vec The vector-scalar multiplication to be negated.
// \return The negation of the vector-scalar multiplication.
//
// This operator implements a performance optimized treatment of the negation of a vector-scalar
// multiplication expression.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) operator-( const VecScalarMultExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (*vec).leftOperand() * ( -(*vec).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar multiplication
//        expression and a scalar value (\f$ \vec{a}=(\vec{b}*s1)*s2 \f$).
// \ingroup math
//
// \param vec The left-hand side vector-scalar multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// vector-scalar multiplication expression and a scalar value.
*/
template< typename VT  // Vector base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( const VecScalarMultExpr<VT>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return (*vec).leftOperand() * ( (*vec).rightOperand() * scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar multiplication
//        expression and a scalar value (\f$ \vec{a}=s2*(\vec{b}*s1) \f$).
// \ingroup math
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param vec The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a vector-scalar multiplication expression.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename VT  // Vector base type of the expression
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( ST scalar, const VecScalarMultExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return (*vec).leftOperand() * ( scalar * (*vec).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a vector-scalar multiplication expression by
//        a scalar value (\f$ \vec{a}=(\vec{b}*s1)/s2 \f$).
// \ingroup math
//
// \param vec The left-hand side vector-scalar multiplication.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result vector.
//
// This operator implements a performance optimized treatment of the division of a
// vector-scalar multiplication expression by a scalar value.
*/
template< typename VT  // Vector base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> &&
                      ( IsInvertible_v<ST> ||
                        IsInvertible_v< RightOperand_t< VectorType_t<VT> > > ) >* = nullptr >
inline decltype(auto) operator/( const VecScalarMultExpr<VT>& vec, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != ST(0), "Division by zero detected" );

   return (*vec).leftOperand() * ( (*vec).rightOperand() / scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar multiplication
//        expression and a dense vector (\f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$).
// \ingroup math
//
// \param lhs The left-hand side vector-scalar multiplication.
// \param rhs The right-hand side dense vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a vector-scalar multiplication and a dense vector. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1  // Vector base type of the left-hand side expression
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag of the right-hand side dense vector
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT1>& lhs, const DenseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs) ) * (*lhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector and a vector-scalar
//        multiplication expression (\f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$).
// \ingroup math
//
// \param lhs The left-hand side dense vector.
// \param rhs The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense vector and a vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , bool TF         // Transpose flag of the left-hand side dense vector
        , typename VT2 >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseVector<VT1,TF>& lhs, const VecScalarMultExpr<VT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs) * (*rhs).leftOperand() ) * (*rhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar multiplication
//        expression and a sparse vector (\f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side vector-scalar multiplication.
// \param rhs The right-hand side sparse vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a vector-scalar multiplication and a sparse vector. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1  // Vector base type of the left-hand side expression
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF >     // Transpose flag of the left-hand side sparse vector
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT1>& lhs, const SparseVector<VT2,TF>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs) ) * (*lhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector and a vector-scalar
//        multiplication expression (\f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector.
// \param rhs The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse vector and a vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{b}*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*s1 \f$.
*/
template< typename VT1    // Type of the left-hand side sparse vector
        , bool TF         // Transpose flag of the left-hand side sparse vector
        , typename VT2 >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const SparseVector<VT1,TF>& lhs, const VecScalarMultExpr<VT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs) * (*rhs).leftOperand() ) * (*rhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of two vector-scalar multiplication
//        expressions (\f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$).
// \ingroup math
//
// \param lhs The left-hand side vector-scalar multiplication.
// \param rhs The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of two vector-scalar multiplication expressions. It restructures the expression
// \f$ \vec{a}=(\vec{b}*s1)*(\vec{c}*s2) \f$ to the expression \f$ \vec{a}=(\vec{b}*\vec{c})*(s1*s2) \f$.
*/
template< typename VT1    // Vector base type of the left-hand side expression
        , typename VT2 >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT1>& lhs, const VecScalarMultExpr<VT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs).leftOperand() ) * ( (*lhs).rightOperand() * (*rhs).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and a vector-scalar
//        multiplication expression (\f$ \vec{a}=B*(\vec{c}*s1) \f$).
// \ingroup math
//
// \param lhs The left-hand side dense matrix.
// \param rhs The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense matrix and a vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=B*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order of the left-hand side dense matrix
        , typename VT >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseMatrix<MT,SO>& mat, const VecScalarMultExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*mat) * (*vec).leftOperand() ) * (*vec).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose vector-scalar
//        multiplication expression and a dense matrix (\f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side transpose vector-scalar multiplication.
// \param rhs The right-hand side dense matrix.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a transpose vector-scalar multiplication and a dense matrix. It restructures the
// expression \f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$ to the expression \f$ \vec{a}^T=(\vec{b}^T*C)*s1 \f$.
*/
template< typename VT  // Vector base type of the left-hand side expression
        , typename MT  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT>& vec, const DenseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*vec).leftOperand() * (*mat) ) * (*vec).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix and a vector-scalar
//        multiplication expression (\f$ \vec{a}=B*(\vec{c}*s1) \f$).
// \ingroup math
//
// \param lhs The left-hand side sparse matrix.
// \param rhs The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse matrix and a vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=B*(\vec{c}*s1) \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order of the left-hand side sparse matrix
        , typename VT >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const SparseMatrix<MT,SO>& mat, const VecScalarMultExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*mat) * (*vec).leftOperand() ) * (*vec).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a transpose vector-scalar
//        multiplication expression and a sparse matrix (\f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side transpose vector-scalar multiplication.
// \param rhs The right-hand side sparse matrix.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a transpose vector-scalar multiplication and a sparse matrix. It restructures the
// expression \f$ \vec{a}^T=(\vec{b}^T*s1)*C \f$ to the expression \f$ \vec{a}^T=(\vec{b}^T*C)*s1 \f$.
*/
template< typename VT  // Vector base type of the left-hand side expression
        , typename MT  // Type of the right-hand side sparse matrix
        , bool SO >    // Storage order of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT>& vec, const SparseMatrix<MT,SO>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*vec).leftOperand() * (*mat) ) * (*vec).rightOperand();
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
/*!\brief Calculation of the transpose of the given vector-scalar multiplication.
// \ingroup math
//
// \param vector The vector-scalar multiplication expression to be transposed.
// \return The transpose of the expression.
//
// This operator implements the performance optimized treatment of the transpose of a
// vector-scalar multiplication. It restructures the expression \f$ a=trans(b*s1) \f$ to
// the expression \f$ a=trans(b)*s1 \f$.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) trans( const VecScalarMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return trans( (*vector).leftOperand() ) * (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculation of the complex conjugate of the given vector-scalar multiplication.
// \ingroup math
//
// \param vector The vector-scalar multiplication expression to be conjugated.
// \return The complex conjugate of the expression.
//
// This operator implements the performance optimized treatment of the complex conjugate of a
// vector-scalar multiplication. It restructures the expression \f$ a=conj(b*s1) \f$ to the
// expression \f$ a=conj(b)*s1 \f$.
*/
template< typename VT >  // Vector base type of the expression
inline decltype(auto) conj( const VecScalarMultExpr<VT>& vector )
{
   BLAZE_FUNCTION_TRACE;

   return conj( (*vector).leftOperand() ) * (*vector).rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
