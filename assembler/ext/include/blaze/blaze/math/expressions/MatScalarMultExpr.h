//=================================================================================================
/*!
//  \file blaze/math/expressions/MatScalarMultExpr.h
//  \brief Header file for the MatScalarMultExpr base class
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

#ifndef _BLAZE_MATH_EXPRESSIONS_MATSCALARMULTEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_MATSCALARMULTEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
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
/*!\brief Base class for all matrix/scalar multiplication expression templates.
// \ingroup math
//
// The MatScalarMultExpr class serves as a tag for all expression templates that implement a
// matrix/scalar multiplication. All classes, that represent a matrix/scalar multiplication
// and that are used within the expression template environment of the Blaze library have
// to derive publicly from this class in order to qualify as matrix/scalar multiplication
// expression template. Only in case a class is derived publicly from the MatScalarMultExpr
// base class, the IsMatScalarMultExpr type trait recognizes the class as valid matrix/scalar
// multiplication expression template.
*/
template< typename MT >  // Matrix base type of the expression
struct MatScalarMultExpr
   : public MultExpr<MT>
{};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING UNARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unary minus operator for the negation of a matrix-scalar multiplication
//        (\f$ A = -(B*s) \f$).
// \ingroup math
//
// \param mat The matrix-scalar multiplication to be negated.
// \return The negation of the matrix-scalar multiplication.
//
// This operator implements a performance optimized treatment of the negation of a matrix-scalar
// multiplication expression.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) operator-( const MatScalarMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return (*mat).leftOperand() * ( -(*mat).rightOperand() );
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
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a scalar value (\f$ A=(B*s1)*s2 \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-scalar multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// matrix-scalar multiplication expression and a scalar value.
*/
template< typename MT  // Matrix base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( const MatScalarMultExpr<MT>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   return (*mat).leftOperand() * ( (*mat).rightOperand() * scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a scalar value and a matrix-scalar
//        multiplication expression (\f$ A=s2*(B*s1) \f$).
// \ingroup math
//
// \param scalar The left-hand side scalar value for the multiplication.
// \param mat The right-hand side matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the multiplication of a
// scalar value and a matrix-scalar multiplication expression.
*/
template< typename ST  // Type of the left-hand side scalar
        , typename MT  // Matrix base type of the expression
        , EnableIf_t< IsScalar_v<ST> >* = nullptr >
inline decltype(auto) operator*( ST scalar, const MatScalarMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return (*mat).leftOperand() * ( scalar * (*mat).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division operator for the division of a matrix-scalar multiplication expression by
//        a scalar value (\f$ A=(B*s1)/s2 \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-scalar multiplication.
// \param scalar The right-hand side scalar value for the division.
// \return The scaled result matrix.
//
// This operator implements a performance optimized treatment of the division of a
// matrix-scalar multiplication expression by a scalar value.
*/
template< typename MT  // Matrix base type of the expression
        , typename ST  // Type of the right-hand side scalar
        , EnableIf_t< IsScalar_v<ST> &&
                      ( IsInvertible_v<ST> ||
                        IsInvertible_v< RightOperand_t< MatrixType_t<MT> > > ) >* = nullptr >
inline decltype(auto) operator/( const MatScalarMultExpr<MT>& mat, ST scalar )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_USER_ASSERT( scalar != ST(0), "Division by zero detected" );

   return (*mat).leftOperand() * ( (*mat).rightOperand() / scalar );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a dense vector (\f$ \vec{a}=(B*s1)*\vec{c} \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-scalar multiplication.
// \param vec The right-hand side dense vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a matrix-scalar multiplication and a dense vector. It restructures the expression
// \f$ \vec{a}=(B*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Type of the right-hand side dense vector
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT>& mat, const DenseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*mat).leftOperand() * (*vec) ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense vector and a
//        matrix-scalar multiplication expression (\f$ \vec{a}^T=\vec{c}^T*(B*s1) \f$).
// \ingroup math
//
// \param vec The left-hand side dense vector.
// \param mat The right-hand side matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense vector and a matrix-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{c}^T*(B*s1) \f$ to the expression \f$ \vec{a}^T=(\vec{c}^T*B)*s1 \f$.
*/
template< typename VT    // Type of the left-hand side dense vector
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseVector<VT,true>& vec, const MatScalarMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*vec) * (*mat).leftOperand() ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a sparse vector (\f$ \vec{a}=(B*s1)*\vec{c} \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-scalar multiplication.
// \param vec The right-hand side sparse vector.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a matrix-scalar multiplication and a sparse vector. It restructures the expression
// \f$ \vec{a}=(B*s1)*\vec{c} \f$ to the expression \f$ \vec{a}=(B*\vec{c})*s1 \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Type of the right-hand side sparse vector
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT>& mat, const SparseVector<VT,false>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*mat).leftOperand() * (*vec) ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse vector and a
//        matrix-scalar multiplication expression (\f$ \vec{a}^T=\vec{c}^T*(B*s1) \f$).
// \ingroup math
//
// \param vec The left-hand side sparse vector.
// \param mat The right-hand side matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse vector and a matrix-scalar multiplication. It restructures the expression
// \f$ \vec{a}=\vec{c}^T*(B*s1) \f$ to the expression \f$ \vec{a}^T=(\vec{c}^T*B)*s1 \f$.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const SparseVector<VT,true>& vec, const MatScalarMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*vec) * (*mat).leftOperand() ) * (*mat).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a vector-scalar multiplication expression
//        (\f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$).
// \ingroup math
//
// \param mat The left-hand side matrix-scalar multiplication.
// \param vec The right-hand side vector-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a matrix-
// scalar multiplication and a vector-scalar multiplication. It restructures the expression
// \f$ \vec{a}=(B*s1)*(\vec{c}*s2) \f$ to the expression \f$ \vec{a}=(B*\vec{c})*(s1*s2) \f$.
*/
template< typename MT    // Matrix base type of the left-hand side expression
        , typename VT >  // Vector base type of the right-hand side expression
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT>& mat, const VecScalarMultExpr<VT>& vec )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*mat).leftOperand() * (*vec).leftOperand() ) * ( (*mat).rightOperand() * (*vec).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a vector-scalar multiplication
//        expression and a matrix-scalar multiplication expression
//        (\f$ \vec{a}^T=\vec{b}^T*(C*s1) \f$).
// \ingroup math
//
// \param vec The left-hand side vector-scalar multiplication.
// \param mat The right-hand side matrix-scalar multiplication.
// \return The scaled result vector.
//
// This operator implements the performance optimized treatment of the multiplication of a vector-
// scalar multiplication and a matrix-scalar multiplication. It restructures the expression
// \f$ \vec{a}=(\vec{b}^T*s1)*(C*s2) \f$ to the expression \f$ \vec{a}^T=(\vec{b}^T*C)*(s1*s2) \f$.
*/
template< typename VT    // Vector base type of the left-hand side expression
        , typename MT >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const VecScalarMultExpr<VT>& vec, const MatScalarMultExpr<MT>& mat )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*vec).leftOperand() * (*mat).leftOperand() ) * ( (*vec).rightOperand() * (*mat).rightOperand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a dense matrix (\f$ A=(B*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side matrix-scalar multiplication.
// \param rhs The right-hand side dense matrix.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication
// of a matrix-scalar multiplication and a dense matrix. It restructures the expression
// \f$ A=(B*s1)*C \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1  // Matrix base type of the left-hand side expression
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO >     // Storage order of the right-hand side dense matrix
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT1>& lhs, const DenseMatrix<MT2,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs) ) * (*lhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a dense matrix and a matrix-scalar
//        multiplication expression (\f$ A=(B*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side dense matrix.
// \param rhs The right-hand side matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication
// of a dense matrix and a matrix-scalar multiplication. It restructures the expression
// \f$ A=B*(C*s1) \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , bool SO         // Storage order of the left-hand side dense matrix
        , typename MT2 >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const DenseMatrix<MT1,SO>& lhs, const MatScalarMultExpr<MT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs) * (*rhs).leftOperand() ) * (*rhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a matrix-scalar multiplication
//        expression and a sparse matrix (\f$ A=(B*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side matrix-scalar multiplication.
// \param rhs The right-hand side sparse matrix.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication
// of a matrix-scalar multiplication and a sparse matrix. It restructures the expression
// \f$ A=(B*s1)*C \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1  // Matrix base type of the left-hand side expression
        , typename MT2  // Type of the right-hand side sparse matrix
        , bool SO >     // Storage order of the right-hand side sparse matrix
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT1>& lhs, const SparseMatrix<MT2,SO>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs) ) * (*lhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of a sparse matrix and a matrix-scalar
//        multiplication expression (\f$ A=(B*s1)*C \f$).
// \ingroup math
//
// \param lhs The left-hand side sparse  matrix.
// \param rhs The right-hand side matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication
// of a sparse matrix and a matrix-scalar multiplication. It restructures the expression
// \f$ A=B*(C*s1) \f$ to the expression \f$ A=(B*C)*s1 \f$.
*/
template< typename MT1    // Type of the left-hand side sparse matrix
        , bool SO         // Storage order of the left-hand side sparse matrix
        , typename MT2 >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const SparseMatrix<MT1,SO>& lhs, const MatScalarMultExpr<MT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs) * (*rhs).leftOperand() ) * (*rhs).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication operator for the multiplication of two matrix-scalar multiplication
//        expressions (\f$ A=(B*s1)*(C*s2) \f$).
// \ingroup math
//
// \param lhs The left-hand side matrix-scalar multiplication.
// \param rhs The right-hand side matrix-scalar multiplication.
// \return The scaled result matrix.
//
// This operator implements the performance optimized treatment of the multiplication of
// two matrix-scalar multiplication expressions. It restructures the expression
// \f$ A=(B*s1)*(C*s2) \f$ to the expression \f$ A=(B*C)*(s1*s2) \f$.
*/
template< typename MT1    // Matrix base type of the left-hand side expression
        , typename MT2 >  // Matrix base type of the right-hand side expression
inline decltype(auto)
   operator*( const MatScalarMultExpr<MT1>& lhs, const MatScalarMultExpr<MT2>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   return ( (*lhs).leftOperand() * (*rhs).leftOperand() ) * ( (*lhs).rightOperand() * (*rhs).rightOperand() );
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
/*!\brief Calculation of the transpose of the given matrix-scalar multiplication.
// \ingroup math
//
// \param matrix The matrix-scalar multiplication expression to be transposed.
// \return The transpose of the expression.
//
// This operator implements the performance optimized treatment of the transpose of a
// matrix-scalar multiplication. It restructures the expression \f$ A=trans(B*s1) \f$ to
// the expression \f$ A=trans(B)*s1 \f$.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) trans( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return trans( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Calculation of the complex conjugate of the given matrix-scalar multiplication.
// \ingroup math
//
// \param matrix The matrix-scalar multiplication expression to be conjugated.
// \return The complex conjugate of the expression.
//
// This operator implements the performance optimized treatment of the complex conjugate of a
// matrix-scalar multiplication. It restructures the expression \f$ a=conj(b*s1) \f$ to the
// expression \f$ a=conj(b)*s1 \f$.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) conj( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   return conj( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-symmetric matrix-scalar multiplication expression as symmetric.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid symmetric matrix specification.
//
// This function implements the application of the \a declsym operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declsym(B*s1) \f$ to the expression
// \f$ A=declsym(B)*s1 \f$. In case the given matrix is not a square matrix, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declsym( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid symmetric matrix specification" );
   }

   return declsym( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-Hermitian matrix-scalar multiplication expression as Hermitian.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid Hermitian matrix specification.
//
// This function implements the application of the declherm() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declherm(B*s1) \f$ to the expression
// \f$ A=declherm(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declherm( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid Hermitian matrix specification" );
   }

   return declherm( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-lower matrix-scalar multiplication expression as lower.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid lower matrix specification.
//
// This function implements the application of the decllow() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=decllow(B*s1) \f$ to the expression
// \f$ A=decllow(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) decllow( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid lower matrix specification" );
   }

   return decllow( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-unilower matrix-scalar multiplication expression as unilower.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid unilower matrix specification.
//
// This function implements the application of the declunilow() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declunilow(B*s1) \f$ to the expression
// \f$ A=declunilow(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declunilow( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid unilower matrix specification" );
   }

   return declunilow( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-strictly-lower matrix-scalar multiplication expression as
//        strictly lower.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid strictly lower matrix specification.
//
// This function implements the application of the declstrlow() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declstrlow(B*s1) \f$ to the expression
// \f$ A=declstrlow(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declstrlow( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid strictly lower matrix specification" );
   }

   return declstrlow( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-upper matrix-scalar multiplication expression as upper.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid upper matrix specification.
//
// This function implements the application of the declupp() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declupp(B*s1) \f$ to the expression
// \f$ A=declupp(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declupp( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid upper matrix specification" );
   }

   return declupp( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-uniupper matrix-scalar multiplication expression as uniupper.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid uniupper matrix specification.
//
// This function implements the application of the decluniupp() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=decluniupp(B*s1) \f$ to the expression
// \f$ A=decluniupp(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) decluniupp( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid uniupper matrix specification" );
   }

   return decluniupp( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-strictly-upper matrix-scalar multiplication expression as
//        strictly upper.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid strictly upper matrix specification.
//
// This function implements the application of the declstrupp() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=declstrupp(B*s1) \f$ to the expression
// \f$ A=declstrupp(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) declstrupp( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid strictly upper matrix specification" );
   }

   return declstrupp( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Declares the given non-diagonal matrix-scalar multiplication expression as diagonal.
// \ingroup math
//
// \param matrix The input matrix-scalar multiplication expression.
// \return The redeclared expression.
// \exception std::invalid_argument Invalid diagonal matrix specification.
//
// This function implements the application of the decldiag() operation on a matrix-scalar
// multiplication. It restructures the expression \f$ A=decldiag(B*s1) \f$ to the expression
// \f$ A=decldiag(B)*s1 \f$. In case the given matrix is not a square matrix,
// a \a std::invalid_argument exception is thrown.
*/
template< typename MT >  // Matrix base type of the expression
inline decltype(auto) decldiag( const MatScalarMultExpr<MT>& matrix )
{
   BLAZE_FUNCTION_TRACE;

   if( !isSquare( *matrix ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid diagonal matrix specification" );
   }

   return decldiag( (*matrix).leftOperand() ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
