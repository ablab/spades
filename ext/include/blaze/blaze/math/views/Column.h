//=================================================================================================
/*!
//  \file blaze/math/views/Column.h
//  \brief Header file for the implementation of the Column view
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

#ifndef _BLAZE_MATH_VIEWS_COLUMN_H_
#define _BLAZE_MATH_VIEWS_COLUMN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/UniformVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DeclExpr.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatKronExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatMultExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/MatNoAliasExpr.h>
#include <blaze/math/expressions/MatNoSIMDExpr.h>
#include <blaze/math/expressions/MatRepeatExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/VecExpandExpr.h>
#include <blaze/math/expressions/VecTVecMapExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/functors/Bind2nd.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/column/ColumnData.h>
#include <blaze/math/views/column/Dense.h>
#include <blaze/math/views/column/Sparse.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given matrix.
// \ingroup column
//
// \param matrix The matrix containing the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column<3UL>( D );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column<4UL>( S );
   \endcode

// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than or equal to the total number
// of the columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto col3 = column<3UL>( D, unchecked );
   auto col4 = column<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Column_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant matrix.
// \ingroup column
//
// \param matrix The constant matrix containing the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column<3UL>( D );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column<4UL>( S );
   \endcode

// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than or equal to the total number
// of the columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto col3 = column<3UL>( D, unchecked );
   auto col4 = column<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Column_<const MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given temporary matrix.
// \ingroup column
//
// \param matrix The temporary matrix containing the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// matrix. In case the column is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // Column index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Matrix<MT,SO>&& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Column_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given matrix.
// \ingroup column
//
// \param matrix The matrix containing the column.
// \param index The index of the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column( S, 4UL );
   \endcode

// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than or equal to the total number
// of the columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto col3 = column( D, 3UL, unchecked );
   auto col4 = column( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Matrix<MT,SO>& matrix, size_t index, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Column_<MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given constant matrix.
// \ingroup column
//
// \param matrix The constant matrix containing the column.
// \param index The index of the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 3rd column of the dense matrix D
   auto col3 = column( D, 3UL );

   // Creating a view on the 4th column of the sparse matrix S
   auto col4 = column( S, 4UL );
   \endcode

// By default, the provided column arguments are checked at runtime. In case the column is not
// properly specified (i.e. if the specified index is greater than or equal to the total number
// of the columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto col3 = column( D, 3UL, unchecked );
   auto col4 = column( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const Matrix<MT,SO>& matrix, size_t index, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Column_<const MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific column of the given temporary matrix.
// \ingroup column
//
// \param matrix The temporary matrix containing the column.
// \param index The index of the column.
// \param args Optional column arguments.
// \return View on the specified column of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given temporary
// matrix. In case the column is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( Matrix<MT,SO>&& matrix, size_t index, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Column_<MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix addition.
// \ingroup column
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime column arguments.
// \return View on the specified column of the addition.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// addition.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMatAddExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).leftOperand(), args... ) +
          column<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix subtraction.
// \ingroup column
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime column arguments.
// \return View on the specified column of the subtraction.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// subtraction.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMatSubExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).leftOperand(), args... ) -
          column<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given Schur product.
// \ingroup column
//
// \param matrix The constant Schur product.
// \param args The runtime column arguments.
// \return View on the specified column of the Schur product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given Schur
// product.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const SchurExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).leftOperand(), args... ) *
          column<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/matrix multiplication.
// \ingroup column
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime column arguments.
// \return View on the specified column of the multiplication.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix/matrix
// multiplication.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMatMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return (*matrix).leftOperand() * column<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given Kronecker product.
// \ingroup column
//
// \param matrix The constant Kronecker product.
// \param args The runtime column arguments.
// \return View on the specified column of the Kronecker product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given Kronecker
// product.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMatKronExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs...> cd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   const size_t columns( (*matrix).rightOperand().columns() );

   return kron( column( (*matrix).leftOperand(), cd.column()/columns, unchecked ),
                column( (*matrix).rightOperand(), cd.column()%columns, unchecked ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer product.
// \ingroup column
//
// \param matrix The constant outer product.
// \param args Optional column arguments.
// \return View on the specified column of the outer product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given outer
// product.
*/
template< size_t I            // Column index
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const VecTVecMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return (*matrix).leftOperand() * (*matrix).rightOperand()[I];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer product.
// \ingroup column
//
// \param matrix The constant outer product.
// \param index The index of the column.
// \param args Optional column arguments.
// \return View on the specified column of the outer product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given outer
// product.
*/
template< typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const VecTVecMultExpr<MT>& matrix, size_t index, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return (*matrix).leftOperand() * (*matrix).rightOperand()[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar multiplication.
// \ingroup column
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime column arguments.
// \return View on the specified column of the multiplication.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// multiplication.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatScalarMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).leftOperand(), args... ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix/scalar division.
// \ingroup column
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime column arguments.
// \return View on the specified column of the division.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix/scalar
// division.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatScalarDivExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).leftOperand(), args... ) / (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given unary matrix map operation.
// \ingroup column
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified column of the unary map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given unary
// matrix map operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( column<CCAs...>( (*matrix).operand(), args... ), (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given binary matrix map operation.
// \ingroup column
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified column of the binary map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given binary
// matrix map operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatMatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( column<CCAs...>( (*matrix).leftOperand(), args... ),
               column<CCAs...>( (*matrix).rightOperand(), args... ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer map operation.
// \ingroup column
//
// \param matrix The constant outer map operation.
// \param args Optional column arguments.
// \return View on the specified column of the outer map operation.
// \exception std::invalid_argument Invalid column access index.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given outer
// map operation.
*/
template< size_t I            // Column index
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const VecTVecMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return map( (*matrix).leftOperand(),
               blaze::bind2nd( (*matrix).operation(), (*matrix).rightOperand()[I] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given outer map operation.
// \ingroup column
//
// \param matrix The constant outer map operation.
// \param index The index of the column.
// \param args Optional column arguments.
// \return View on the specified column of the outer map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given outer
// map operation.
*/
template< typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) column( const VecTVecMapExpr<MT>& matrix, size_t index, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return map( (*matrix).leftOperand(),
               blaze::bind2nd( (*matrix).operation(), (*matrix).rightOperand()[index] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix evaluation operation.
// \ingroup column
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime column arguments.
// \return View on the specified column of the evaluation operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// evaluation operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatEvalExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( column<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix serialization operation.
// \ingroup column
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime column arguments.
// \return View on the specified column of the serialization operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// serialization operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatSerialExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( column<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix no-alias operation.
// \ingroup column
//
// \param matrix The constant matrix no-alias operation.
// \param args The runtime column arguments.
// \return View on the specified column of the no-alias operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// no-alias operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatNoAliasExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( column<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix no-SIMD operation.
// \ingroup column
//
// \param matrix The constant matrix no-SIMD operation.
// \param args The runtime column arguments.
// \return View on the specified column of the no-SIMD operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// no-SIMD operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatNoSIMDExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( column<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix declaration operation.
// \ingroup column
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime column arguments.
// \return View on the specified column of the declaration operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// declaration operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const DeclExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column<CCAs...>( (*matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix transpose operation.
// \ingroup column
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime column arguments.
// \return View on the specified column of the transpose operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given matrix
// transpose operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatTransExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return trans( row<CCAs...>( (*matrix).operand(), args... ) );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given column-major vector expansion operation.
// \ingroup column
//
// \param matrix The constant vector expansion operation.
// \param args The runtime column arguments.
// \return View on the specified column of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column-major
// vector expansion operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , size_t... CEAs      // Compile time expansion arguments
        , typename... RCAs    // Runtime column arguments
        , EnableIf_t< IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) column( const VecExpandExpr<MT,CEAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      const ColumnData<CCAs...> cd( args... );
      if( (*matrix).columns() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return subvector( (*matrix).operand(), 0UL, (*matrix).rows(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given row-major vector expansion operation.
// \ingroup column
//
// \param matrix The constant vector expansion operation.
// \param args The runtime column arguments.
// \return View on the specified column of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given row-major
// vector expansion operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , size_t... CEAs      // Compile time expansion arguments
        , typename... RCAs    // Runtime column arguments
        , EnableIf_t< !IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) column( const VecExpandExpr<MT,CEAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs...> cd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   using ET = ElementType_t< MatrixType_t<MT> >;

   return UniformVector<ET,columnVector>( (*matrix).rows(), (*matrix).operand()[cd.column()] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given matrix repeat operation.
// \ingroup column
//
// \param matrix The constant matrix repeat operation.
// \param args The runtime column arguments.
// \return View on the specified column of the repeat operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column-major
// matrix repeat operation.
*/
template< size_t... CCAs      // Compile time column arguments
        , typename MT         // Matrix base type of the expression
        , size_t... CRAs      // Compile time repeater arguments
        , typename... RCAs >  // Runtime column arguments
inline decltype(auto) column( const MatRepeatExpr<MT,CRAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs...> cd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).columns() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return repeat( column( (*matrix).operand()
                        , cd.column() % (*matrix).operand().columns()
                        , unchecked )
                , (*matrix).template repetitions<0UL>() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMN OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given column.
// \ingroup column
//
// \param column The column to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void reset( Column<MT,SO,DF,SF,CCAs...>& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary column.
// \ingroup column
//
// \param column The temporary column to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void reset( Column<MT,SO,DF,SF,CCAs...>&& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given column.
// \ingroup column
//
// \param column The column to be cleared.
// \return void
//
// Clearing a column is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void clear( Column<MT,SO,DF,SF,CCAs...>& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary column.
// \ingroup column
//
// \param column The temporary column to be cleared.
// \return void
//
// Clearing a column is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline void clear( Column<MT,SO,DF,SF,CCAs...>&& column )
{
   column.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense column is in default state.
// \ingroup column
//
// \param column The dense column to be tested for its default state.
// \return \a true in case the given dense column is component-wise zero, \a false otherwise.
//
// This function checks whether the dense column is in default state. For instance, in case the
// column is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all column elements are 0 and \a false in case any column element is not 0.
// The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO            // Storage order
        , bool SF            // Symmetry flag
        , size_t... CCAs >   // Compile time column arguments
inline bool isDefault( const Column<MT,SO,true,SF,CCAs...>& column )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<column.size(); ++i )
      if( !isDefault<RF>( column[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse column is in default state.
// \ingroup column
//
// \param column The sparse column to be tested for its default state.
// \return \a true in case the given column is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse column is in default state. For instance, in case the
// column is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all column elements are 0 and \a false in case any column element is not 0.
// The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( column( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO            // Storage order
        , bool SF            // Symmetry flag
        , size_t... CCAs >   // Compile time column arguments
inline bool isDefault( const Column<MT,SO,false,SF,CCAs...>& column )
{
   using blaze::isDefault;

   for( const auto& element : column )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given column are intact.
// \ingroup column
//
// \param column The column to be tested.
// \return \a true in case the given column's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the column are intact, i.e. if its state is
// valid. In case the invariants are intact, the function returns \a true, else it will return
// \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,columnMajor> A;
   // ... Resizing and initialization
   if( isIntact( column( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool isIntact( const Column<MT,SO,DF,SF,CCAs...>& column ) noexcept
{
   return ( column.column() < column.operand().columns() &&
            isIntact( column.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given columns represent the same observable state.
// \ingroup column
//
// \param a The first column to be tested for its state.
// \param b The second column to be tested for its state.
// \return \a true in case the two columns share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given columns refer to exactly the
// same range of the same matrix. In case both columns represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side column
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side column
        , size_t... CCAs1    // Compile time column arguments of the left-hand side column
        , typename MT2       // Type of the matrix of the right-hand side column
        , bool SF2           // Symmetry flag of the right-hand side column
        , size_t... CCAs2 >  // Compile time column arguments of the right-hand side column
inline bool isSame( const Column<MT1,SO,DF,SF1,CCAs1...>& a,
                    const Column<MT2,SO,DF,SF2,CCAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) && ( a.column() == b.column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool trySet( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return trySet( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The value to be set to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return trySet( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryAdd( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryAdd( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The value to be added to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryAdd( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool trySub( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return trySub( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The value to be subtracted from the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return trySub( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryMult( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryMult( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The factor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryMult( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryDiv( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryDiv( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The divisor for the elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryDiv( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param count The number of bits to shift the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline bool tryShift( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, int count )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryShift( column.operand(), index, column.column(), count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param count The number of bits to shift the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE bool
   tryShift( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, int count )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryShift( column.operand(), index, column.column(), size, 1UL, count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryBitand( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryBitand( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitand( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryBitand( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryBitor( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryBitor( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryBitor( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a column.
// \ingroup column
//
// \param sv The target column.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
inline bool tryBitxor( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < column.size(), "Invalid vector access index" );

   return tryBitxor( column.operand(), index, column.column(), value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a column.
// \ingroup column
//
// \param column The target column.
// \param index The index of the first element of the range to be modified.
// \param size The number of elements of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const Column<MT,SO,DF,SF,CCAs...>& column, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*column).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*column).size(), "Invalid range size" );

   return tryBitxor( column.operand(), index, column.column(), size, 1UL, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Column<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Column<MT,SO,SF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAddAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Column<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return trySubAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Column<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryMultAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Column<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryDivAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector of bits to shift.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryShiftAssign( const Column<MT,SO,SF,SF,CCAs...>& lhs,
                            const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryShiftAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector for the bitwise AND operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitandAssign( const Column<MT,SO,SF,SF,CCAs...>& lhs,
                             const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitandAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector for the bitwise OR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitorAssign( const Column<MT,SO,SF,SF,CCAs...>& lhs,
                            const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitorAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to a column.
// \ingroup column
//
// \param lhs The target left-hand side column.
// \param rhs The right-hand side vector for the bitwise XOR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT     // Type of the matrix
        , bool SO         // Storage order
        , bool DF         // Density flag
        , bool SF         // Symmetry flag
        , size_t... CCAs  // Compile time column arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitxorAssign( const Column<MT,SO,SF,SF,CCAs...>& lhs,
                             const Vector<VT,false>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitxorAssign( lhs.operand(), *rhs, index, lhs.column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column.
// \ingroup column
//
// \param c The column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given column. It returns a
// column object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Compile time column arguments
inline decltype(auto) derestrict( Column<MT,SO,DF,SF,I>& c )
{
   return column<I>( derestrict( c.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column.
// \ingroup column
//
// \param c The temporary column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column. It
// returns a column object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Compile time column arguments
inline decltype(auto) derestrict( Column<MT,SO,DF,SF,I>&& c )
{
   return column<I>( derestrict( c.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column.
// \ingroup column
//
// \param c The column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given column. It returns a
// column object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Column<MT,SO,DF,SF>& c )
{
   return column( derestrict( c.operand() ), c.column(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column.
// \ingroup column
//
// \param c The temporary column to be derestricted.
// \return Column without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column. It
// returns a column object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Column<MT,SO,DF,SF>&& c )
{
   return column( derestrict( c.operand() ), c.column(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given column.
// \ingroup column
//
// \param c The given column.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given column.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline decltype(auto) unview( Column<MT,SO,DF,SF,CCAs...>& c )
{
   return c.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given constant column.
// \ingroup column
//
// \param c The given constant column.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given constant column.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CCAs >  // Compile time column arguments
inline decltype(auto) unview( const Column<MT,SO,DF,SF,CCAs...>& c )
{
   return c.operand();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct Size< Column<MT,SO,DF,SF,CRAs...>, 0UL >
   : public Size<MT,0UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAXSIZE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct MaxSize< Column<MT,SO,DF,SF,CRAs...>, 0UL >
   : public MaxSize<MT,0UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESTRICTED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF, bool SF, size_t... CCAs >
struct IsRestricted< Column<MT,SO,DF,SF,CCAs...> >
   : public IsRestricted<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASCONSTDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct HasConstDataAccess< Column<MT,SO,true,SF,CCAs...> >
   : public HasConstDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HASMUTABLEDATAACCESS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct HasMutableDataAccess< Column<MT,SO,true,SF,CCAs...> >
   : public HasMutableDataAccess<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISALIGNED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct IsAligned< Column<MT,SO,true,SF,CCAs...> >
   : public BoolConstant< IsAligned_v<MT> && ( IsColumnMajorMatrix_v<MT> || IsSymmetric_v<MT> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SF, size_t... CCAs >
struct IsContiguous< Column<MT,true,true,SF,CCAs...> >
   : public IsContiguous<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISPADDED SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool SF, size_t... CCAs >
struct IsPadded< Column<MT,SO,true,SF,CCAs...> >
   : BoolConstant< IsPadded_v<MT> && ( IsColumnMajorMatrix_v<MT> || IsSymmetric_v<MT> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISOPPOSEDVIEW SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool DF, size_t... CCAs >
struct IsOpposedView< Column<MT,false,DF,false,CCAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
