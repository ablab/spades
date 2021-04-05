//=================================================================================================
/*!
//  \file blaze/math/views/Row.h
//  \brief Header file for the implementation of the Row view
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

#ifndef _BLAZE_MATH_VIEWS_ROW_H_
#define _BLAZE_MATH_VIEWS_ROW_H_


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
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/row/Dense.h>
#include <blaze/math/views/row/RowData.h>
#include <blaze/math/views/row/Sparse.h>
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
/*!\brief Creating a view on a specific row of the given matrix.
// \ingroup row
//
// \param matrix The matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row<3UL>( D );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row<4UL>( S );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row<3UL>( D, unchecked );
   auto row4 = row<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given constant matrix.
// \ingroup row
//
// \param matrix The constant matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row<3UL>( D );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row<4UL>( S );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row<3UL>( D, unchecked );
   auto row4 = row<4UL>( S, unchecked );
   \endcode
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Matrix<MT,SO>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Row_<const MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given temporary matrix.
// \ingroup row
//
// \param matrix The temporary matrix containing the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary
// matrix. In case the row is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // Row index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>&& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given matrix.
// \ingroup row
//
// \param matrix The matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row( D, 3UL );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row( S, 4UL );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row( D, 3UL, unchecked );
   auto row4 = row( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given constant matrix.
// \ingroup row
//
// \param matrix The constant matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the 3rd row of the dense matrix D
   auto row3 = row( D, 3UL );

   // Creating a view on the 4th row of the sparse matrix S
   auto row4 = row( S, 4UL );
   \endcode

// By default, the provided row arguments are checked at runtime. In case the row is not properly
// specified (i.e. if the specified index is greater than or equal to the total number of the rows
// in the given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped
// by providing the optional \a blaze::unchecked argument.

   \code
   auto row3 = row( D, 3UL, unchecked );
   auto row4 = row( S, 4UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const Matrix<MT,SO>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Row_<const MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific row of the given temporary matrix.
// \ingroup row
//
// \param matrix The temporary matrix containing the row.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the matrix.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given temporary
// matrix. In case the row is not properly specified (i.e. if the specified index is greater
// than or equal to the total number of the rows in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( Matrix<MT,SO>&& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Row_<MT>;
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
/*!\brief Creating a view on a specific row of the given matrix/matrix addition.
// \ingroup row
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime row arguments.
// \return View on the specified row of the addition.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// addition.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatAddExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) +
          row<CRAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/matrix subtraction.
// \ingroup row
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime row arguments.
// \return View on the specified row of the subtraction.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// subtraction.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatSubExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) -
          row<CRAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given Schur product.
// \ingroup row
//
// \param matrix The constant Schur product.
// \param args The runtime row arguments.
// \return View on the specified row of the Schur product.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given Schur product.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const SchurExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) *
          row<CRAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/matrix multiplication.
// \ingroup row
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime row arguments.
// \return View on the specified row of the multiplication.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix/matrix
// multiplication.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given Kronecker product.
// \ingroup row
//
// \param matrix The constant Kronecker product.
// \param args The runtime row arguments.
// \return View on the specified row of the Kronecker product.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given Kronecker
// product.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatKronExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs...> rd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   const size_t rows( (*matrix).rightOperand().rows() );

   return kron( row( (*matrix).leftOperand(), rd.row()/rows, unchecked ),
                row( (*matrix).rightOperand(), rd.row()%rows, unchecked ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer product.
// \ingroup row
//
// \param matrix The constant outer product.
// \param args Optional row arguments.
// \return View on the specified row of the outer product.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given outer product.
*/
template< size_t I            // Row index
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return (*matrix).leftOperand()[I] * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer product.
// \ingroup row
//
// \param matrix The constant outer product.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the outer product.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given outer product.
*/
template< typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMultExpr<MT>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return (*matrix).leftOperand()[index] * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/scalar multiplication.
// \ingroup row
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime row arguments.
// \return View on the specified row of the multiplication.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix/scalar
// multiplication.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatScalarMultExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix/scalar division.
// \ingroup row
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime row arguments.
// \return View on the specified row of the division.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix/scalar
// division.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatScalarDivExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).leftOperand(), args... ) / (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given unary matrix map operation.
// \ingroup row
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime row arguments.
// \return View on the specified row of the unary map operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given unary matrix
// map operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( row<CRAs...>( (*matrix).operand(), args... ), (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given binary matrix map operation.
// \ingroup row
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime row arguments.
// \return View on the specified row of the binary map operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given binary matrix
// map operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatMatMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( row<CRAs...>( (*matrix).leftOperand(), args... ),
               row<CRAs...>( (*matrix).rightOperand(), args... ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer map operation.
// \ingroup row
//
// \param matrix The constant outer map operation.
// \param args Optional row arguments.
// \return View on the specified row of the outer map operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given outer map
// operation.
*/
template< size_t I            // Row index
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMapExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= I ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return map( (*matrix).rightOperand(),
               blaze::bind2nd( (*matrix).operation(), (*matrix).leftOperand()[I] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given outer map operation.
// \ingroup row
//
// \param matrix The constant outer map operation.
// \param index The index of the row.
// \param args Optional row arguments.
// \return View on the specified row of the outer map operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given outer map
// operation.
*/
template< typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Optional row arguments
inline decltype(auto) row( const VecTVecMapExpr<MT>& matrix, size_t index, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= index ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return map( (*matrix).rightOperand(),
               blaze::bind2nd( (*matrix).operation(), (*matrix).leftOperand()[index] ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix evaluation operation.
// \ingroup row
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime row arguments.
// \return View on the specified row of the evaluation operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// evaluation operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatEvalExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( row<CRAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix serialization operation.
// \ingroup row
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime row arguments.
// \return View on the specified row of the serialization operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// serialization operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatSerialExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( row<CRAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix no-alias operation.
// \ingroup row
//
// \param matrix The constant matrix no-alias operation.
// \param args The runtime row arguments.
// \return View on the specified row of the no-alias operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// no-alias operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatNoAliasExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( row<CRAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix no-SIMD operation.
// \ingroup row
//
// \param matrix The constant matrix no-SIMD operation.
// \param args The runtime row arguments.
// \return View on the specified row of the no-SIMD operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// no-SIMD operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatNoSIMDExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( row<CRAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix declaration operation.
// \ingroup row
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime row arguments.
// \return View on the specified row of the declaration operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// declaration operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const DeclExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return row<CRAs...>( (*matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix transpose operation.
// \ingroup row
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime row arguments.
// \return View on the specified row of the transpose operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given matrix
// transpose operation.
*/
template< size_t... CRAs      // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatTransExpr<MT>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return trans( column<CRAs...>( (*matrix).operand(), args... ) );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given row-major vector expansion operation.
// \ingroup row
//
// \param matrix The constant vector expansion operation.
// \param args The runtime row arguments.
// \return View on the specified row of the expansion operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given row-major
// vector expansion operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename... RRAs  // Runtime row arguments
        , EnableIf_t< IsRowMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) row( const VecExpandExpr<MT,CEAs...>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      const RowData<CRAs...> rd( args... );
      if( (*matrix).rows() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return subvector( (*matrix).operand(), 0UL, (*matrix).columns(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given column-major vector expansion operation.
// \ingroup row
//
// \param matrix The constant vector expansion operation.
// \param args The runtime row arguments.
// \return View on the specified row of the expansion operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given column-major
// vector expansion operation.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename... RRAs  // Runtime row arguments
        , EnableIf_t< !IsRowMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) row( const VecExpandExpr<MT,CEAs...>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs...> rd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   using ET = ElementType_t< MatrixType_t<MT> >;

   return UniformVector<ET,rowVector>( (*matrix).columns(), (*matrix).operand()[rd.row()] );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given matrix repeat operation.
// \ingroup row
//
// \param matrix The constant matrix repeat operation.
// \param args The runtime row arguments.
// \return View on the specified row of the repeat operation.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given row-major
// matrix repeat operation.
*/
template< size_t... CRAs1     // Compile time row arguments
        , typename MT         // Matrix base type of the expression
        , size_t... CRAs2     // Compile time repeater arguments
        , typename... RRAs >  // Runtime row arguments
inline decltype(auto) row( const MatRepeatExpr<MT,CRAs2...>& matrix, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const RowData<CRAs1...> rd( args... );

   if( isChecked( args... ) ) {
      if( (*matrix).rows() <= rd.row() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid row access index" );
      }
   }

   return repeat( row( (*matrix).operand()
                     , rd.row() % (*matrix).operand().rows()
                     , unchecked )
                , (*matrix).template repetitions<1UL>() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROW OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given row.
// \ingroup row
//
// \param row The row to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Row<MT,SO,DF,SF,CRAs...>& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary row.
// \ingroup row
//
// \param row The temporary row to be resetted.
// \return void
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void reset( Row<MT,SO,DF,SF,CRAs...>&& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given row.
// \ingroup row
//
// \param row The row to be cleared.
// \return void
//
// Clearing a row is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Row<MT,SO,DF,SF,CRAs...>& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary row.
// \ingroup row
//
// \param row The temporary row to be cleared.
// \return void
//
// Clearing a row is equivalent to resetting it via the reset() function.
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline void clear( Row<MT,SO,DF,SF,CRAs...>&& row )
{
   row.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense row is in default state.
// \ingroup row
//
// \param row The dense row to be tested for its default state.
// \return \a true in case the given dense row is component-wise zero, \a false otherwise.
//
// This function checks whether the dense row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all row elements are 0 and \a false in case any row element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
inline bool isDefault( const Row<MT,SO,true,SF,CRAs...>& row )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<row.size(); ++i )
      if( !isDefault<RF>( row[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse row is in default state.
// \ingroup row
//
// \param row The sparse row to be tested for its default state.
// \return \a true in case the given row is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse row is in default state. For instance, in case the
// row is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all row elements are 0 and \a false in case any row element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( row( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
inline bool isDefault( const Row<MT,SO,false,SF,CRAs...>& row )
{
   using blaze::isDefault;

   for( const auto& element : row )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given row are intact.
// \ingroup row
//
// \param row The row to be tested.
// \return \a true in case the given row's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the row are intact, i.e. if its state is valid.
// In case the invariants are intact, the function returns \a true, else it will return \a false.
// The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( row( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT       // Type of the matrix
        , bool SO           // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , size_t... CRAs >  // Compile time row arguments
inline bool isIntact( const Row<MT,SO,DF,SF,CRAs...>& row ) noexcept
{
   return ( row.row() < row.operand().rows() &&
            isIntact( row.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given rows represent the same observable state.
// \ingroup row
//
// \param a The first row to be tested for its state.
// \param b The second row to be tested for its state.
// \return \a true in case the two rows share a state, \a false otherwise.
//
// This overload of the isSame() function tests if the two given rows refer to exactly the same
// range of the same matrix. In case both rows represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1       // Type of the matrix of the left-hand side row
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF1           // Symmetry flag of the left-hand side row
        , size_t... CRAs1    // Compile time row arguments of the left-hand side row
        , typename MT2       // Type of the matrix of the right-hand side row
        , bool SF2           // Symmetry flag of the right-hand side row
        , size_t... CRAs2 >  // Compile time row arguments of the right-hand side row
inline bool isSame( const Row<MT1,SO,DF,SF1,CRAs1...>& a,
                    const Row<MT2,SO,DF,SF2,CRAs2...>& b ) noexcept
{
   return ( isSame( a.operand(), b.operand() ) && ( a.row() == b.row() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool trySet( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return trySet( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return trySet( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryAdd( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryAdd( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryAdd( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool trySub( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return trySub( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return trySub( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryMult( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryMult( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryMult( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryDiv( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryDiv( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryDiv( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs >  // Compile time row arguments
inline bool tryShift( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, int count )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryShift( row.operand(), row.row(), index, count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs >  // Compile time row arguments
BLAZE_ALWAYS_INLINE bool
   tryShift( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, int count )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryShift( row.operand(), row.row(), index, 1UL, size, count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryBitand( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryBitand( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitand( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryBitand( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryBitor( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryBitor( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryBitor( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
inline bool tryBitxor( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < row.size(), "Invalid vector access index" );

   return tryBitxor( row.operand(), row.row(), index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a row.
// \ingroup row
//
// \param row The target row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename ET >   // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const Row<MT,SO,DF,SF,CRAs...>& row, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*row).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*row).size(), "Invalid range size" );

   return tryBitxor( row.operand(), row.row(), index, 1UL, size, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                       const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryAddAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAddAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool trySubAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return trySubAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryMultAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                           const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryMultAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryDivAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryDivAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryShiftAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                            const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryShiftAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitandAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                             const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitandAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitorAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                            const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitorAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to a row.
// \ingroup row
//
// \param lhs The target left-hand side row.
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
        , size_t... CRAs  // Compile time row arguments
        , typename VT >   // Type of the right-hand side vector
inline bool tryBitxorAssign( const Row<MT,SO,DF,SF,CRAs...>& lhs,
                             const Vector<VT,true>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitxorAssign( lhs.operand(), *rhs, lhs.row(), index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row.
// \ingroup row
//
// \param r The row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given row. It returns a row
// object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Row index
inline decltype(auto) derestrict( Row<MT,SO,DF,SF,I>& r )
{
   return row<I>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row.
// \ingroup row
//
// \param r The temporary row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row. It
// returns a row object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF      // Symmetry flag
        , size_t I >   // Row index
inline decltype(auto) derestrict( Row<MT,SO,DF,SF,I>&& r )
{
   return row<I>( derestrict( r.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given row.
// \ingroup row
//
// \param r The row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given row. It returns a row
// object that does provide the same interface but does not have any restrictions on the data
// access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Row<MT,SO,DF,SF>& r )
{
   return row( derestrict( r.operand() ), r.row(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary row.
// \ingroup row
//
// \param r The temporary row to be derestricted.
// \return Row without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary row. It
// returns a row object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool SO      // Storage order
        , bool DF      // Density flag
        , bool SF >    // Symmetry flag
inline decltype(auto) derestrict( Row<MT,SO,DF,SF>&& r )
{
   return row( derestrict( r.operand() ), r.row(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given row.
// \ingroup row
//
// \param r The given row.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given row.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
inline decltype(auto) unview( Row<MT,SO,DF,SF,CRAs...>& r )
{
   return r.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given constant row.
// \ingroup row
//
// \param r The given constant row.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given constant row.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT        // Type of the matrix
        , bool SO            // Storage order
        , bool DF            // Density flag
        , bool SF            // Symmetry flag
        , size_t... CRAs >   // Compile time row arguments
inline decltype(auto) unview( const Row<MT,SO,DF,SF,CRAs...>& r )
{
   return r.operand();
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
struct Size< Row<MT,SO,DF,SF,CRAs...>, 0UL >
   : public Size<MT,1UL>
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
struct MaxSize< Row<MT,SO,DF,SF,CRAs...>, 0UL >
   : public MaxSize<MT,1UL>
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
template< typename MT, bool SO, bool DF, bool SF, size_t... CRAs >
struct IsRestricted< Row<MT,SO,DF,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct HasConstDataAccess< Row<MT,SO,true,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct HasMutableDataAccess< Row<MT,SO,true,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct IsAligned< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< IsAligned_v<MT> && ( IsRowMajorMatrix_v<MT> || IsSymmetric_v<MT> ) >
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISCONTIGUOUS SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SF, size_t... CRAs >
struct IsContiguous< Row<MT,true,true,SF,CRAs...> >
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
template< typename MT, bool SO, bool SF, size_t... CRAs >
struct IsPadded< Row<MT,SO,true,SF,CRAs...> >
   : public BoolConstant< IsPadded_v<MT> && ( IsRowMajorMatrix_v<MT> || IsSymmetric_v<MT> ) >
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
template< typename MT, bool DF, size_t... CRAs >
struct IsOpposedView< Row<MT,false,DF,false,CRAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
