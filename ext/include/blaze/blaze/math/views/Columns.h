//=================================================================================================
/*!
//  \file blaze/math/views/Columns.h
//  \brief Header file for the implementation of the Columns view
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <utility>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
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
#include <blaze/math/expressions/MatReduceExpr.h>
#include <blaze/math/expressions/MatRepeatExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/TVecMatMultExpr.h>
#include <blaze/math/expressions/VecExpandExpr.h>
#include <blaze/math/expressions/VecTVecMapExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsColumnMajorMatrix.h>
#include <blaze/math/typetraits/IsColumns.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Forward.h>
#include <blaze/math/views/column/ColumnData.h>
#include <blaze/math/views/columns/BaseTemplate.h>
#include <blaze/math/views/columns/Dense.h>
#include <blaze/math/views/columns/Sparse.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegerSequence.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsPointer.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns<4UL,2UL>( S );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns<1UL,3UL>( D, unchecked );
   auto columns2 = columns<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_< MT, index_sequence<I,Is...> >;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given constant matrix.
// \ingroup columns
//
// \param matrix The constant matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   const blaze::DynamicMatrix<double,columnMajor> D( ... );
   const blaze::CompressedMatrix<double,columnMajor> S( ... );

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns<1UL,3UL>( D );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns<4UL,2UL>( S );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns<1UL,3UL>( D, unchecked );
   auto columns2 = columns<4UL,2UL>( S, unchecked );
   \endcode
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const Matrix<MT,SO>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Columns_< const MT, index_sequence<I,Is...> >;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given temporary matrix.
// \ingroup columns
//
// \param matrix The temporary matrix containing the columns.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given temporary
// matrix. In case any column is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>&& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_< MT, index_sequence<I,Is...> >;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1.data(), indices1.size(), unchecked );
   auto columns2 = columns( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT>;
   return ReturnType( *matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given constant matrix.
// \ingroup columns
//
// \param matrix The constant matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1.data(), indices1.size() );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2.data(), indices2.size() );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1.data(), indices1.size(), unchecked );
   auto columns2 = columns( S, indices2.data(), indices2.size(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const Matrix<MT,SO>& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Columns_<const MT>;
   return ReturnType( *matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given temporary matrix.
// \ingroup columns
//
// \param matrix The temporary matrix containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given temporary
// matrix. In case any column is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>&& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT>;
   return ReturnType( *matrix, indices, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns( D, []( size_t i ){ return 2UL*i + 1UL; }, 2UL );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns( S, []( size_t i ){ return 4UL - 2UL*i; }, 2UL );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, []( size_t i ){ return 2UL*i + 1UL; }, 2UL, unchecked );
   auto columns2 = columns( S, []( size_t i ){ return 4UL - 2UL*i; }, 2UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT,P>;
   return ReturnType( *matrix, p, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given constant matrix.
// \ingroup columns
//
// \param matrix The constant matrix containing the columns.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given constant
// matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns( D, []( size_t i ){ return 2UL*i + 1UL; }, 2UL );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, []( size_t i ){ return 4UL - 2UL*i; }, 2UL );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, []( size_t i ){ return 2UL*i + 1UL; }, 2UL, unchecked );
   auto columns2 = columns( S, []( size_t i ){ return 4UL - 2UL*i; }, 2UL, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const Matrix<MT,SO>& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Columns_<const MT,P>;
   return ReturnType( *matrix, p, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given temporary matrix.
// \ingroup columns
//
// \param matrix The temporary matrix containing the columns.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args Optional arguments.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given temporary
// matrix. In case any column is not properly specified (i.e. if any specified index is greater
// than or equal to the total number of columns in the given matrix) a \a std::invalid_argument
// exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( Matrix<MT,SO>&& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Columns_<MT,P>;
   return ReturnType( *matrix, p, n, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The sequence of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;
   using blaze::index_sequence;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns( D, index_sequence<1UL,3UL>() );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns( S, index_sequence<4UL,2UL>() );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, index_sequence<1UL,3UL>(), unchecked );
   auto columns2 = columns( S, index_sequence<4UL,2UL>(), unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , size_t... Is        // Column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, index_sequence<Is...> indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   MAYBE_UNUSED( indices );

   return columns<Is...>( std::forward<MT>( matrix ), args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The list of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   auto columns1 = columns( D, { 1UL, 3UL } );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   auto columns2 = columns( S, { 4UL, 2UL } );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, { 1UL, 3UL }, unchecked );
   auto columns2 = columns( S, { 4UL, 2UL }, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, initializer_list<T> indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.begin(), indices.size(), args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The array of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::array<size_t,2UL> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1 );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::array<size_t,2UL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2 );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1, unchecked );
   auto columns2 = columns( S, indices2, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , size_t N            // Number of indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const std::array<T,N>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), N, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The vector of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   const std::vector<size_t> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1 );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   const std::vector<size_t> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2 );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1, unchecked );
   auto columns2 = columns( S, indices2, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const std::vector<T>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param indices The vector of column indices.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.

   \code
   using blaze::columnMajor;

   blaze::DynamicMatrix<double,columnMajor> D;
   blaze::CompressedMatrix<double,columnMajor> S;
   // ... Resizing and initialization

   // Creating a view on the 1st and 3rd column of the dense matrix D
   blaze::SmallArray<size_t,2UL> indices1{ 1UL, 3UL };
   auto columns1 = columns( D, indices1 );

   // Creating a view on the 4th and 2nd column of the sparse matrix S
   blaze::SmallArray<size_t,2UL> indices2{ 4UL, 2UL };
   auto columns2 = columns( S, indices2 );
   \endcode

// By default, the provided column indices are checked at runtime. In case any column is not
// properly specified (i.e. if any specified index is greater than or equal to the total number
// of columns in the given matrix) a \a std::invalid_argument exception is thrown. The checks
// can be skipped by providing the optional \a blaze::unchecked argument.

   \code
   auto columns1 = columns( D, indices1, unchecked );
   auto columns2 = columns( S, indices2, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , typename T          // Type of the column indices
        , size_t N            // Number of preallocated elements
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const SmallArray<T,N>& indices, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), indices.data(), indices.size(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns of the given matrix.
// \ingroup columns
//
// \param matrix The matrix containing the columns.
// \param pair The pair of arguments for the element selection.
// \param args Optional arguments.
// \return View on the specified columns of the matrix.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing a selection of columns of the given matrix.
// In case any column is not properly specified (i.e. if any specified index is greater than or
// equal to the total number of columns in the given matrix) a \a std::invalid_argument exception
// is thrown.
*/
template< typename MT         // Type of the matrix
        , typename T1         // First type of the pair of arguments
        , typename T2         // Second type of the pair of arguments
        , typename... RRAs >  // Optional arguments
inline decltype(auto) columns( MT&& matrix, const std::pair<T1,T2>& pair, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( std::forward<MT>( matrix ), pair.first, pair.second, args... );
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
/*!\brief Creating a view on a selection of columns on the given matrix/matrix addition.
// \ingroup columns
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the addition.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix addition.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatMatAddExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).leftOperand(), args... ) +
          columns<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/matrix subtraction.
// \ingroup columns
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the subtraction.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix subtraction.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatMatSubExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).leftOperand(), args... ) -
          columns<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given Schur product.
// \ingroup columns
//
// \param matrix The constant Schur product.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the Schur product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given Schur product.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const SchurExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).leftOperand(), args... ) %
          columns<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/matrix multiplication.
// \ingroup columns
//
// \param matrix The constant matrix/matrix multiplication.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/matrix multiplication.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatMatMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return (*matrix).leftOperand() * columns<CCAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given Kronecker product.
// \ingroup columns
//
// \param matrix The constant Kronecker product.
// \param args Optional arguments.
// \return View on the specified selection of columns on the Kronecker product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given Kronecker product.
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Matrix base type of the expression
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const MatMatKronExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   decltype(auto) lhs( (*matrix).leftOperand()  );
   decltype(auto) rhs( (*matrix).rightOperand() );

   const size_t M( rhs.rows()    );
   const size_t N( rhs.columns() );

   const auto lhsRows( [M]( size_t i ){ return i / M; } );
   const auto rhsRows( [M]( size_t i ){ return i % M; } );

   const auto lhsColumns( [N]( size_t i ) {
      constexpr size_t indices[] = { I, Is... };
      return indices[i] / N;
   } );

   const auto rhsColumns( [N]( size_t i ) {
      constexpr size_t indices[] = { I, Is... };
      return indices[i] % N;
   } );

   return columns( rows( lhs, lhsRows, (*matrix).rows(), args... ), lhsColumns, sizeof...(Is)+1UL, args... ) %
          columns( rows( rhs, rhsRows, (*matrix).rows(), args... ), rhsColumns, sizeof...(Is)+1UL, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given Kronecker product.
// \ingroup columns
//
// \param matrix The constant Kronecker product.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified selection of columns on the Kronecker product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given Kronecker product.
*/
template< typename MT         // Matrix base type of the expression
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const MatMatKronExpr<MT>& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   decltype(auto) lhs( (*matrix).leftOperand()  );
   decltype(auto) rhs( (*matrix).rightOperand() );

   const size_t M( rhs.rows()    );
   const size_t N( rhs.columns() );

   const auto lhsRows( [M]( size_t i ){ return i / M; } );
   const auto rhsRows( [M]( size_t i ){ return i % M; } );

   SmallArray<size_t,128UL> lhsColumns;
   lhsColumns.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      lhsColumns.pushBack( indices[i] / N );
   }

   SmallArray<size_t,128UL> rhsColumns;
   rhsColumns.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      rhsColumns.pushBack( indices[i] % N );
   }

   return columns( rows( lhs, lhsRows, (*matrix).rows(), args... ), lhsColumns, n, args... ) %
          columns( rows( rhs, rhsRows, (*matrix).rows(), args... ), rhsColumns, n, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given Kronecker product.
// \ingroup columns
//
// \param matrix The constant Kronecker product.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified selection of columns on the Kronecker product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given Kronecker product.
*/
template< typename MT         // Matrix base type of the expression
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional arguments
inline decltype(auto) columns( const MatMatKronExpr<MT>& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   decltype(auto) lhs( (*matrix).leftOperand()  );
   decltype(auto) rhs( (*matrix).rightOperand() );

   const size_t M( rhs.rows()    );
   const size_t N( rhs.columns() );

   const auto lhsRows( [M]( size_t i ){ return i / M; } );
   const auto rhsRows( [M]( size_t i ){ return i % M; } );

   const auto lhsColumns( [p,N]( size_t i ) { return p(i) / N; } );
   const auto rhsColumns( [p,N]( size_t i ) { return p(i) % N; } );

   return columns( rows( lhs, lhsRows, (*matrix).rows(), args... ), lhsColumns, n, args... ) %
          columns( rows( rhs, rhsRows, (*matrix).rows(), args... ), rhsColumns, n, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given outer product.
// \ingroup columns
//
// \param matrix The constant outer product.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the outer product.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given outer product.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const VecTVecMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return (*matrix).leftOperand() * elements<CCAs...>( (*matrix).rightOperand(), args... );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/scalar multiplication.
// \ingroup columns
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the multiplication.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/scalar multiplication.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatScalarMultExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).leftOperand(), args... ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix/scalar division.
// \ingroup columns
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the division.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix/scalar division.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatScalarDivExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).leftOperand(), args... ) / (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given unary matrix map operation.
// \ingroup columns
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the unary map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given unary matrix map operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( columns<CCAs...>( (*matrix).operand(), args... ), (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given binary matrix map operation.
// \ingroup columns
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the binary map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given binary matrix map operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatMatMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( columns<CCAs...>( (*matrix).leftOperand(), args... ),
               columns<CCAs...>( (*matrix).rightOperand(), args... ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given outer map operation.
// \ingroup columns
//
// \param matrix The constant outer map operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the outer map operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given outer map operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const VecTVecMapExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return map( (*matrix).leftOperand(),
                  elements<CCAs...>( (*matrix).rightOperand(), args... ), (*matrix).operation() );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix evaluation operation.
// \ingroup columns
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the evaluation operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix evaluation operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatEvalExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( columns<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix serialization operation.
// \ingroup columns
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the serialization operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix serialization operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatSerialExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( columns<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix no-alias operation.
// \ingroup columns
//
// \param matrix The constant matrix no-alias operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the no-alias operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix no-alias operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatNoAliasExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( columns<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix no-SIMD operation.
// \ingroup columns
//
// \param matrix The constant matrix no-SIMD operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the no-SIMD operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix no-SIMD operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatNoSIMDExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( columns<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix declaration operation.
// \ingroup columns
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the declaration operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix declaration operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const DeclExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns<CCAs...>( (*matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix transpose operation.
// \ingroup columns
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns on the transpose operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix transpose operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) >* = nullptr >
inline decltype(auto) columns( const MatTransExpr<MT>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return trans( rows<CCAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given column-major vector expansion
//        operation.
// \ingroup columns
//
// \param matrix The constant vector expansion operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns of the
// given column-major vector expansion operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) > 0UL ) &&
                      IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) columns( const VecExpandExpr<MT,CEAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      constexpr size_t indices[] = { CCAs... };
      for( size_t i=0UL; i<sizeof...(CCAs); ++i ) {
         if( (*matrix).columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   return expand< sizeof...( CCAs ) >( (*matrix).operand() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given column-major vector expansion
//        operation.
// \ingroup columns
//
// \param matrix The constant vector expansion operation.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The runtime column arguments.
// \return View on the specified selection of columns of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns of the
// given column-major vector expansion operation.
*/
template< typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename T        // Type of the column indices or index producer
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) columns( const VecExpandExpr<MT,CEAs...>& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( (*matrix).columns() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   return expand( (*matrix).operand(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given column-major vector expansion
//        operation.
// \ingroup columns
//
// \param matrix The constant vector expansion operation.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The runtime column arguments.
// \return View on the specified selection of columns of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns of the
// given column-major vector expansion operation.
*/
template< typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename P        // Type of the index producer
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) columns( const VecExpandExpr<MT,CEAs...>& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( (*matrix).columns() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   return expand( (*matrix).operand(), n );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given row-major vector expansion operation.
// \ingroup columns
//
// \param matrix The constant vector expansion operation.
// \param args The runtime column arguments.
// \return View on the specified selection of columns of the expansion operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns of the
// given row-major vector expansion operation.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Matrix base type of the expression
        , size_t... CEAs    // Compile time expansion arguments
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< ( sizeof...( CCAs ) + sizeof...( RCAs ) > 0UL ) &&
                      !IsColumnMajorMatrix_v<MT> >* = nullptr >
inline decltype(auto) columns( const VecExpandExpr<MT,CEAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return expand<CEAs...>( elements<CCAs...>( (*matrix).operand(), args... ), (*matrix).expansion() );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix repeat operation.
// \ingroup columns
//
// \param matrix The constant matrix repeat operation.
// \param args Optional arguments.
// \return View on the specified selection of columns on the matrix repeat operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix repeat operation.
*/
template< size_t I            // First column index
        , size_t... Is        // Remaining column indices
        , typename MT         // Matrix base type of the expression
        , size_t... CRAs      // Compile time repeater arguments
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const MatRepeatExpr<MT,CRAs...>& matrix, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( (*matrix).columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   auto lambda = [columns=(*matrix).operand().columns()]( size_t i ) {
      constexpr size_t indices[] = { I, Is... };
      return indices[i] % columns;
   };

   return repeat( columns( (*matrix).operand(), std::move(lambda), sizeof...(Is)+1UL, unchecked )
                , (*matrix).template repetitions<0UL>(), 1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix repeat operation.
// \ingroup columns
//
// \param matrix The constant matrix repeat operation.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified selection of columns on the matrix repeat operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix repeat operation.
*/
template< typename MT         // Matrix base type of the expression
        , size_t... CRAs      // Compile time repeater arguments
        , typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const MatRepeatExpr<MT,CRAs...>& matrix, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( (*matrix).columns() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices( indices, indices+n );

   for( size_t& index : newIndices ) {
      index = index % (*matrix).operand().columns();
   }

   return repeat( columns( (*matrix).operand(), newIndices, unchecked )
                , (*matrix).template repetitions<0UL>(), 1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of columns on the given matrix repeat operation.
// \ingroup columns
//
// \param matrix The constant matrix repeat operation.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args Optional arguments.
// \return View on the specified selection of columns on the matrix repeat operation.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified selection of columns on the
// given matrix repeat operation.
*/
template< typename MT         // Matrix base type of the expression
        , size_t... CRAs      // Compile time repeater arguments
        , typename P          // Type of the index producer
        , typename... RCAs >  // Optional column arguments
inline decltype(auto) columns( const MatRepeatExpr<MT,CRAs...>& matrix, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( (*matrix).columns() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   auto lambda = [columns=(*matrix).operand().columns(),p]( size_t i ) {
      return p(i) % columns;
   };

   return repeat( columns( (*matrix).operand(), std::move(lambda), n, unchecked )
                , (*matrix).template repetitions<0UL>(), 1UL );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMNS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
//
// This function returns an expression representing the specified columns of the given column
// selection.
*/
template< size_t I          // First required column index
        , size_t... Is      // Remaining required column indices
        , typename MT       // Type of the matrix
        , typename... RCAs  // Optional column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > &&
                      RemoveReference_t<MT>::compileTimeArgs >* = nullptr >
inline decltype(auto) columns( MT&& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return columns( c.operand(), subsequence<I,Is...>( RemoveReference_t<MT>::idces() ), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified columns of the given column
// selection.
*/
template< size_t I          // First required column index
        , size_t... Is      // Remaining required column indices
        , typename MT       // Type of the matrix
        , typename... RCAs  // Optional column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > &&
                      !RemoveReference_t<MT>::compileTimeArgs >* = nullptr >
inline decltype(auto) columns( MT&& c, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      constexpr size_t indices[] = { I, Is... };
      for( size_t i=0UL; i<sizeof...(Is)+1UL; ++i ) {
         if( c.columns() <= indices[i] ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   return columns( c.operand(), { c.idx(I), c.idx(Is)... }, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< typename MT       // Type of the matrix
        , typename T        // Type of the column indices
        , typename... RCAs  // Optional column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > >* = nullptr >
inline decltype(auto) columns( MT&& c, T* indices, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( c.columns() <= size_t( indices[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( c.idx( indices[i] ) );
   }

   return columns( c.operand(), newIndices.data(), newIndices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on specific columns of the given column selection.
// \ingroup columns
//
// \param c The selection of columns containing the columns.
// \param p Callable producing the indices.
// \param n The total number of indices.
// \param args The optional column arguments.
// \return View on the specified columns of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< typename MT       // Type of the matrix
        , typename P        // Type of the index producer
        , typename... RCAs  // Optional column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > && !IsPointer_v<P> >* = nullptr >
inline decltype(auto) columns( MT&& c, P p, size_t n, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<n; ++i ) {
         if( c.columns() <= size_t( p(i) ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
         }
      }
   }

   SmallArray<size_t,128UL> newIndices;
   newIndices.reserve( n );

   for( size_t i=0UL; i<n; ++i ) {
      newIndices.pushBack( c.idx( p(i) ) );
   }

   return columns( c.operand(), newIndices.data(), newIndices.size(), unchecked );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ELEMENTS)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given vector/matrix multiplication.
// \ingroup columns
//
// \param vector The constant vector/matrix multiplication.
// \param args The runtime element arguments.
// \return View on the specified elements of the multiplication.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified elements of the given
// vector/matrix multiplication.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const TVecMatMultExpr<VT>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return (*vector).leftOperand() * columns<CEAs...>( (*vector).rightOperand(), args... );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a selection of elements of the given column-wise matrix reduction
//        operation.
// \ingroup columns
//
// \param vector The constant column-wise matrix reduction operation.
// \param args The runtime element arguments.
// \return View on the specified elements of the multiplication.
// \exception std::invalid_argument Invalid element access index.
//
// This function returns an expression representing the specified elements of the given
// column-wise matrix reduction operation.
*/
template< size_t... CEAs      // Compile time element arguments
        , typename VT         // Vector base type of the expression
        , typename... REAs >  // Runtime element arguments
inline decltype(auto) elements( const MatReduceExpr<VT,columnwise>& vector, REAs... args )
{
   BLAZE_FUNCTION_TRACE;

   try {
      return reduce<columnwise>( columns<CEAs...>( (*vector).operand(), args... ), (*vector).operation() );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (ROW)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific row of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the row.
// \param args The runtime row arguments.
// \return View on the specified row of the column selection.
// \exception std::invalid_argument Invalid row access index.
//
// This function returns an expression representing the specified row of the given column selection.
*/
template< size_t... CRAs    // Compile time row arguments
        , typename MT       // Type of the matrix
        , typename... RRAs  // Runtime row arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > >* = nullptr >
inline decltype(auto) row( MT&& columns, RRAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return elements( row<CRAs...>( columns.operand(), args... ), columns.idces() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (COLUMN)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the column.
// \param args The optional column arguments.
// \return View on the specified column of the column selection.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< size_t I          // Column index
        , typename MT       // Type of the matrix
        , typename... RCAs  // Optional column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > &&
                      RemoveReference_t<MT>::compileTimeArgs >* = nullptr >
inline decltype(auto) column( MT&& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return column< RemoveReference_t<MT>::idx(I) >( columns.operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific column of the given column selection.
// \ingroup columns
//
// \param columns The selection of columns containing the column.
// \param args The runtime column arguments.
// \return View on the specified column of the column selection.
// \exception std::invalid_argument Invalid column access index.
//
// This function returns an expression representing the specified column of the given column
// selection.
*/
template< size_t... CCAs    // Compile time column arguments
        , typename MT       // Type of the matrix
        , typename... RCAs  // Runtime column arguments
        , EnableIf_t< IsColumns_v< RemoveReference_t<MT> > &&
                      ( sizeof...( CCAs ) == 0UL || !RemoveReference_t<MT>::compileTimeArgs ) >* = nullptr >
inline decltype(auto) column( MT&& columns, RCAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const ColumnData<CCAs...> cd( args... );

   if( isChecked( args... ) ) {
      if( columns.columns() <= cd.column() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid column access index" );
      }
   }

   return column( columns.operand(), columns.idx( cd.column() ), args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given column selection.
// \ingroup columns
//
// \param columns The column selection to be resetted.
// \return void
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary column selection.
// \ingroup columns
//
// \param columns The temporary column selection to be resetted.
// \return void
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>&& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column of the given column selection.
// \ingroup columns
//
// \param columns The column selection to be resetted.
// \param i The index of the column to be resetted.
// \return void
//
// This function resets the values in the specified column of the given column selection to their
// default value. Note that the capacity of the column remains unchanged.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void reset( Columns<MT,SO,DF,SF,CCAs...>& columns, size_t i )
{
   columns.reset( i );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given column selection.
// \ingroup columns
//
// \param columns The column selection to be cleared.
// \return void
//
// Clearing a column selection is equivalent to resetting it via the reset() function.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void clear( Columns<MT,SO,DF,SF,CCAs...>& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary column selection.
// \ingroup columns
//
// \param columns The column selection to be cleared.
// \return void
//
// Clearing a column selection is equivalent to resetting it via the reset() function.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline void clear( Columns<MT,SO,DF,SF,CCAs...>&& columns )
{
   columns.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense column selection is in default state.
// \ingroup columns
//
// \param columns The dense column selection to be tested for its default state.
// \return \a true in case the given dense column selection is component-wise zero, \a false otherwise.
//
// This function checks whether the dense column selection is in default state. For instance, in
// case the column selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all column elements are 0 and \a false in case any column
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF   // Relaxation flag
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool isDefault( const Columns<MT,SO,true,SF,CCAs...>& columns )
{
   using blaze::isDefault;

   if( SO == false ) {
      for( size_t i=0UL; i<columns.rows(); ++i )
         for( size_t j=0UL; j<columns.columns(); ++j )
            if( !isDefault<RF>( columns(i,j) ) )
               return false;
   }
   else {
      for( size_t j=0UL; j<columns.columns(); ++j )
         for( size_t i=0UL; i<columns.rows(); ++i )
            if( !isDefault<RF>( columns(i,j) ) )
               return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse column selection is in default state.
// \ingroup columns
//
// \param columns The sparse column selection to be tested for its default state.
// \return \a true in case the given sparse column selection is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse column selection is in default state. For instance, in
// case the column selection is instantiated for a built-in integral or floating point data type,
// the function returns \a true in case all column elements are 0 and \a false in case any column
// element is not 0. The following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isDefault( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF   // Relaxation flag
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool isDefault( const Columns<MT,SO,false,SF,CCAs...>& columns )
{
   using blaze::isDefault;

   for( size_t j=0UL; j<columns.columns(); ++j ) {
      for( auto element=columns.cbegin(j); element!=columns.cend(j); ++element )
         if( !isDefault<RF>( element->value() ) ) return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given column selection are intact.
// \ingroup columns
//
// \param columns The column selection to be tested.
// \return \a true in case the given column selection's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the column selection are intact, i.e. if its
// state is valid. In case the invariants are intact, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<double,columnMajor> A;
   // ... Resizing and initialization
   if( isIntact( columns( A, { 2UL, 4UL, 6UL, 8UL } ) ) ) { ... }
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool isIntact( const Columns<MT,SO,DF,SF,CCAs...>& columns ) noexcept
{
   return ( columns.rows() == columns.operand().rows() &&
            columns.columns() <= columns.operand().columns() &&
            isIntact( columns.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given column selection and matrix represent the same observable state.
// \ingroup columns
//
// \param a The column selection to be tested for its state.
// \param b The matrix to be tested for its state.
// \return \a true in case the column selection and matrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to all columns
// of the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side column selection
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool isSame( const Columns<MT,SO1,DF,SF,CCAs...>& a, const Matrix<MT,SO2>& b ) noexcept
{
   if( !isSame( a.operand(), *b ) || ( a.rows() != (*b).rows() ) || ( a.columns() != (*b).columns() ) )
      return false;

   for( size_t j=0UL; j<a.columns(); ++j ) {
      if( a.idx(j) != j )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given matrix and column selection represent the same observable state.
// \ingroup columns
//
// \param a The matrix to be tested for its state.
// \param b The column selection to be tested for its state.
// \return \a true in case the matrix and column selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to all columns
// of the given matrix in ascending and consecutive order and by that represents the same observable
// state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side matrix
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , bool SO2 >        // Storage order of the right-hand side column selection
inline bool isSame( const Matrix<MT,SO1>& a, const Columns<MT,SO2,DF,SF,CCAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given column selection and submatrix represent the same observable state.
// \ingroup columns
//
// \param a The column selection to be tested for its state.
// \param b The submatrix to be tested for its state.
// \return \a true in case the column selection and submatrix share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to same columns
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT       // Type of the matrix
        , bool SO1          // Storage order of the left-hand side column selection
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , AlignmentFlag AF  // Alignment flag
        , bool SO2          // Storage order of the right-hand side submatrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool isSame( const Columns<MT,SO1,DF,SF,CCAs...>& a, const Submatrix<MT,AF,SO2,DF,CSAs...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || ( a.rows() != (*b).rows() ) || ( a.columns() != (*b).columns() ) )
      return false;

   for( size_t j=0UL; j<a.columns(); ++j ) {
      if( a.idx(j) != b.column()+j )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given submatrix and column selection represent the same observable state.
// \ingroup columns
//
// \param a The submatrix to be tested for its state.
// \param b The column selection to be tested for its state.
// \return \a true in case the submatrix and column selection share a state, \a false otherwise.
//
// This overload of the isSame function tests if the given column selection refers to same columns
// as the given submatrix in ascending and consecutive order and by that represents the same
// observable state. In this case, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT         // Type of the matrix
        , AlignmentFlag AF    // Alignment flag
        , bool SO1            // Storage order of the left-hand side submatrix
        , bool DF             // Density flag
        , size_t... CSAs      // Compile time submatrix arguments
        , bool SO2            // Storage order of the right-hand side column selection
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool isSame( const Submatrix<MT,AF,SO1,DF,CSAs...>& a, const Columns<MT,SO2,DF,SF,CCAs...>& b ) noexcept
{
   return isSame( b, a );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given column selections represent the same observable state.
// \ingroup columns
//
// \param a The first selection of columns to be tested for its state.
// \param b The second selection of columns to be tested for its state.
// \return \a true in case the two column selections share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given column selections refer to exactly
// the same range of the same matrix. In case both selections represent the same observable
// state, the function returns \a true, otherwise it returns \a false.
*/
template< typename MT          // Type of the matrix
        , bool SO              // Storage order
        , bool DF              // Density flag
        , bool SF              // Symmetry flag
        , typename... CCAs1    // Compile time column arguments of the left-hand side column selection
        , typename... CCAs2 >  // Compile time column arguments of the right-hand side column selection
inline bool isSame( const Columns<MT,SO,DF,SF,CCAs1...>& a,
                    const Columns<MT,SO,DF,SF,CCAs2...>& b ) noexcept
{
   if( !isSame( a.operand(), b.operand() ) || a.rows() != b.rows() || a.columns() != b.columns() )
      return false;

   for( size_t i=0UL; i<a.columns(); ++i ) {
      if( a.idx(i) != b.idx(i) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense column selection.
// \ingroup columns
//
// \param c The dense column selection to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given dense column selection by means of the specified matrix type
// or matrix inversion algorithm \c IF (see the InversionFlag documentation):

   \code
   invert<asLower>( A );     // Inversion of a lower triangular matrix
   invert<asUniUpper>( A );  // Inversion of an upper unitriangular matrix
   invert<byLU>( A );        // Inversion by means of an LU decomposition
   invert<byLLH>( A );       // Inversion by means of a Cholesky decomposition
   ...
   \endcode

// The matrix inversion fails if ...
//
//  - ... the given column selection is not a square matrix;
//  - ... the given column selection is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a c may already have been modified.
*/
template< InversionFlag IF    // Inversion algorithm
        , typename MT         // Type of the dense matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline auto invert( Columns<MT,SO,true,SF,CCAs...>& c )
   -> DisableIf_t< HasMutableDataAccess_v<MT> >
{
   using RT = ResultType_t< Columns<MT,SO,true,SF,CCAs...> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION  ( RT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( RT );

   RT tmp( c );
   invert<IF>( tmp );
   c = tmp;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be set.
// \param j The column index of the element to be set.
// \param value The value to be set to the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool trySet( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return trySet( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The value to be set to the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !trySet( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be added to the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryAdd( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryAdd( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The value to be added to the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryAdd( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The value to be subtracted from the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool trySub( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return trySub( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The value to be subtracted to the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !trySub( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The factor for the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryMult( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryMult( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The factor for the elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryMult( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The divisor for the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryDiv( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryDiv( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The divisor for the elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryDiv( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param count The number of bits to shift the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline bool tryShift( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, int count )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryShift( c.operand(), c.idx(i), j, count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param count The number of bits to shift the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
BLAZE_ALWAYS_INLINE bool
   tryShift( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, int count )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryShift( c.operand(), row, c.idx(j), m, 1UL, count ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The bit pattern to be used on the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryBitand( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryBitand( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitand( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryBitand( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The bit pattern to be used on the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryBitor( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryBitor( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryBitor( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param i The row index of the element to be modified.
// \param j The column index of the element to be modified.
// \param value The bit pattern to be used on the element.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
inline bool tryBitxor( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < c.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < c.columns(), "Invalid column access index" );

   return tryBitxor( c.operand(), c.idx(i), j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a column selection.
// \ingroup columns
//
// \param c The target column selection.
// \param row The index of the first row of the range to be modified.
// \param column The index of the first column of the range to be modified.
// \param m The number of rows of the range to be modified.
// \param n The number of columns of the range to be modified.
// \param value The bit pattern to be used on the range of elements.
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
        , typename... CCAs  // Compile time column arguments
        , typename ET >     // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const Columns<MT,SO,DF,SF,CCAs...>& c, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*c).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*c).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*c).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*c).columns(), "Invalid number of columns" );

   const size_t jend( column + n );

   for( size_t j=column; j<jend; ++j ) {
      if( !tryBitxor( c.operand(), row, c.idx(j), m, 1UL, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a column vector to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a row vector to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to the band of a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be assigned.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                       const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !trySet( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a matrix to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                       const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a column vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryAddAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be added.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryAdd( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryAddAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryAddAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return trySubAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be subtracted.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool trySubAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !trySub( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool trySubAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !trySubAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryMultAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a row vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector to be multiplied.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to the band
//        of a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector to be multiplied.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryMultAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                           const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryMult( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the Schur product assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix for the Schur product.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool trySchurAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryMultAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a column vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector divisor.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryDivAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector divisor.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector divisor.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryDivAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                          const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryDiv( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a column vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector of bits to shift.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryShiftAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryShiftAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector of bits to shift.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryShiftAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryShift( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a vector to the band of a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector of bits to shift.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryShiftAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryShift( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a matrix to a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix of bits to shift.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryShiftAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryShiftAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector for the bitwise AND operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitandAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryBitandAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector for the bitwise AND operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitandAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitand( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a vector to the band of
//        a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector for the bitwise AND operation.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryBitandAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitand( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix for the bitwise AND operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryBitandAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                             const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryBitandAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector for the bitwise OR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryBitorAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector for the bitwise OR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitor( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to the band of
//        a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector for the bitwise OR operation.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryBitorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                            const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitor( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix for the bitwise OR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryBitorAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryBitorAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a column vector to a
//        column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side column vector for the bitwise XOR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitxorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   return tryBitxorAssign( lhs.operand(), *rhs, row, lhs.idx( column ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a row vector to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side row vector for the bitwise XOR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT >     // Type of the right-hand side vector
inline bool tryBitxorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitxor( lhs.operand(), row, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to the band of
//        a column selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side vector for the bitwise XOR operation.
// \param band The index of the band the right-hand side vector is assigned to.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
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
        , typename... CCAs  // Compile time column arguments
        , typename VT       // Type of the right-hand side vector
        , bool TF >         // Transpose flag of the right-hand side vector
inline bool tryBitxorAssign( const Columns<MT,SO,DF,SF,CCAs...>& lhs,
                             const Vector<VT,TF>& rhs, ptrdiff_t band, size_t row, size_t column )
{
   MAYBE_UNUSED( band );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<(*rhs).size(); ++i ) {
      if( !tryBitxor( lhs.operand(), row+i, lhs.idx( column+i ), (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a matrix to a column
//        selection.
// \ingroup columns
//
// \param lhs The target left-hand side column selection.
// \param rhs The right-hand side matrix for the bitwise XOR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1      // Type of the matrix
        , bool SO1          // Storage order
        , bool DF           // Density flag
        , bool SF           // Symmetry flag
        , typename... CCAs  // Compile time column arguments
        , typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline bool tryBitxorAssign( const Columns<MT1,SO1,DF,SF,CCAs...>& lhs,
                             const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<(*rhs).columns(); ++j ) {
      if( !tryBitxorAssign( lhs.operand(), blaze::column( *rhs, j, unchecked ), row, lhs.idx( column+j ) ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given column selection.
// \ingroup columns
//
// \param c The column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given column selection.
// It returns a column selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF,CCAs...>& c )
{
   return columns( derestrict( c.operand() ), c.idces(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary column selection.
// \ingroup columns
//
// \param r The temporary column selection to be derestricted.
// \return Column selection without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary column selection.
// It returns a column selection that does provide the same interface but does not have any
// restrictions on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline decltype(auto) derestrict( Columns<MT,SO,DF,SF,CCAs...>&& c )
{
   return columns( derestrict( c.operand() ), c.idces(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given column selection.
// \ingroup columns
//
// \param c The given column selection.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given column selection.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline decltype(auto) unview( Columns<MT,SO,DF,SF,CCAs...>& c )
{
   return c.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given constant column selection.
// \ingroup columns
//
// \param c The given constant column selection.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given constant column
// selection.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool DF             // Density flag
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
inline decltype(auto) unview( const Columns<MT,SO,DF,SF,CCAs...>& c )
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
template< typename MT, bool SO, bool DF, bool SF, typename... CCAs >
struct Size< Columns<MT,SO,DF,SF,CCAs...>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is, typename... CCAs >
struct Size< Columns<MT,SO,DF,SF,index_sequence<I,Is...>,CCAs...>, 1UL >
   : public Ptrdiff_t<1UL+sizeof...(Is)>
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
template< typename MT, bool SO, bool DF, bool SF, typename... CCAs >
struct MaxSize< Columns<MT,SO,DF,SF,CCAs...>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO, bool DF, bool SF, size_t I, size_t... Is, typename... CCAs >
struct MaxSize< Columns<MT,SO,DF,SF,index_sequence<I,Is...>,CCAs...>, 1UL >
   : public Ptrdiff_t<1UL+sizeof...(Is)>
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
template< typename MT, bool SO, bool DF, bool SF, typename... CCAs >
struct IsRestricted< Columns<MT,SO,DF,SF,CCAs...> >
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
template< typename MT, bool SO, bool SF, typename... CCAs >
struct HasConstDataAccess< Columns<MT,SO,true,SF,CCAs...> >
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
template< typename MT, bool SO, bool SF, typename... CCAs >
struct HasMutableDataAccess< Columns<MT,SO,true,SF,CCAs...> >
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
template< typename MT, bool SO, bool SF, typename... CCAs >
struct IsAligned< Columns<MT,SO,true,SF,CCAs...> >
   : public IsAligned<MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
