//=================================================================================================
/*!
//  \file blaze/math/views/Band.h
//  \brief Header file for the implementation of the Band view
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

#ifndef _BLAZE_MATH_VIEWS_BAND_H_
#define _BLAZE_MATH_VIEWS_BAND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DeclExpr.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/MatEvalExpr.h>
#include <blaze/math/expressions/MatMapExpr.h>
#include <blaze/math/expressions/MatMatAddExpr.h>
#include <blaze/math/expressions/MatMatKronExpr.h>
#include <blaze/math/expressions/MatMatMapExpr.h>
#include <blaze/math/expressions/MatMatSubExpr.h>
#include <blaze/math/expressions/MatNoAliasExpr.h>
#include <blaze/math/expressions/MatNoSIMDExpr.h>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/expressions/MatScalarDivExpr.h>
#include <blaze/math/expressions/MatScalarMultExpr.h>
#include <blaze/math/expressions/MatSerialExpr.h>
#include <blaze/math/expressions/MatTransExpr.h>
#include <blaze/math/expressions/SchurExpr.h>
#include <blaze/math/expressions/VecExpandExpr.h>
#include <blaze/math/expressions/VecTVecMapExpr.h>
#include <blaze/math/expressions/VecTVecMultExpr.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsBand.h>
#include <blaze/math/typetraits/IsOpposedView.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSubmatrix.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/math/views/band/BandData.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/band/Dense.h>
#include <blaze/math/views/band/Sparse.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Forward.h>
#include <blaze/math/views/subvector/SubvectorData.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/mpl/Min.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the band.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band<1L>( D );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band<-1L>( S );
   \endcode

// By default, the provided band arguments are checked at runtime. // In case the band is not
// properly specified (i.e. if the specified index does not correspond to a valid band in the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.

   \code
   auto ub1 = band<1L>( D, unchecked );
   auto lb1 = band<-1L>( S, unchecked );
   \endcode
*/
template< ptrdiff_t I         // Band index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( Matrix<MT,SO>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Band_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the band.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band<1L>( D );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band<-1L>( S );
   \endcode

// By default, the provided band arguments are checked at runtime. // In case the band is not
// properly specified (i.e. if the specified index does not correspond to a valid band in the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.

   \code
   auto ub1 = band<1L>( D, unchecked );
   auto lb1 = band<-1L>( S, unchecked );
   \endcode
*/
template< ptrdiff_t I         // Band index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( const Matrix<MT,SO>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Band_<const MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the band.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given temporary
// matrix. In case the band is not properly specified (i.e. if the specified index does not
// correspond to a valid band in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< ptrdiff_t I         // Band index
        , typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( Matrix<MT,SO>&& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Band_<MT,I>;
   return ReturnType( *matrix, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the band.
// \param index The band index.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto ub1 = band( D, 1L );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto lb1 = band( S, -1L );
   \endcode

// By default, the provided band arguments are checked at runtime. // In case the band is not
// properly specified (i.e. if the specified index does not correspond to a valid band in the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.

   \code
   auto ub1 = band( D, 1L, unchecked );
   auto lb1 = band( S, -1L, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( Matrix<MT,SO>& matrix, ptrdiff_t index, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Band_<MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the band.
// \param index The band index.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given constant
// matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the upper secondary diagonal of the dense matrix D
   auto lb1 = band( D, 1L );

   // Creating a view on the lower secondary diagonal of the sparse matrix S
   auto ub1 = band( S, -1L );
   \endcode

// By default, the provided band arguments are checked at runtime. // In case the band is not
// properly specified (i.e. if the specified index does not correspond to a valid band in the
// given matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.

   \code
   auto ub1 = band( D, 1L, unchecked );
   auto lb1 = band( S, -1L, unchecked );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( const Matrix<MT,SO>& matrix, ptrdiff_t index, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = const Band_<const MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on a specific band of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the band.
// \param index The band index.
// \param args Optional band arguments.
// \return View on the specified band of the matrix.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given temporary
// matrix. In case the band is not properly specified (i.e. if the specified index does not
// correspond to a valid band in the given matrix) a \a std::invalid_argument exception is
// thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RBAs >  // Optional band arguments
inline decltype(auto) band( Matrix<MT,SO>&& matrix, ptrdiff_t index, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Band_<MT>;
   return ReturnType( *matrix, index, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given matrix.
// \ingroup band
//
// \param matrix The matrix containing the diagonal.
// \param args Optional diagonal arguments.
// \return View on the diagonal of the matrix.
//
// This function returns an expression representing the diagonal of the given matrix.

   \code
   using blaze::rowMajor;

   blaze::DynamicMatrix<double,rowMajor> D;
   blaze::CompressedMatrix<double,rowMajor> S;
   // ... Resizing and initialization

   // Creating a view on the diagonal of the dense matrix D
   auto diag1 = diagonal( D );

   // Creating a view on the diagonal of the sparse matrix S
   auto diag2 = diagonal( S );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RDAs >  // Optional diagonal arguments
inline decltype(auto) diagonal( Matrix<MT,SO>& matrix, RDAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( *matrix, unchecked, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given constant matrix.
// \ingroup band
//
// \param matrix The constant matrix containing the diagonal.
// \param args Optional diagonal arguments.
// \return View on the diagonal of the matrix.
//
// This function returns an expression representing the diagonal of the given constant matrix.

   \code
   using blaze::rowMajor;

   const blaze::DynamicMatrix<double,rowMajor> D( ... );
   const blaze::CompressedMatrix<double,rowMajor> S( ... );

   // Creating a view on the diagonal of the dense matrix D
   auto diag1 = diagonal( D );

   // Creating a view on the diagonal of the sparse matrix S
   auto diag2 = diagonal( S );
   \endcode
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RDAs >  // Optional diagonal arguments
inline decltype(auto) diagonal( const Matrix<MT,SO>& matrix, RDAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( *matrix, unchecked, args... );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Creating a view on the diagonal of the given temporary matrix.
// \ingroup band
//
// \param matrix The temporary matrix containing the diagonal.
// \param args Optional diagonal arguments.
// \return View on the diagonal of the matrix.
//
// This function returns an expression representing the diagonal of the given temporary matrix.
// In case the diagonal is not properly specified (i.e. in case the given matrix has zero rows
// or zero columns) a \a std::invalid_argument exception is thrown.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , typename... RDAs >  // Optional diagonal arguments
inline decltype(auto) diagonal( Matrix<MT,SO>&& matrix, RDAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<0L>( *matrix, unchecked, args... );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix addition.
// \ingroup band
//
// \param matrix The constant matrix/matrix addition.
// \param args The runtime band arguments.
// \return View on the specified band of the addition.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// addition.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatMatAddExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).leftOperand(), args... ) +
          band<CBAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/matrix subtraction.
// \ingroup band
//
// \param matrix The constant matrix/matrix subtraction.
// \param args The runtime band arguments.
// \return View on the specified band of the subtraction.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix/matrix
// subtraction.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatMatSubExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).leftOperand(), args... ) -
          band<CBAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given Schur product.
// \ingroup band
//
// \param matrix The constant Schur product.
// \param args The runtime band arguments.
// \return View on the specified band of the Schur product.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given Schur product.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const SchurExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).leftOperand(), args... ) *
          band<CBAs...>( (*matrix).rightOperand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given Kronecker product.
// \ingroup band
//
// \param matrix The constant Kronecker product.
// \param args The runtime band arguments.
// \return View on the specified band of the Kronecker product.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given Kronecker
// product.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatMatKronExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const BandData<CBAs...> bd( args... );

   if( isChecked( args... ) ) {
      if( ( bd.band() > 0L && bd.column() >= (*matrix).columns() ) ||
          ( bd.band() < 0L && bd.row() >= (*matrix).rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }

   const size_t row   ( bd.band() <  0L ? -bd.band() : 0UL );
   const size_t column( bd.band() >= 0L ?  bd.band() : 0UL );
   const size_t n     ( min( (*matrix).rows() - row, (*matrix).columns() - column ) );

   return diagonal( submatrix( *matrix, row, column, n, n, unchecked ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given outer product.
// \ingroup band
//
// \param matrix The constant Kronecker product.
// \param args The runtime band arguments.
// \return View on the specified band of the Kronecker product.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given outer product.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const VecTVecMultExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const BandData<CBAs...> bd( args... );

   if( isChecked( args... ) ) {
      if( ( bd.band() > 0L && bd.column() >= (*matrix).columns() ) ||
          ( bd.band() < 0L && bd.row() >= (*matrix).rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }

   decltype(auto) leftOperand ( (*matrix).leftOperand()  );
   decltype(auto) rightOperand( (*matrix).rightOperand() );

   const size_t row   ( bd.band() <  0L ? -bd.band() : 0UL );
   const size_t column( bd.band() >= 0L ?  bd.band() : 0UL );
   const size_t size  ( min( leftOperand.size() - row, rightOperand.size() - column ) );

   return transTo<defaultTransposeFlag>( subvector( leftOperand , row   , size, unchecked ) ) *
          transTo<defaultTransposeFlag>( subvector( rightOperand, column, size, unchecked ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar multiplication.
// \ingroup band
//
// \param matrix The constant matrix/scalar multiplication.
// \param args The runtime band arguments.
// \return View on the specified band of the multiplication.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// multiplication.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatScalarMultExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).leftOperand(), args... ) * (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix/scalar division.
// \ingroup band
//
// \param matrix The constant matrix/scalar division.
// \param args The runtime band arguments.
// \return View on the specified band of the division.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix/scalar
// division.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatScalarDivExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).leftOperand(), args... ) / (*matrix).rightOperand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given unary matrix map operation.
// \ingroup band
//
// \param matrix The constant unary matrix map operation.
// \param args The runtime band arguments.
// \return View on the specified band of the unary map operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given unary matrix
// map operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatMapExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( band<CBAs...>( (*matrix).operand(), args... ), (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given binary matrix map operation.
// \ingroup band
//
// \param matrix The constant binary matrix map operation.
// \param args The runtime band arguments.
// \return View on the specified band of the binary map operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given binary matrix
// map operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatMatMapExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return map( band<CBAs...>( (*matrix).leftOperand(), args... ),
               band<CBAs...>( (*matrix).rightOperand(), args... ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given outer map operation.
// \ingroup band
//
// \param matrix The constant outer map operation.
// \param args The runtime band arguments.
// \return View on the specified band of the outer map operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given outer map
// operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const VecTVecMapExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const BandData<CBAs...> bd( args... );

   if( isChecked( args... ) ) {
      if( ( bd.band() > 0L && bd.column() >= (*matrix).columns() ) ||
          ( bd.band() < 0L && bd.row() >= (*matrix).rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }

   decltype(auto) leftOperand ( (*matrix).leftOperand()  );
   decltype(auto) rightOperand( (*matrix).rightOperand() );

   const size_t row   ( bd.band() <  0L ? -bd.band() : 0UL );
   const size_t column( bd.band() >= 0L ?  bd.band() : 0UL );
   const size_t size  ( min( leftOperand.size() - row, rightOperand.size() - column ) );

   return map( transTo<defaultTransposeFlag>( subvector( leftOperand , row   , size, unchecked ) ),
               transTo<defaultTransposeFlag>( subvector( rightOperand, column, size, unchecked ) ),
               (*matrix).operation() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix evaluation operation.
// \ingroup band
//
// \param matrix The constant matrix evaluation operation.
// \param args The runtime band arguments.
// \return View on the specified band of the evaluation operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// evaluation operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatEvalExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return eval( band<CBAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix serialization operation.
// \ingroup band
//
// \param matrix The constant matrix serialization operation.
// \param args The runtime band arguments.
// \return View on the specified band of the serialization operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// serialization operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatSerialExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return serial( band<CBAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix no-alias operation.
// \ingroup band
//
// \param matrix The constant matrix no-alias operation.
// \param args The runtime band arguments.
// \return View on the specified band of the no-alias operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// no-alias operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatNoAliasExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return noalias( band<CBAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix no-SIMD operation.
// \ingroup band
//
// \param matrix The constant matrix no-SIMD operation.
// \param args The runtime band arguments.
// \return View on the specified band of the no-SIMD operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// no-SIMD operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatNoSIMDExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return nosimd( band<CBAs...>( (*matrix).operand(), args... ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix declaration operation.
// \ingroup band
//
// \param matrix The constant matrix declaration operation.
// \param args The runtime band arguments.
// \return View on the specified band of the declaration operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// declaration operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const DeclExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<CBAs...>( (*matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix transpose operation.
// \ingroup band
//
// \param matrix The constant matrix transpose operation.
// \param args The runtime band arguments.
// \return View on the specified band of the transpose operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// transpose operation.
*/
template< ptrdiff_t I         // Band index
        , typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatTransExpr<MT>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band<-I>( (*matrix).operand(), args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix transpose operation.
// \ingroup band
//
// \param matrix The constant matrix transpose operation.
// \param index The band index.
// \param args Optional band arguments.
// \return View on the specified band of the transpose operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// transpose operation.
*/
template< typename MT         // Type of the matrix
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatTransExpr<MT>& matrix, ptrdiff_t index, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   return band( (*matrix).operand(), -index, args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given vector expansion operation.
// \ingroup band
//
// \param matrix The constant vector expansion operation.
// \param args The runtime band arguments.
// \return View on the specified band of the expansion operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given vector
// expansion operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , size_t... CEAs      // Compile time expansion arguments
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const VecExpandExpr<MT,CEAs...>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const BandData<CBAs...> bd( args... );

   if( isChecked( args... ) ) {
      if( ( bd.band() > 0L && bd.column() >= (*matrix).columns() ) ||
          ( bd.band() < 0L && bd.row() >= (*matrix).rows() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid band access index" );
      }
   }

   using VT = VectorType_t< RemoveReference_t< decltype( (*matrix).operand() ) > >;

   constexpr bool TF( TransposeFlag_v<VT> );

   const size_t index( TF ? bd.column() : bd.row() );
   const size_t size ( min( (*matrix).rows() - bd.row(), (*matrix).columns() - bd.column() ) );

   return subvector( transTo<defaultTransposeFlag>( (*matrix).operand() ), index, size, unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific band of the given matrix repeat operation.
// \ingroup band
//
// \param matrix The constant matrix repeat operation.
// \param args The runtime band arguments.
// \return View on the specified band of the repeat operation.
// \exception std::invalid_argument Invalid band access index.
//
// This function returns an expression representing the specified band of the given matrix
// repeat operation.
*/
template< ptrdiff_t... CBAs   // Compile time band arguments
        , typename MT         // Type of the matrix
        , size_t... CRAs      // Compile time repeater arguments
        , typename... RBAs >  // Runtime band arguments
inline decltype(auto) band( const MatRepeatExpr<MT,CRAs...>& matrix, RBAs... args )
{
   BLAZE_FUNCTION_TRACE;

   using ReturnType = Band_< MatrixType_t<MT>, CBAs... >;
   return ReturnType( *matrix, args... );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL RESTRUCTURING FUNCTIONS (SUBVECTOR)
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given band.
// \ingroup band
//
// \param b The band containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the band.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given band.
*/
template< AlignmentFlag AF  // Alignment flag
        , size_t I          // Index of the first subvector element
        , size_t N          // Size of the subvector
        , typename VT       // Type of the vector
        , typename... RSAs  // Optional subvector arguments
        , EnableIf_t< IsBand_v< RemoveReference_t<VT> > &&
                      RemoveReference_t<VT>::compileTimeArgs >* = nullptr >
inline decltype(auto) subvector( VT&& b, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   constexpr ptrdiff_t I2 = RemoveReference_t<VT>::band();

   constexpr size_t row   ( ( I2 >= 0L ? 0UL : -I2 ) + I );
   constexpr size_t column( ( I2 >= 0L ?  I2 : 0UL ) + I );

   constexpr auto check( getCheck( args... ) );

   try {
      return diagonal( submatrix<AF,row,column,N,N>( b.operand(), check ), check );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creating a view on a specific subvector of the given band.
// \ingroup band
//
// \param b The band containing the subvector.
// \param args The optional subvector arguments.
// \return View on the specified subvector of the band.
// \exception std::invalid_argument Invalid subvector specification.
//
// This function returns an expression representing the specified subvector of the given band.
*/
template< AlignmentFlag AF  // Alignment flag
        , size_t... CSAs    // Compile time subvector arguments
        , typename VT       // Type of the vector
        , typename... RSAs  // Optional subvector arguments
        , EnableIf_t< IsBand_v< RemoveReference_t<VT> > &&
                      ( sizeof...( CSAs ) == 0UL || !RemoveReference_t<VT>::compileTimeArgs ) >* = nullptr >
inline decltype(auto) subvector( VT&& b, RSAs... args )
{
   BLAZE_FUNCTION_TRACE;

   const SubvectorData<CSAs...> sd( args... );

   const size_t row   ( b.row() + sd.offset() );
   const size_t column( b.column() + sd.offset() );

   constexpr auto check( getCheck( args... ) );

   try {
      return diagonal( submatrix<AF>( b.operand(), row, column, sd.size(), sd.size(), check ), check );
   }
   catch( ... ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BAND OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given band.
// \ingroup band
//
// \param band The band to be resetted.
// \return void
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void reset( Band<MT,TF,DF,MF,CBAs...>& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Resetting the given temporary band.
// \ingroup band
//
// \param band The temporary band to be resetted.
// \return void
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void reset( Band<MT,TF,DF,MF,CBAs...>&& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given band.
// \ingroup band
//
// \param band The band to be cleared.
// \return void
//
// Clearing a band is equivalent to resetting it via the reset() function.
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void clear( Band<MT,TF,DF,MF,CBAs...>& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Clearing the given temporary band.
// \ingroup band
//
// \param band The temporary band to be cleared.
// \return void
//
// Clearing a band is equivalent to resetting it via the reset() function.
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline void clear( Band<MT,TF,DF,MF,CBAs...>&& band )
{
   band.reset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given dense band is in default state.
// \ingroup band
//
// \param band The dense band to be tested for its default state.
// \return \a true in case the given dense band is component-wise zero, \a false otherwise.
//
// This function checks whether the dense band is in default state. For instance, in case the
// band is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all band elements are 0 and \a false in case any band element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( band( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF    // Relaxation flag
        , typename MT          // Type of the dense matrix
        , bool TF              // Transpose flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool isDefault( const Band<MT,TF,true,MF,CBAs...>& band )
{
   using blaze::isDefault;

   for( size_t i=0UL; i<band.size(); ++i )
      if( !isDefault<RF>( band[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given sparse band is in default state.
// \ingroup band
//
// \param band The sparse band to be tested for its default state.
// \return \a true in case the given band is component-wise zero, \a false otherwise.
//
// This function checks whether the sparse band is in default state. For instance, in case the
// band is instantiated for a built-in integral or floating point data type, the function returns
// \a true in case all band elements are 0 and \a false in case any band element is not 0. The
// following example demonstrates the use of the \a isDefault function:

   \code
   blaze::CompressedMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isDefault( band( A, 0UL ) ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF    // Relaxation flag
        , typename MT          // Type of the sparse matrix
        , bool TF              // Transpose flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool isDefault( const Band<MT,TF,false,MF,CBAs...>& band )
{
   using blaze::isDefault;

   for( const auto& element : band )
      if( !isDefault<RF>( element.value() ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the invariants of the given band are intact.
// \ingroup band
//
// \param band The band to be tested.
// \return \a true in case the given band's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the band are intact, i.e. if its state is valid.
// In case the invariants are intact, the function returns \a true, else it will return \a false.
// The following example demonstrates the use of the \a isIntact() function:

   \code
   blaze::DynamicMatrix<int,rowMajor> A;
   // ... Resizing and initialization
   if( isIntact( band( A, 0UL ) ) ) { ... }
   \endcode
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool isIntact( const Band<MT,TF,DF,MF,CBAs...>& band ) noexcept
{
   const ptrdiff_t index( band.band() );

   return ( ( index >= 0L || size_t( -index ) < band.operand().rows()    ) &&
            ( index <= 0L || size_t(  index ) < band.operand().columns() ) &&
            isIntact( band.operand() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two regular bands.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of two regular
// bands. In case both bands represent the same observable state, the function returns \a true,
// otherwise it returns \a false.
*/
template< typename MT1          // Type of the matrix of the left-hand side band
        , bool TF               // Transpose flag
        , bool DF               // Density flag
        , bool MF               // Multiplication flag
        , ptrdiff_t... CBAs1    // Compile time band arguments of the left-hand side band
        , typename MT2          // Type of the matrix of the right-hand side band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the right-hand side band
inline auto isSame_backend( const Band<MT1,TF,DF,MF,CBAs1...>& a,
                            const Band<MT2,TF,DF,MF,CBAs2...>& b ) noexcept
   -> DisableIf_t< IsSubmatrix_v<MT1> || IsSubmatrix_v<MT2>, bool >
{
   return ( isSame( a.operand(), b.operand() ) && ( a.band() == b.band() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the left band being a band on a submatrix.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the left
// band being a band on a submatrix. In case both bands represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1          // Type of the submatrix of the left-hand side band
        , bool TF               // Transpose flag
        , bool DF               // Density flag
        , bool MF               // Multiplication flag
        , ptrdiff_t... CBAs1    // Compile time band arguments of the left-hand side band
        , typename MT2          // Type of the matrix of the right-hand side band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the right-hand side band
inline auto isSame_backend( const Band<MT1,TF,DF,MF,CBAs1...>& a,
                            const Band<MT2,TF,DF,MF,CBAs2...>& b ) noexcept
   -> EnableIf_t< IsSubmatrix_v<MT1> && !IsSubmatrix_v<MT2>, bool >
{
   return ( isSame( a.operand().operand(), b.operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() ) &&
            ( a.column() + a.operand().column() == b.column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for the right band being a band on a submatrix.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of the right
// band being a band on a submatrix. In case both bands represent the same observable state,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename MT1          // Type of the matrix of the left-hand side band
        , bool TF               // Transpose flag
        , bool DF               // Density flag
        , bool MF               // Multiplication flag
        , ptrdiff_t... CBAs1    // Compile time band arguments of the left-hand side band
        , typename MT2          // Type of the submatrix of the right-hand side band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the right-hand side band
inline auto isSame_backend( const Band<MT1,TF,DF,MF,CBAs1...>& a,
                            const Band<MT2,TF,DF,MF,CBAs2...>& b ) noexcept
   -> EnableIf_t< !IsSubmatrix_v<MT1> && IsSubmatrix_v<MT2>, bool >
{
   return ( isSame( a.operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() == b.row() + b.operand().row() ) &&
            ( a.column() == b.column() + b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the isSame() function for two bands on submatrices.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This backend implementation of the isSame() function handles the special case of both bands
// being bands on submatrices. In case both bands represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1          // Type of the submatrix of the left-hand side band
        , bool TF               // Transpose flag
        , bool DF               // Density flag
        , bool MF               // Multiplication flag
        , ptrdiff_t... CBAs1    // Compile time band arguments of the left-hand side band
        , typename MT2          // Type of the submatrix of the right-hand side band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the right-hand side band
inline auto isSame_backend( const Band<MT1,TF,DF,MF,CBAs1...>& a,
                            const Band<MT2,TF,DF,MF,CBAs2...>& b ) noexcept
   -> EnableIf_t< IsSubmatrix_v<MT1> && IsSubmatrix_v<MT2>, bool >
{
   return ( isSame( a.operand().operand(), b.operand().operand() ) &&
            ( a.size() == b.size() ) &&
            ( a.row() + a.operand().row() == b.row() + b.operand().row() ) &&
            ( a.column() + a.operand().column() == b.column() + b.operand().column() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the two given bands represent the same observable state.
// \ingroup band
//
// \param a The first band to be tested for its state.
// \param b The second band to be tested for its state.
// \return \a true in case the two bands share a state, \a false otherwise.
//
// This overload of the isSame function tests if the two given bands refer to exactly the same
// range of the same matrix. In case both bands represent the same observable state, the function
// returns \a true, otherwise it returns \a false.
*/
template< typename MT1          // Type of the matrix of the left-hand side band
        , bool TF               // Transpose flag
        , bool DF               // Density flag
        , bool MF               // Multiplication flag
        , ptrdiff_t... CBAs1    // Compile time band arguments of the left-hand side band
        , typename MT2          // Type of the matrix of the right-hand side band
        , ptrdiff_t... CBAs2 >  // Compile time band arguments of the right-hand side band
inline bool isSame( const Band<MT1,TF,DF,MF,CBAs1...>& a,
                    const Band<MT2,TF,DF,MF,CBAs2...>& b ) noexcept
{
   return isSame_backend( a, b );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be set.
// \param value The value to be set to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool trySet( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return trySet( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !trySet( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The value to be added to the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryAdd( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryAdd( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryAdd( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The value to be subtracted from the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool trySub( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return trySub( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !trySub( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The factor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryMult( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryMult( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryMult( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryMult( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The divisor for the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryDiv( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryDiv( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by scaling a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryDiv( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryDiv( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param count The number of bits to shift the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline bool tryShift( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, int count )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryShift( band.operand(), band.row()+index, band.column()+index, count );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by shifting a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
BLAZE_ALWAYS_INLINE bool
   tryShift( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, int count )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryShift( band.operand(), band.row()+i, band.column()+i, count ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryBitand( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryBitand( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise AND on a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitand( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryBitand( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryBitor( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryBitor( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryBitor( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a band.
// \ingroup band
//
// \param band The target band.
// \param index The index of the element to be modified.
// \param value The bit pattern to be used on the element.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
inline bool tryBitxor( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index < band.size(), "Invalid vector access index" );

   return tryBitxor( band.operand(), band.row()+index, band.column()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a band.
// \ingroup band
//
// \param band The target band.
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
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename ET >      // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const Band<MT,TF,DF,MF,CBAs...>& band, size_t index, size_t size, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( index <= (*band).size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + size <= (*band).size(), "Invalid range size" );

   const size_t iend( index + size );

   for( size_t i=index; i<iend; ++i ) {
      if( !tryBitxor( band.operand(), band.row()+i, band.column()+i, value ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be assigned.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                       const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be added.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryAddAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryAddAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be subtracted.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool trySubAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return trySubAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the multiplication assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector to be multiplied.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryMultAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                           const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryMultAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the division assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector divisor.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryDivAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                          const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryDivAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the shift assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector of bits to shift.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryShiftAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                            const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryShiftAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise AND assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector for the bitwise AND operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryBitandAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                             const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitandAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector for the bitwise OR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryBitorAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                            const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitorAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to a band.
// \ingroup band
//
// \param lhs The target left-hand side band.
// \param rhs The right-hand side vector for the bitwise XOR operation.
// \param index The index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT        // Type of the matrix
        , bool TF            // Transpose flag
        , bool DF            // Density flag
        , bool MF            // Multiplication flag
        , ptrdiff_t... CBAs  // Compile time band arguments
        , typename VT >      // Type of the right-hand side vector
inline bool tryBitxorAssign( const Band<MT,TF,DF,MF,CBAs...>& lhs,
                             const Vector<VT,TF>& rhs, size_t index )
{
   BLAZE_INTERNAL_ASSERT( index <= lhs.size(), "Invalid vector access index" );
   BLAZE_INTERNAL_ASSERT( index + (*rhs).size() <= lhs.size(), "Invalid vector size" );

   return tryBitxorAssign( lhs.operand(), *rhs, lhs.band(), lhs.row()+index, lhs.column()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given band.
// \ingroup band
//
// \param b The band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given band. It returns a
// band object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT    // Type of the matrix
        , bool TF        // Transpose flag
        , bool DF        // Density flag
        , bool MF        // Multiplication flag
        , ptrdiff_t I >  // Band index
inline decltype(auto) derestrict( Band<MT,TF,DF,MF,I>& b )
{
   return band<I>( derestrict( b.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary band.
// \ingroup band
//
// \param b The temporary band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary band. It
// returns a band object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT    // Type of the matrix
        , bool TF        // Transpose flag
        , bool DF        // Density flag
        , bool MF        // Multiplication flag
        , ptrdiff_t I >  // Band index
inline decltype(auto) derestrict( Band<MT,TF,DF,MF,I>&& b )
{
   return band<I>( derestrict( b.operand() ), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given band.
// \ingroup band
//
// \param b The band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given band. It returns a
// band object that does provide the same interface but does not have any restrictions on the
// data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) derestrict( Band<MT,TF,DF,MF>& b )
{
   return band( derestrict( b.operand() ), b.band(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Removal of all restrictions on the data access to the given temporary band.
// \ingroup band
//
// \param b The temporary band to be derestricted.
// \return Band without access restrictions.
//
// This function removes all restrictions on the data access to the given temporary band. It
// returns a band object that does provide the same interface but does not have any restrictions
// on the data access.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the matrix
        , bool TF      // Transpose flag
        , bool DF      // Density flag
        , bool MF >    // Multiplication flag
inline decltype(auto) derestrict( Band<MT,TF,DF,MF>&& b )
{
   return band( derestrict( b.operand() ), b.band(), unchecked );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given band.
// \ingroup band
//
// \param b The given band.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given band.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline decltype(auto) unview( Band<MT,TF,DF,MF,CBAs...>& b )
{
   return b.operand();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the underlying matrix of the given constant band.
// \ingroup band
//
// \param b The given constant band.
// \return Reference to the underlying matrix.
//
// This function returns a reference to the underlying matrix of the given constant band.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT          // Type of the matrix
        , bool TF              // Transpose flag
        , bool DF              // Density flag
        , bool MF              // Multiplication flag
        , ptrdiff_t... CBAs >  // Compile time band arguments
inline decltype(auto) unview( const Band<MT,TF,DF,MF,CBAs...>& b )
{
   return b.operand();
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
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t I >
struct Size< Band<MT,TF,DF,MF,I>, 0UL >
   : public If_t< ( Size_v<MT,0UL> >= 0L && Size_v<MT,1UL> >= 0L )
                , Min_t< Ptrdiff_t< Size_v<MT,0UL> - ( I >= 0L ? 0L : -I ) >
                       , Ptrdiff_t< Size_v<MT,1UL> - ( I >= 0L ? I : 0L ) > >
                , Ptrdiff_t<-1L> >
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
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t I >
struct MaxSize< Band<MT,TF,DF,MF,I>, 0UL >
   : public If_t< ( MaxSize_v<MT,0UL> >= 0L && MaxSize_v<MT,1UL> >= 0L )
                , Min_t< Ptrdiff_t< MaxSize_v<MT,0UL> - ( I >= 0L ? 0L : -I ) >
                       , Ptrdiff_t< MaxSize_v<MT,1UL> - ( I >= 0L ? I : 0L ) > >
                , Ptrdiff_t<-1L> >
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
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t... CBAs >
struct IsRestricted< Band<MT,TF,DF,MF,CBAs...> >
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
template< typename MT, bool TF, bool MF, ptrdiff_t... CBAs >
struct HasConstDataAccess< Band<MT,TF,true,MF,CBAs...> >
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
template< typename MT, bool TF, bool MF, ptrdiff_t... CBAs >
struct HasMutableDataAccess< Band<MT,TF,true,MF,CBAs...> >
   : public HasMutableDataAccess<MT>
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
template< typename MT, bool TF, bool DF, bool MF, ptrdiff_t... CBAs >
struct IsOpposedView< Band<MT,TF,DF,MF,CBAs...> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
