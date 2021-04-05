//=================================================================================================
/*!
//  \file blaze/math/adaptors/StrictlyLowerMatrix.h
//  \brief Header file for the implementation of a strictly lower triangular matrix adaptor
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

#ifndef _BLAZE_MATH_ADAPTORS_STRICTLYLOWERMATRIX_H_
#define _BLAZE_MATH_ADAPTORS_STRICTLYLOWERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/adaptors/strictlylowermatrix/BaseTemplate.h>
#include <blaze/math/adaptors/strictlylowermatrix/Dense.h>
#include <blaze/math/adaptors/strictlylowermatrix/Sparse.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/Forward.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/DeclDiagTrait.h>
#include <blaze/math/traits/DeclHermTrait.h>
#include <blaze/math/traits/DeclLowTrait.h>
#include <blaze/math/traits/DeclStrLowTrait.h>
#include <blaze/math/traits/DeclStrUppTrait.h>
#include <blaze/math/traits/DeclSymTrait.h>
#include <blaze/math/traits/DeclUniLowTrait.h>
#include <blaze/math/traits/DeclUniUppTrait.h>
#include <blaze/math/traits/DeclUppTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasConstDataAccess.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsAdaptor.h>
#include <blaze/math/typetraits/IsAligned.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/MaxSize.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/YieldsDiagonal.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/MaybeUnused.h>


namespace blaze {

//=================================================================================================
//
//  STRICTLYLOWERMATRIX OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name StrictlyLowerMatrix operators */
//@{
template< typename MT, bool SO, bool DF >
void reset( StrictlyLowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
void reset( StrictlyLowerMatrix<MT,SO,DF>& m, size_t i );

template< typename MT, bool SO, bool DF >
void clear( StrictlyLowerMatrix<MT,SO,DF>& m );

template< RelaxationFlag RF, typename MT, bool SO, bool DF >
bool isDefault( const StrictlyLowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
bool isIntact( const StrictlyLowerMatrix<MT,SO,DF>& m );

template< typename MT, bool SO, bool DF >
void swap( StrictlyLowerMatrix<MT,SO,DF>& a, StrictlyLowerMatrix<MT,SO,DF>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be resetted.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( StrictlyLowerMatrix<MT,SO,DF>& m )
{
   m.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the specified row/column of the given strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be resetted.
// \param i The index of the row/column to be resetted.
// \return void
//
// This function resets the values in the specified row/column of the given strictly lower matrix
// to their default value. In case the given matrix is a \a rowMajor matrix the function resets
// the values in row \a i, if it is a \a columnMajor matrix the function resets the values in
// column \a i. Note that the capacity of the row/column remains unchanged.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void reset( StrictlyLowerMatrix<MT,SO,DF>& m, size_t i )
{
   m.reset( i );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be cleared.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void clear( StrictlyLowerMatrix<MT,SO,DF>& m )
{
   m.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given resizable strictly lower matrix is in default state.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the resizable strictly lower triangular matrix is in default
// state.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the adapted matrix
        , bool SO            // Storage order of the adapted matrix
        , bool DF >          // Density flag
inline bool isDefault_backend( const StrictlyLowerMatrix<MT,SO,DF>& m, TrueType )
{
   return ( m.rows() == 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the given fixed-size strictly lower matrix is in default state.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the fixed-size strictly lower triangular matrix is in default
// state.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the adapted matrix
        , bool SO            // Storage order of the adapted matrix
        , bool DF >          // Density flag
inline bool isDefault_backend( const StrictlyLowerMatrix<MT,SO,DF>& m, FalseType )
{
   if( SO ) {
      for( size_t j=0UL; j<m.columns(); ++j ) {
         for( size_t i=j+1UL; i<m.rows(); ++i ) {
            if( !isDefault<RF>( m(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t i=1UL; i<m.rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !isDefault<RF>( m(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given strictly lower matrix is in default state.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be tested for its default state.
// \return \a true in case the given matrix is in default state, \a false otherwise.
//
// This function checks whether the strictly lower triangular matrix is in default state. The
// following example demonstrates the use of the \a isDefault function:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyLowerMatrix;

   StrictlyLowerMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isDefault( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( A ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the adapted matrix
        , bool SO            // Storage order of the adapted matrix
        , bool DF >          // Density flag
inline bool isDefault( const StrictlyLowerMatrix<MT,SO,DF>& m )
{
   return isDefault_backend<RF>( m, typename IsResizable<MT>::Type() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given strictly lower matrix are intact.
// \ingroup strictly_lower_matrix
//
// \param m The strictly lower matrix to be tested.
// \return \a true in case the given matrix's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the strictly lower matrix are intact, i.e.
// if its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   using blaze::DynamicMatrix;
   using blaze::StrictlyLowerMatrix;

   StrictlyLowerMatrix< DynamicMatrix<int> > A;
   // ... Resizing and initialization
   if( isIntact( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline bool isIntact( const StrictlyLowerMatrix<MT,SO,DF>& m )
{
   return m.isIntact();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two matrices.
// \ingroup strictly_lower_matrix
//
// \param a The first matrix to be swapped.
// \param b The second matrix to be swapped.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline void swap( StrictlyLowerMatrix<MT,SO,DF>& a, StrictlyLowerMatrix<MT,SO,DF>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a single element of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool trySet( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( i < (*mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < (*mat).columns(), "Invalid column access index" );

   MAYBE_UNUSED( mat );

   return ( i > j || isDefault( value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by setting a range of elements of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
// \param row The index of the first row of the range to be multiplied.
// \param column The index of the first column of the range to be multiplied.
// \param m The number of rows of the range to be multiplied.
// \param n The number of columns of the range to be multiplied.
// \param value The value to be set to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySet( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   BLAZE_INTERNAL_ASSERT( row <= (*mat).rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= (*mat).columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + m <= (*mat).rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + n <= (*mat).columns(), "Invalid number of columns" );

   MAYBE_UNUSED( mat );

   return ( m == 0UL ) ||
          ( n == 0UL ) ||
          ( row >= column + n ) ||
          isDefault( value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a single element of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryAdd( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return trySet( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by adding to a range of elements of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
// \param row The index of the first row of the range to be multiplied.
// \param column The index of the first column of the range to be multiplied.
// \param m The number of rows of the range to be multiplied.
// \param n The number of columns of the range to be multiplied.
// \param value The value to be added to the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryAdd( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   return trySet( mat, row, column, m, n, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a single element of a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool trySub( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return trySet( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by subtracting from a range of elements of a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
// \param row The index of the first row of the range to be multiplied.
// \param column The index of the first column of the range to be multiplied.
// \param m The number of rows of the range to be multiplied.
// \param n The number of columns of the range to be multiplied.
// \param value The value to be subtracted from the range of elements.
// \return \a true in case the operation would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   trySub( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   return trySet( mat, row, column, m, n, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a single element of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryBitor( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return trySet( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise OR on a range of elements of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitor( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   return trySet( mat, row, column, m, n, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a single element of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
inline bool tryBitxor( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t i, size_t j, const ET& value )
{
   return tryAdd( mat, i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by a bitwise XOR on a range of elements of a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param mat The target strictly lower matrix.
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
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename ET >  // Type of the element
BLAZE_ALWAYS_INLINE bool
   tryBitxor( const StrictlyLowerMatrix<MT,SO,DF>& mat, size_t row, size_t column, size_t m, size_t n, const ET& value )
{
   return tryAdd( mat, row, column, m, n, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side dense vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side dense vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   MAYBE_UNUSED( lhs );

   if( column < row )
      return true;

   const size_t iend( min( column - row + 1UL, (*rhs).size() ) );

   for( size_t i=0UL; i<iend; ++i ) {
      if( !isDefault( (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side dense vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side dense vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                       const DenseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t ibegin( ( row <= column )?( 0UL ):( row - column ) );

   for( size_t i=ibegin; i<(*rhs).size(); ++i ) {
      if( !isDefault( (*rhs)[i] ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense vector to the band of
//        a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side dense vector to be assigned.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag of the right-hand side dense vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const DenseVector<VT,TF>& rhs,
                       ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs, row, column );

   if( band >= 0L ) {
      for( size_t i=0UL; i<(*rhs).size(); ++i ) {
         if( !isDefault( (*rhs)[i] ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side sparse vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side sparse vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );

   MAYBE_UNUSED( lhs );

   if( column < row )
      return true;

   const auto last( (*rhs).lowerBound( column - row + 1UL ) );

   for( auto element=(*rhs).begin(); element!=last; ++element ) {
      if( !isDefault( element->value() ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side sparse vector to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT    // Type of the adapted matrix
        , bool SO        // Storage order of the adapted matrix
        , bool DF        // Density flag
        , typename VT >  // Type of the right-hand side sparse vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                       const SparseVector<VT,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const auto last( (*rhs).end() );
   auto element( (*rhs).lowerBound( ( row <= column )?( 0UL ):( row - column ) ) );

   for( ; element!=last; ++element ) {
      if( !isDefault( element->value() ) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse vector to the band of
//        a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side sparse vector to be assigned.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side sparse vector
        , bool TF >    // Transpose flag of the right-hand side sparse vector
inline bool tryAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const SparseVector<VT,TF>& rhs,
                       ptrdiff_t band, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).size() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).size() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs, row, column );

   if( band >= 0L ) {
      for( const auto& element : *rhs ) {
         if( !isDefault( element.value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side dense matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool tryAssign( const StrictlyLowerMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t M( (*rhs).rows()    );
   const size_t N( (*rhs).columns() );

   if( row >= column + N )
      return true;

   const size_t iend( min( column + N - row, M ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row + i >= column );
      const size_t jbegin( ( containsDiagonal )?( row + i - column ):( 0UL ) );

      for( size_t j=jbegin; j<N; ++j ) {
         if( !isDefault( (*rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a dense matrix to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side dense matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side dense matrix
inline bool tryAssign( const StrictlyLowerMatrix<MT1,SO,DF>& lhs,
                       const DenseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t M( (*rhs).rows()    );
   const size_t N( (*rhs).columns() );

   if( row >= column + N )
      return true;

   const size_t jbegin( ( row <= column )?( 0UL ):( row - column ) );

   for( size_t j=jbegin; j<N; ++j )
   {
      const size_t iend( min( column + j - row + 1UL, M ) );

      for( size_t i=0UL; i<iend; ++i ) {
         if( !isDefault( (*rhs)(i,j) ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side sparse matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline bool tryAssign( const StrictlyLowerMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,false>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t M( (*rhs).rows()    );
   const size_t N( (*rhs).columns() );

   if( row >= column + N )
      return true;

   const size_t iend( min( column + N - row, M ) );

   for( size_t i=0UL; i<iend; ++i )
   {
      const bool containsDiagonal( row + i >= column );
      const size_t index( ( containsDiagonal )?( row + i - column ):( 0UL ) );

      const auto last( (*rhs).end(i) );
      auto element( (*rhs).lowerBound( i, index ) );

      for( ; element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the assignment of a sparse matrix to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side sparse matrix to be assigned.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1    // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename MT2 >  // Type of the right-hand side sparse matrix
inline bool tryAssign( const StrictlyLowerMatrix<MT1,SO,DF>& lhs,
                       const SparseMatrix<MT2,true>& rhs, size_t row, size_t column )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT2 );

   BLAZE_INTERNAL_ASSERT( row <= lhs.rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( column <= lhs.columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( row + (*rhs).rows() <= lhs.rows(), "Invalid number of rows" );
   BLAZE_INTERNAL_ASSERT( column + (*rhs).columns() <= lhs.columns(), "Invalid number of columns" );

   MAYBE_UNUSED( lhs );

   const size_t M( (*rhs).rows()    );
   const size_t N( (*rhs).columns() );

   if( row >= column + N )
      return true;

   const size_t jbegin( ( row < column )?( 0UL ):( row - column ) );

   for( size_t j=jbegin; j<N; ++j )
   {
      const size_t index( column + j - row + 1UL );
      const auto last( (*rhs).lowerBound( min( index, M ), j ) );

      for( auto element=(*rhs).begin(j); element!=last; ++element ) {
         if( !isDefault( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side vector to be added.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a vector to the band of
//        a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryAddAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the addition assignment of a matrix to a strictly lower
//        matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryAddAssign( const StrictlyLowerMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side vector to be subtracted.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool trySubAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                          const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a vector to the band of
//        a strictly lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool trySubAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                          ptrdiff_t band, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the subtraction assignment of a matrix to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool trySubAssign( const StrictlyLowerMatrix<MT1,SO1,DF>& lhs,
                          const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side vector for the bitwise OR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryBitorAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                            const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a vector to the band
//        of a strictly lower matrix.
// \ingroup strictly_lower_matrix
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryBitorAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                            ptrdiff_t band, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise OR assignment of a matrix to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryBitorAssign( const StrictlyLowerMatrix<MT1,SO1,DF>& lhs,
                            const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
// \param rhs The right-hand side vector for the bitwise XOR operation.
// \param row The row index of the first element to be modified.
// \param column The column index of the first element to be modified.
// \return \a true in case the assignment would be successful, \a false if not.
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryBitxorAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs,
                             const Vector<VT,TF>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a vector to the band
//        of a strictly lower matrix.
// \ingroup strictly_lower_matrix
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF      // Density flag
        , typename VT  // Type of the right-hand side vector
        , bool TF >    // Transpose flag of the right-hand side vector
inline bool tryBitxorAssign( const StrictlyLowerMatrix<MT,SO,DF>& lhs, const Vector<VT,TF>& rhs,
                            ptrdiff_t band, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, band, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Predict invariant violations by the bitwise XOR assignment of a matrix to a strictly
//        lower matrix.
// \ingroup strictly_lower_matrix
//
// \param lhs The target left-hand side strictly lower matrix.
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
template< typename MT1  // Type of the adapted matrix
        , bool SO1      // Storage order of the adapted matrix
        , bool DF       // Density flag
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline bool tryBitxorAssign( const StrictlyLowerMatrix<MT1,SO1,DF>& lhs,
                             const Matrix<MT2,SO2>& rhs, size_t row, size_t column )
{
   return tryAssign( lhs, *rhs, row, column );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a reference to the instance without the access restrictions to the upper part.
// \ingroup math_shims
//
// \param m The strictly lower matrix to be derestricted.
// \return Reference to the matrix without access restrictions.
//
// This function returns a reference to the given strictly lower matrix instance that has no
// access restrictions to the upper part of the matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in the violation of invariants, erroneous results and/or in compilation errors.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
inline MT& derestrict( StrictlyLowerMatrix<MT,SO,DF>& m )
{
   return m.matrix_;
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
template< typename MT, bool SO, bool DF >
struct Size< StrictlyLowerMatrix<MT,SO,DF>, 0UL >
   : public Size<MT,0UL>
{};

template< typename MT, bool SO, bool DF >
struct Size< StrictlyLowerMatrix<MT,SO,DF>, 1UL >
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
template< typename MT, bool SO, bool DF >
struct MaxSize< StrictlyLowerMatrix<MT,SO,DF>, 0UL >
   : public MaxSize<MT,0UL>
{};

template< typename MT, bool SO, bool DF >
struct MaxSize< StrictlyLowerMatrix<MT,SO,DF>, 1UL >
   : public MaxSize<MT,1UL>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSQUARE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsSquare< StrictlyLowerMatrix<MT,SO,DF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsUniform< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsUniform<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsSymmetric< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsUniform<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsHermitian< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsUniform<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsStrictlyLower< StrictlyLowerMatrix<MT,SO,DF> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsStrictlyUpper< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsUniform<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsAdaptor< StrictlyLowerMatrix<MT,SO,DF> >
   : public TrueType
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
template< typename MT, bool SO, bool DF >
struct IsRestricted< StrictlyLowerMatrix<MT,SO,DF> >
   : public TrueType
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
template< typename MT, bool SO >
struct HasConstDataAccess< StrictlyLowerMatrix<MT,SO,true> >
   : public TrueType
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
template< typename MT, bool SO, bool DF >
struct IsAligned< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsAligned<MT>
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
template< typename MT, bool SO, bool DF >
struct IsContiguous< StrictlyLowerMatrix<MT,SO,DF> >
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
template< typename MT, bool SO, bool DF >
struct IsPadded< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsPadded<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsResizable< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsResizable<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct IsShrinkable< StrictlyLowerMatrix<MT,SO,DF> >
   : public IsShrinkable<MT>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REMOVEADAPTOR SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct RemoveAdaptor< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = MT;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct AddTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( IsStrictlyLower_v<T1> && IsStrictlyLower_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename AddTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SubTraitEval1< T1, T2
                    , EnableIf_t< IsMatrix_v<T1> &&
                                  IsMatrix_v<T2> &&
                                  ( ( IsUniLower_v<T1> && IsUniLower_v<T2> ) ||
                                    ( IsStrictlyLower_v<T1> && IsStrictlyLower_v<T2> ) ) &&
                                  !( IsIdentity_v<T1> && IsIdentity_v<T2> ) &&
                                  !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename SubTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHURTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SchurTraitEval1< T1, T2
                      , EnableIf_t< IsMatrix_v<T1> &&
                                    IsMatrix_v<T2> &&
                                    ( ( IsStrictlyLower_v<T1> && !IsUpper_v<T2> ) ||
                                      ( !IsUpper_v<T1> && IsStrictlyLower_v<T2> ) ) &&
                                    !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename SchurTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsScalar_v<T2> &&
                                   ( IsStrictlyLower_v<T1> && !IsUniform_v<T1> ) > >
{
   using Type = StrictlyLowerMatrix< typename MultTraitEval2<T1,T2>::Type >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsScalar_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( IsStrictlyLower_v<T2> && !IsUniform_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename MultTraitEval2<T1,T2>::Type >;
};

template< typename T1, typename T2 >
struct MultTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( ( IsStrictlyLower_v<T1> && IsLower_v<T2> ) ||
                                     ( IsLower_v<T1> && IsStrictlyLower_v<T2> ) ) &&
                                   !( IsIdentity_v<T1> || IsIdentity_v<T2> ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename MultTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  KRONTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct KronTraitEval1< T1, T2
                     , EnableIf_t< IsMatrix_v<T1> &&
                                   IsMatrix_v<T2> &&
                                   ( ( IsStrictlyLower_v<T1> ) ||
                                     ( IsLower_v<T1> && IsStrictlyLower_v<T2> ) ) &&
                                   !( IsZero_v<T1> || IsZero_v<T2> ) > >
{
   using Type = StrictlyLowerMatrix< typename KronTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct DivTraitEval1< T1, T2
                    , EnableIf_t< IsStrictlyLower_v<T1> && IsScalar_v<T2> > >
{
   using Type = StrictlyLowerMatrix< typename DivTraitEval2<T1,T2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
struct UnaryMapTraitEval1< T, OP
                         , EnableIf_t< YieldsStrictlyLower_v<OP,T> &&
                                       !YieldsZero_v<OP,T> > >
{
   using Type = StrictlyLowerMatrix< typename UnaryMapTraitEval2<T,OP>::Type, StorageOrder_v<T> >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2, typename OP >
struct BinaryMapTraitEval1< T1, T2, OP
                          , EnableIf_t< YieldsStrictlyLower_v<OP,T1,T2> &&
                                        !YieldsDiagonal_v<OP,T1,T2> > >
{
   using Type = StrictlyLowerMatrix< typename BinaryMapTraitEval2<T1,T2,OP>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSYMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclSymTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = ZeroMatrix< typename MT::ElementType, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLHERMTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclHermTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = ZeroMatrix< typename MT::ElementType, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclLowTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = StrictlyLowerMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNILOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclUniLowTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = INVALID_TYPE;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRLOWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclStrLowTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = StrictlyLowerMatrix<MT,SO,DF>;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclUppTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = ZeroMatrix< typename MT::ElementType, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLUNIUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclUniUppTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = INVALID_TYPE;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLSTRUPPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclStrUppTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = ZeroMatrix< typename MT::ElementType, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DECLDIAGTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, bool SO, bool DF >
struct DeclDiagTrait< StrictlyLowerMatrix<MT,SO,DF> >
{
   using Type = ZeroMatrix< typename MT::ElementType, SO >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct HighType< StrictlyLowerMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   using Type = StrictlyLowerMatrix< typename HighType<MT1,MT2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT1, bool SO1, bool DF1, typename MT2, bool SO2, bool DF2 >
struct LowType< StrictlyLowerMatrix<MT1,SO1,DF1>, StrictlyLowerMatrix<MT2,SO2,DF2> >
{
   using Type = StrictlyLowerMatrix< typename LowType<MT1,MT2>::Type >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBMATRIXTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I, size_t N >
struct SubmatrixTraitEval1< MT, I, I, N, N
                          , EnableIf_t< I != inf && N != inf &&
                                        IsStrictlyLower_v<MT> &&
                                        !IsZero_v<MT> > >
{
   using Type = StrictlyLowerMatrix< typename SubmatrixTraitEval2<MT,I,I,N,N>::Type >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
