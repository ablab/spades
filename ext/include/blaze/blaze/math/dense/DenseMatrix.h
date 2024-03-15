//=================================================================================================
/*!
//  \file blaze/math/dense/DenseMatrix.h
//  \brief Header file for utility functions for dense matrices
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

#ifndef _BLAZE_MATH_DENSE_DENSEMATRIX_H_
#define _BLAZE_MATH_DENSE_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Triangular.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/Epsilon.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/clapack/potrf.h>
#include <blaze/math/ReductionFlag.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/shims/IsFinite.h>
#include <blaze/math/shims/IsInf.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/IsReal.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsRowMajorMatrix.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/RemoveAdaptor.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/views/Check.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Limits.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/NumericCast.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseMatrix operators */
//@{
template< typename T1, typename T2 >
auto operator==( const DenseMatrix<T1,false>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2 >
auto operator==( const DenseMatrix<T1,true>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2, bool SO >
auto operator==( T1 scalar, const DenseMatrix<T2,SO>& mat )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2, bool SO >
auto operator!=( const DenseMatrix<T1,SO>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2, bool SO >
auto operator!=( T1 scalar, const DenseMatrix<T2,SO>& mat )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename MT, bool SO, typename ST >
auto operator+=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator+=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator-=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator-=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator*=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator*=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator/=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator/=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO >
MT& operator<<=( DenseMatrix<MT,SO>& mat, int count );

template< typename MT, bool SO >
MT& operator<<=( DenseMatrix<MT,SO>&& mat, int count );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator<<=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator<<=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT, bool SO >
MT& operator>>=( DenseMatrix<MT,SO>& mat, int count );

template< typename MT, bool SO >
MT& operator>>=( DenseMatrix<MT,SO>&& mat, int count );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator>>=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator>>=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT, bool SO, typename ST >
auto operator&=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator&=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator&=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator&=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT, bool SO, typename ST >
auto operator|=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator|=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator|=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator|=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT, bool SO, typename ST >
auto operator^=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator^=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator^=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs );

template< typename MT1, bool SO1, typename MT2, bool SO2 >
MT1& operator^=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a row-major dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side row-major dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline auto operator==( const DenseMatrix<T1,false>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >
{
   using CT1 = CompositeType_t<T1>;

   // Evaluation of the dense matrix operand
   CT1 A( *mat );

   // In order to compare the matrix and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<A.rows(); ++i ) {
      for( size_t j=0; j<A.columns(); ++j ) {
         if( !equal( A(i,j), scalar ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a column-major dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side column-major dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1    // Type of the left-hand side dense matrix
        , typename T2 >  // Type of the right-hand side scalar
inline auto operator==( const DenseMatrix<T1,true>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >
{
   using CT1 = CompositeType_t<T1>;

   // Evaluation of the dense matrix operand
   CT1 A( *mat );

   // In order to compare the matrix and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t j=0; j<A.columns(); ++j ) {
      for( size_t i=0; i<A.rows(); ++i ) {
         if( !equal( A(i,j), scalar ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a scalar value and a dense matrix.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the comparison.
// \param mat The right-hand side dense matrix for the comparison.
// \return \a true if all elements of the matrix are equal to the scalar, \a false if not.
//
// If all values of the matrix are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline auto operator==( T1 scalar, const DenseMatrix<T2,SO>& mat )
   -> EnableIf_t< IsScalar_v<T1>, bool >
{
   return ( mat == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
//
// If one value of the matrix is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense matrix
        , typename T2  // Type of the right-hand side scalar
        , bool SO >    // Storage order
inline auto operator!=( const DenseMatrix<T1,SO>& mat, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >
{
   return !( mat == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a scalar value and a dense matrix.
// \ingroup dense_matrix
//
// \param scalar The left-hand side scalar value for the comparison.
// \param mat The right-hand side dense matrix for the comparison.
// \return \a true if at least one element of the matrix is different from the scalar, \a false if not.
//
// If one value of the matrix is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense matrix
        , bool SO >    // Storage order
inline auto operator!=( T1 scalar, const DenseMatrix<T2,SO>& mat )
   -> EnableIf_t< IsScalar_v<T1>, bool >
{
   return !( mat == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a dense matrix and a scalar value
//        (\f$ A+=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid addition to restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator+=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryAdd( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid addition to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left + scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a temporary dense matrix and a scalar
//        value (\f$ A+=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid addition to restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator+=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator+=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a dense matrix and a scalar value
//        (\f$ A-=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid subtraction from restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator-=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !trySub( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subtraction from restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left - scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a temporary dense matrix and a
//        scalar value (\f$ A-=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid subtraction from restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator-=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator-=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a dense matrix and
//        a scalar value (\f$ A*=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryMult( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left * scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a temporary dense matrix
//        and a scalar value (\f$ A*=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator*=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense matrix by a scalar value
//        (\f$ A/=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( isDivisor( scalar ), "Division by zero detected" );

   if( IsRestricted_v<MT> ) {
      if( !tryDiv( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left / scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a temporary dense matrix by a scalar
//        value (\f$ A/=s \f$).
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator/=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the uniform left-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param mat The dense matrix for the uniform left-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid left-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline MT& operator<<=( DenseMatrix<MT,SO>& mat, int count )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryShift( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), count ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid left-shift of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left << count );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the uniform left-shift of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param mat The temporary dense matrix for the uniform left-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid left-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline MT& operator<<=( DenseMatrix<MT,SO>&& mat, int count )
{
   return operator<<=( *mat, count );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the elementwise left-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid left-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator<<=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   if( IsRestricted_v<MT1> ) {
      if( !tryShiftAssign( *lhs, *rhs, 0UL, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid left-shift of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left << (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the elementwise left-shift of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side temporary dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid left-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator<<=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return operator<<=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the uniform right-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param mat The dense matrix for the uniform right-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid right-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline MT& operator>>=( DenseMatrix<MT,SO>& mat, int count )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryShift( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), count ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-shift of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left >> count );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the uniform right-shift of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param mat The temporary dense matrix for the uniform right-shift operation.
// \param count The number of bits to shift all matrix elements.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid right-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
inline MT& operator>>=( DenseMatrix<MT,SO>&& mat, int count )
{
   return operator>>=( *mat, count );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the elementwise right-shift of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid right-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator>>=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   if( IsRestricted_v<MT1> ) {
      if( !tryShiftAssign( *lhs, *rhs, 0UL, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-shift of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left >> (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the elementwise right-shift of a temporary dense.
// \ingroup dense_matrix
//
// \param lhs The left-hand side temporary dense matrix to be shifted.
// \param rhs The right-hand side dense matrix of bits to shift.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid right-shift of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator>>=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return operator>>=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise AND of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator&=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryBitand( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise AND of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left & scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a temporary dense matrix and a
//        scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise AND of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator&=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator&=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise AND operation.
// \param rhs The right-hand side dense matrix for the bitwise AND operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise AND of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator&=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   if( IsRestricted_v<MT1> ) {
      if( !tryBitandAssign( *lhs, *rhs, 0UL, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise AND of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left & (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side temporary dense matrix for the bitwise AND operation.
// \param rhs The right-hand side dense matrix for the bitwise AND operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise AND of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator&=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return operator&=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise OR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator|=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryBitor( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise OR of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left | scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a temporary dense matrix and a
//        scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise OR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator|=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator|=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise OR operation.
// \param rhs The right-hand side dense matrix for the bitwise OR operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise OR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator|=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   if( IsRestricted_v<MT1> ) {
      if( !tryBitorAssign( *lhs, *rhs, 0UL, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise OR of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left | (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side temporary dense matrix for the bitwise OR operation.
// \param rhs The right-hand side dense matrix for the bitwise OR operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise OR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator|=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return operator|=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a dense matrix and a scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side dense matrix for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise XOR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator^=( DenseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryBitxor( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise XOR of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *mat ) );

   smpAssign( left, left ^ scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a temporary dense matrix and a
//        scalar value.
// \ingroup dense_matrix
//
// \param mat The left-hand side temporary dense matrix for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return Reference to the left-hand side dense matrix.
// \exception std::invalid_argument Invalid bitwise XOR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side dense matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator^=( DenseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator^=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the bitwise XOR operation.
// \param rhs The right-hand side dense matrix for the bitwise XOR operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise XOR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator^=( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   if( IsRestricted_v<MT1> ) {
      if( !tryBitxorAssign( *lhs, *rhs, 0UL, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise XOR of restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left ^ (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a temporary dense matrix.
// \ingroup dense_matrix
//
// \param lhs The left-hand side temporary dense matrix for the bitwise XOR operation.
// \param rhs The right-hand side dense matrix for the bitwise XOR operation.
// \return Reference to the dense matrix.
// \exception std::invalid_argument Invalid bitwise XOR of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline MT1& operator^=( DenseMatrix<MT1,SO1>&& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return operator^=( *lhs, *rhs );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseMatrix functions */
//@{
template< typename MT, bool SO >
bool isnan( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
bool isinf( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
bool isfinite( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isSymmetric( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isHermitian( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniform( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isZero( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isLower( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniLower( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isStrictlyLower( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUpper( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniUpper( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isStrictlyUpper( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isDiagonal( const DenseMatrix<MT,SO>& dm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isIdentity( const DenseMatrix<MT,SO>& dm );

template< typename MT, bool SO >
bool isPositiveDefinite( const DenseMatrix<MT,SO>& dm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense matrix for not-a-number elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the matrix is not-a-number, \a false otherwise.
//
// This function checks the dense matrix for not-a-number (NaN) elements. If at least one
// element of the matrix is not-a-number, the function returns \a true, otherwise it returns
// \a false.

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isnan( const DenseMatrix<MT,SO>& dm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return false;

   using CT = CompositeType_t<MT>;

   CT A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         const size_t jbegin( IsUpper_v<MT> ? i : 0UL );
         const size_t jend  ( IsLower_v<MT> ? i+1UL : A.columns() );
         for( size_t j=jbegin; j<jend; ++j )
            if( isnan( A(i,j) ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         const size_t ibegin( IsLower_v<MT> ? j : 0UL );
         const size_t iend  ( IsUpper_v<MT> ? j+1UL : A.rows() );
         for( size_t i=ibegin; i<iend; ++i )
            if( isnan( A(i,j) ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense matrix for infinite elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked for infinite elements.
// \return \a true if at least one element of the matrix is infinite, \a false otherwise.
//
// This function checks the dense matrix for infinite (NaN) elements. If at least one element
// of the matrix is infinite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isinf( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isinf( const DenseMatrix<MT,SO>& dm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return false;

   using CT = CompositeType_t<MT>;

   CT A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         const size_t jbegin( IsUpper_v<MT> ? i : 0UL );
         const size_t jend  ( IsLower_v<MT> ? i+1UL : A.columns() );
         for( size_t j=jbegin; j<jend; ++j )
            if( isinf( A(i,j) ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         const size_t ibegin( IsLower_v<MT> ? j : 0UL );
         const size_t iend  ( IsUpper_v<MT> ? j+1UL : A.rows() );
         for( size_t i=ibegin; i<iend; ++i )
            if( isinf( A(i,j) ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense matrix for finite elements.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked for finite elements.
// \return \a true if all elements of the matrix are finite, \a false otherwise.
//
// This function checks if all elements of the dense matrix are finite elements (i.e. normal,
// subnormal or zero elements, but not infinite or NaN). If all elements of the matrix are
// finite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::DynamicMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isfinite( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isfinite( const DenseMatrix<MT,SO>& dm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return true;

   using CT = CompositeType_t<MT>;

   CT A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         const size_t jbegin( IsUpper_v<MT> ? i : 0UL );
         const size_t jend  ( IsLower_v<MT> ? i+1UL : A.columns() );
         for( size_t j=jbegin; j<jend; ++j )
            if( !isfinite( A(i,j) ) ) return false;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         const size_t ibegin( IsLower_v<MT> ? j : 0UL );
         const size_t iend  ( IsUpper_v<MT> ? j+1UL : A.rows() );
         for( size_t i=ibegin; i<iend; ++i )
            if( !isfinite( A(i,j) ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is symmetric.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is symmetric, \a false if not.
//
// This function checks if the given dense matrix is symmetric. The matrix is considered to be
// symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$). The
// following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isSymmetric( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isSymmetric<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a symmetric matrix:

   \code
   if( isSymmetric( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isSymmetric( const DenseMatrix<MT,SO>& dm )
{
   using CT = CompositeType_t<MT>;

   if( IsSymmetric_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsUniform_v<MT> || (*dm).rows() < 2UL )
      return true;

   if( IsTriangular_v<MT> )
      return isDiagonal( *dm );

   CT A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=1UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !equal<RF>( A(i,j), A(j,i) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=1UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !equal<RF>( A(i,j), A(j,i) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is Hermitian.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is Hermitian, \a false if not.
//
// This function checks if the given dense matrix is an Hermitian matrix. The matrix is considered
// to be an Hermitian matrix if it is a square matrix whose conjugate transpose is equal to itself
// (\f$ A = \overline{A^T} \f$), i.e. each matrix element \f$ a_{ij} \f$ is equal to the complex
// conjugate of the element \f$ a_{ji} \f$. The following code example demonstrates the use of the
// function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isHermitian( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isHermitian<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an Hermitian matrix:

   \code
   if( isHermitian( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isHermitian( const DenseMatrix<MT,SO>& dm )
{
   using ET = ElementType_t<MT>;
   using CT = CompositeType_t<MT>;

   if( IsHermitian_v<MT> )
      return true;

   if( !IsScalar_v<ET> || !isSquare( *dm ) )
      return false;

   if( IsBuiltin_v<ET> && IsUniform_v<MT> )
      return true;

   CT A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !equal<RF>( A(i,j), conj( A(j,i) ) ) )
               return false;
         }
         if( !isReal<RF>( A(i,i) ) )
            return false;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !equal<RF>( A(i,j), conj( A(j,i) ) ) )
               return false;
         }
         if( !isReal<RF>( A(j,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given row-major triangular dense matrix is a uniform matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the dense matrix
bool isUniform_backend( const DenseMatrix<MT,false>& dm, TrueType )
{
   BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t ibegin( ( IsStrictlyLower_v<MT> )?( 1UL ):( 0UL ) );
   const size_t iend  ( ( IsStrictlyUpper_v<MT> )?( (*dm).rows()-1UL ):( (*dm).rows() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      if( !IsUpper_v<MT> ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !isDefault<RF>( (*dm)(i,j) ) )
               return false;
         }
      }
      if( !isDefault<RF>( (*dm)(i,i) ) )
         return false;
      if( !IsLower_v<MT> ) {
         for( size_t j=i+1UL; j<(*dm).columns(); ++j ) {
            if( !isDefault<RF>( (*dm)(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given column-major triangular dense matrix is a uniform matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the dense matrix
bool isUniform_backend( const DenseMatrix<MT,true>& dm, TrueType )
{
   BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t jbegin( ( IsStrictlyUpper_v<MT> )?( 1UL ):( 0UL ) );
   const size_t jend  ( ( IsStrictlyLower_v<MT> )?( (*dm).columns()-1UL ):( (*dm).columns() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      if( !IsLower_v<MT> ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !isDefault<RF>( (*dm)(i,j) ) )
               return false;
         }
      }
      if( !isDefault<RF>( (*dm)(j,j) ) )
         return false;
      if( !IsUpper_v<MT> ) {
         for( size_t i=j+1UL; i<(*dm).rows(); ++i ) {
            if( !isDefault<RF>( (*dm)(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given row-major general dense matrix is a uniform matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the dense matrix
bool isUniform_backend( const DenseMatrix<MT,false>& dm, FalseType )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() != 0UL, "Invalid number of columns detected" );

   const auto& cmp( (*dm)(0UL,0UL) );

   for( size_t i=0UL; i<(*dm).rows(); ++i ) {
      for( size_t j=0UL; j<(*dm).columns(); ++j ) {
         if( !equal<RF>( (*dm)(i,j), cmp ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given column-major general dense matrix is a uniform matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the dense matrix
bool isUniform_backend( const DenseMatrix<MT,true>& dm, FalseType )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() != 0UL, "Invalid number of columns detected" );

   const auto& cmp( (*dm)(0UL,0UL) );

   for( size_t j=0UL; j<(*dm).columns(); ++j ) {
      for( size_t i=0UL; i<(*dm).rows(); ++i ) {
         if( !equal<RF>( (*dm)(i,j), cmp ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a uniform matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
//
// This function checks if the given dense matrix is a uniform matrix. The matrix is considered
// to be uniform if all its elements are identical. The following code example demonstrates the
// use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniform( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniform<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a uniform matrix:

   \code
   if( isUniform( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isUniform( const DenseMatrix<MT,SO>& dm )
{
   if( IsUniform_v<MT> ||
       (*dm).rows() == 0UL || (*dm).columns() == 0UL ||
       ( (*dm).rows() == 1UL && (*dm).columns() == 1UL ) )
      return true;

   if( IsUniTriangular_v<MT> )
      return false;

   CompositeType_t<MT> A( *dm );  // Evaluation of the dense matrix operand

   return isUniform_backend<RF>( A, typename IsTriangular<MT>::Type() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a zero matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a zero matrix, \a false if not.
//
// This function checks if the given dense matrix is a zero matrix. The matrix is considered to
// be zero if all its elements are zero. The following code example demonstrates the use of the
// function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isZero( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isZero<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a zero matrix:

   \code
   if( isZero( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isZero( const DenseMatrix<MT,SO>& dm )
{
   const size_t M( (*dm).rows()    );
   const size_t N( (*dm).columns() );

   if( IsZero_v<MT> || M == 0UL || N == 0UL )
      return true;

   if( IsUniTriangular_v<MT> )
      return false;

   if( IsUniform_v<MT> )
      return isZero<RF>( (*dm)(0UL,0UL) );

   CompositeType_t<MT> A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor )
   {
      for( size_t i=0UL; i<M; ++i )
      {
         const size_t jbegin( IsUpper_v<MT>
                              ? ( IsStrictlyUpper_v<MT> ? i+1UL : i )
                              : 0UL );
         const size_t jend  ( IsLower_v<MT> || IsSymmetric_v<MT> || IsHermitian_v<MT>
                              ? ( IsStrictlyLower_v<MT> ? i : i+1UL )
                              : N );

         for( size_t j=jbegin; j<jend; ++j ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<N; ++j )
      {
         const size_t ibegin( IsLower_v<MT>
                              ? ( IsStrictlyLower_v<MT> ? j+1UL : j )
                              : 0UL );
         const size_t iend  ( IsUpper_v<MT> || IsSymmetric_v<MT> || IsHermitian_v<MT>
                              ? ( IsStrictlyUpper_v<MT> ? j : j+1UL )
                              : M );

         for( size_t i=ibegin; i<iend; ++i ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a lower triangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given dense matrix is a lower triangular matrix. The matrix is
// considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially lower triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a lower triangular matrix:

   \code
   if( isLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isLower( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsLower_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsZero_v<MT> || (*dm).rows() < 2UL )
      return true;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( IsUniform_v<MT> )
      return isDefault<RF>( A(0UL,0UL) );

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows()-1UL; ++i ) {
         for( size_t j=i+1UL; j<A.columns(); ++j ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=1UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a lower unitriangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a lower unitriangular matrix, \a false if not.
//
// This function checks if the given dense matrix is a lower unitriangular matrix. The matrix is
// considered to be lower unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a lower unitriangular matrix:

   \code
   if( isUniLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isUniLower( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUniLower_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         if( !isOne<RF>( A(i,i) ) )
            return false;
         for( size_t j=i+1UL; j<A.columns(); ++j ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
         if( !isOne<RF>( A(j,j) ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a strictly lower triangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a strictly lower triangular matrix, \a false if not.
//
// This function checks if the given dense matrix is a strictly lower triangular matrix. The
// matrix is considered to be strictly lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a strictly lower triangular
// matrix:

   \code
   if( isStrictlyLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isStrictlyLower( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsStrictlyLower_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsZero_v<MT> || (*dm).rows() < 2UL )
      return true;

   if( IsUniLower_v<MT> || IsUniUpper_v<MT> )
      return false;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( IsUniform_v<MT> )
      return isDefault<RF>( A(0UL,0UL) );

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=i; j<A.columns(); ++j ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( size_t i=0UL; i<=j; ++i ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is an upper triangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is an upper triangular matrix, \a false if not.
//
// This function checks if the given dense matrix is an upper triangular matrix. The matrix is
// considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially upper triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper triangular matrix:

   \code
   if( isUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isUpper( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUpper_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsZero_v<MT> || (*dm).rows() < 2UL )
      return true;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( IsUniform_v<MT> )
      return isDefault<RF>( A(0UL,0UL) );

   if( SO == rowMajor ) {
      for( size_t i=1UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns()-1UL; ++j ) {
         for( size_t i=j+1UL; i<A.rows(); ++i ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is an upper unitriangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is an upper unitriangular matrix, \a false if not.
//
// This function checks if the given dense matrix is an upper unitriangular matrix. The matrix is
// considered to be upper unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper unitriangular matrix:

   \code
   if( isUniUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isUniUpper( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUniUpper_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
         if( !isOne<RF>( A(i,i) ) )
            return false;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         if( !isOne<RF>( A(j,j) ) )
            return false;
         for( size_t i=j+1UL; i<A.rows(); ++i ) {
            if( !isZero<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a strictly upper triangular matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a strictly upper triangular matrix, \a false if not.
//
// This function checks if the given dense matrix is a strictly upper triangular matrix. The
// matrix is considered to be strictly upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a strictly upper triangular
// matrix:

   \code
   if( isStrictlyUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isStrictlyUpper( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsStrictlyUpper_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsZero_v<MT> || (*dm).rows() < 2UL )
      return true;

   if( IsUniLower_v<MT> || IsUniUpper_v<MT> )
      return false;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( IsUniform_v<MT> )
      return isDefault<RF>( A(0UL,0UL) );

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( size_t j=0UL; j<=i; ++j ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( size_t i=j; i<A.rows(); ++i ) {
            if( !isDefault<RF>( A(i,j) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is diagonal.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is diagonal, \a false if not.
//
// This function tests whether the matrix is diagonal, i.e. if the non-diagonal elements are
// default elements. In case of integral or floating point data types, a diagonal matrix has
// the form

                        \f[\left(\begin{array}{*{5}{c}}
                        aa     & 0      & 0      & \cdots & 0  \\
                        0      & bb     & 0      & \cdots & 0  \\
                        0      & 0      & cc     & \cdots & 0  \\
                        \vdots & \vdots & \vdots & \ddots & 0  \\
                        0      & 0      & 0      & 0      & xx \\
                        \end{array}\right)\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially diagonal. The
// following example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isDiagonal( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDiagonal<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a diagonal matrix:

   \code
   if( isDiagonal( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isDiagonal( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsDiagonal_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( IsZero_v<MT> || (*dm).rows() < 2UL )
      return true;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( IsUniform_v<MT> )
      return isDefault<RF>( A(0UL,0UL) );

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         if( !IsUpper_v<MT> ) {
            for( size_t j=0UL; j<i; ++j ) {
               if( !isDefault<RF>( A(i,j) ) )
                  return false;
            }
         }
         if( !IsLower_v<MT> ) {
            for( size_t j=i+1UL; j<A.columns(); ++j ) {
               if( !isDefault<RF>( A(i,j) ) )
                  return false;
            }
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         if( !IsLower_v<MT> ) {
            for( size_t i=0UL; i<j; ++i ) {
               if( !isDefault<RF>( A(i,j) ) )
                  return false;
            }
         }
         if( !IsUpper_v<MT> ) {
            for( size_t i=j+1UL; i<A.rows(); ++i ) {
               if( !isDefault<RF>( A(i,j) ) )
                  return false;
            }
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is an identity matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is an identity matrix, \a false if not.
//
// This function tests whether the matrix is an identity matrix, i.e. if the diagonal elements
// are 1 and the non-diagonal elements are 0. In case of integral or floating point data types,
// an identity matrix has the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & 0      & 0      & \cdots & 0 \\
                        0      & 1      & 0      & \cdots & 0 \\
                        0      & 0      & 1      & \cdots & 0 \\
                        \vdots & \vdots & \vdots & \ddots & 0 \\
                        0      & 0      & 0      & 0      & 1 \\
                        \end{array}\right)\f]

// The following example demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isIdentity( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isIdentity<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an identity matrix:

   \code
   if( isIdentity( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the dense matrix
        , bool SO >          // Storage order
bool isIdentity( const DenseMatrix<MT,SO>& dm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsIdentity_v<MT> )
      return true;

   if( !isSquare( *dm ) )
      return false;

   if( (*dm).rows() == 0UL )
      return true;

   Tmp A( *dm );  // Evaluation of the dense matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         if( !IsUpper_v<MT> ) {
            for( size_t j=0UL; j<i; ++j ) {
               if( !isZero<RF>( A(i,j) ) )
                  return false;
            }
         }
         if( !IsUniLower_v<MT> && !IsUniUpper_v<MT> && !isOne<RF>( A(i,i) ) ) {
            return false;
         }
         if( !IsLower_v<MT> ) {
            for( size_t j=i+1UL; j<A.columns(); ++j ) {
               if( !isZero<RF>( A(i,j) ) )
                  return false;
            }
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         if( !IsLower_v<MT> ) {
            for( size_t i=0UL; i<j; ++i ) {
               if( !isZero<RF>( A(i,j) ) )
                  return false;
            }
         }
         if( !IsUniLower_v<MT> && !IsUniUpper_v<MT> && !isOne<RF>( A(j,j) ) ) {
            return false;
         }
         if( !IsUpper_v<MT> ) {
            for( size_t i=j+1UL; i<A.rows(); ++i ) {
               if( !isZero<RF>( A(i,j) ) )
                  return false;
            }
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense matrix is a positive definite matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be checked.
// \return \a true if the matrix is a positive definite matrix, \a false if not.
//
// This function tests whether the matrix is a positive definite matrix. The following example
// demonstrates the use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isPositiveDefinite( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a positive definite matrix:

   \code
   if( isPositiveDefinite( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
//
// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
bool isPositiveDefinite( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   if( !isSquare( *dm ) )
      return false;

   if( (*dm).rows() < 2UL )
      return true;

   RemoveAdaptor_t< ResultType_t<MT> > L( *dm );

   char uplo( IsRowMajorMatrix_v<MT> ? 'U' : 'L' );
   int  n   ( numeric_cast<int>( (*L).rows()    ) );
   int  lda ( numeric_cast<int>( (*L).spacing() ) );
   int  info( 0 );

   potrf( uplo, n, (*L).data(), lda, &info );

   return ( info == 0 );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the rank of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix for the rank computation.
// \return The rank of the given dense matrix.
//
// This function computes the rank of the given dense matrix \a dm. The rank is determined as
// the number of singular values greater than a given tolerance. This tolerance is computed as

   \code
   tolerance = max(m,n) * max(s) * epsilon,
   \endcode

// where \c m is the number of rows of \a dm, \c n is the number of columns of \a dm, \c max(s)
// is the maximum singular value of \a dm and epsilon is the difference between 1 and the least
// value greater than 1 that is representable by the floating point type of the singular values.
// Example:

   \code
   blaze::DynamicMatrix<double> A( 5UL, 8UL );
   // ... Initialization
   rank( A );
   \endcode

// \note This function only works for matrices with \c float, \c double, \c complex<float>, or
// \c complex<double> element type. The attempt to call the function with matrices of any other
// element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a call to this function will result in a linker error.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order
size_t rank( const DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   using BT = UnderlyingBuiltin_t<MT>;

   const auto s( evaluate( svd( *dm ) ) );

   const BT tolerance( max( rows(dm), columns(dm) ) * Limits<BT>::epsilon() * max(s) );
   return std::count_if( begin(s), end(s), [tolerance]( auto val ) { return abs(val) > tolerance; } );
}
//*************************************************************************************************

} // namespace blaze

#endif
