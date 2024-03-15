//=================================================================================================
/*!
//  \file blaze/math/dense/DenseVector.h
//  \brief Header file for utility functions for dense vectors
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

#ifndef _BLAZE_MATH_DENSE_DENSEVECTOR_H_
#define _BLAZE_MATH_DENSE_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/shims/IsFinite.h>
#include <blaze/math/shims/IsInf.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name DenseVector operators */
//@{
template< typename T1, typename T2, bool TF >
auto operator==( const DenseVector<T1,TF>& vec, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2, bool TF >
auto operator==( T1 scalar, const DenseVector<T2,TF>& vec )
   -> EnableIf_t< IsScalar_v<T1>, bool >;

template< typename T1, typename T2, bool TF >
auto operator!=( const DenseVector<T1,TF>& vec, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >;

template< typename T1, typename T2, bool TF >
auto operator!=( T1 scalar, const DenseVector<T2,TF>& vec )
   -> EnableIf_t< IsScalar_v<T1>, bool >;

template< typename VT, bool TF, typename ST >
auto operator+=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator+=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator-=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator-=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator*=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator*=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator/=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator/=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF >
VT& operator<<=( DenseVector<VT,TF>& vec, int count );

template< typename VT, bool TF >
VT& operator<<=( DenseVector<VT,TF>&& vec, int count );

template< typename VT1, typename VT2, bool TF >
VT1& operator<<=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT1, typename VT2, bool TF >
VT1& operator<<=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT, bool TF >
VT& operator>>=( DenseVector<VT,TF>& vec, int count );

template< typename VT, bool TF >
VT& operator>>=( DenseVector<VT,TF>&& vec, int count );

template< typename VT1, typename VT2, bool TF >
VT1& operator>>=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT1, typename VT2, bool TF >
VT1& operator>>=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT, bool TF, typename ST >
auto operator&=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator&=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT1, typename VT2, bool TF >
VT1& operator&=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT1, typename VT2, bool TF >
VT1& operator&=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT, bool TF, typename ST >
auto operator|=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator|=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT1, typename VT2, bool TF >
VT1& operator|=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT1, typename VT2, bool TF >
VT1& operator|=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT, bool TF, typename ST >
auto operator^=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator^=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT1, typename VT2, bool TF >
VT1& operator^=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs );

template< typename VT1, typename VT2, bool TF >
VT1& operator^=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if all elements of the vector are equal to the scalar, \a false if not.
//
// If all values of the vector are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline auto operator==( const DenseVector<T1,TF>& vec, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >
{
   using CT1 = CompositeType_t<T1>;

   // Evaluation of the dense vector operand
   CT1 a( *vec );

   // In order to compare the vector and the scalar value, the data values of the lower-order
   // data type are converted to the higher-order data type within the equal function.
   for( size_t i=0; i<a.size(); ++i )
      if( !equal( a[i], scalar ) ) return false;
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of a scalar value and a dense vector.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the comparison.
// \param vec The right-hand side dense vector for the comparison.
// \return \a true if all elements of the vector are equal to the scalar, \a false if not.
//
// If all values of the vector are equal to the scalar value, the equality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side dense vector
        , bool TF >    // Transpose flag
inline auto operator==( T1 scalar, const DenseVector<T2,TF>& vec )
   -> EnableIf_t< IsScalar_v<T1>, bool >
{
   return ( vec == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the comparison.
// \param scalar The right-hand side scalar value for the comparison.
// \return \a true if at least one element of the vector is different from the scalar, \a false if not.
//
// If one value of the vector is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side dense vector
        , typename T2  // Type of the right-hand side scalar
        , bool TF >    // Transpose flag
inline auto operator!=( const DenseVector<T1,TF>& vec, T2 scalar )
   -> EnableIf_t< IsScalar_v<T2>, bool >
{
   return !( vec == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of a scalar value and a dense vector.
// \ingroup dense_vector
//
// \param scalar The left-hand side scalar value for the comparison.
// \param vec The right-hand side dense vector for the comparison.
// \return \a true if at least one element of the vector is different from the scalar, \a false if not.
//
// If one value of the vector is inequal to the scalar value, the inequality test returns \a true,
// otherwise \a false. Note that this function can only be used with built-in, numerical data
// types!
*/
template< typename T1  // Type of the left-hand side scalar
        , typename T2  // Type of the right-hand side vector
        , bool TF >    // Transpose flag
inline auto operator!=( T1 scalar, const DenseVector<T2,TF>& vec )
   -> EnableIf_t< IsScalar_v<T1>, bool >
{
   return !( vec == scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a dense vector and a scalar value
//        (\f$ \vec{a}+=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid addition to restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator+=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryAdd( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid addition to restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left + scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a temporary dense vector and a scalar
//        value (\f$ \vec{v}+=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the addition.
// \param scalar The right-hand side scalar value for the addition.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid addition to restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator+=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator+=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a dense vector and a scalar value
//        (\f$ \vec{a}-=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid subtraction from restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator-=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !trySub( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subtraction from restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left - scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a temporary dense vector and a
//        scalar value (\f$ \vec{v}-=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the subtraction.
// \param scalar The right-hand side scalar value for the subtraction.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid subtraction from restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator-=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator-=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a dense vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryMult( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left * scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a temporary dense vector
//        and a scalar value (\f$ \vec{v}*=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator*=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   BLAZE_USER_ASSERT( isDivisor( scalar ), "Division by zero detected" );

   if( IsRestricted_v<VT> ) {
      if( !tryDiv( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left / scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a temporary dense vector by a scalar
//        value (\f$ \vec{a}/=s \f$).
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator/=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the uniform left-shift of a dense vector.
// \ingroup dense_vector
//
// \param vec The dense vector for the uniform left-shift operation.
// \param count The number of bits to shift all vector elements.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid left-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline VT& operator<<=( DenseVector<VT,TF>& vec, int count )
{
   if( IsRestricted_v<VT> ) {
      if( !tryShift( *vec, 0UL, (*vec).size(), count ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid left-shift of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left << count );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the uniform shift of a temporary dense vector.
// \ingroup dense_vector
//
// \param vec The temporary dense vector for the uniform left-shift operation.
// \param count The number of bits to shift all vector elements.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid left-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline VT& operator<<=( DenseVector<VT,TF>&& vec, int count )
{
   return operator<<=( *vec, count );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the elementwise left-shift of a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector to be shifted.
// \param rhs The right-hand side dense vector of bits to shift.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid left-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator<<=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   if( IsRestricted_v<VT1> ) {
      if( !tryShiftAssign( *lhs, *rhs, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid left-shift of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left << (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Left-shift assignment operator for the elementwise left-shift of a temporary dense
//        vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side temporary dense vector to be shifted.
// \param rhs The right-hand side dense vector of bits to shift.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid left-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator<<=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs )
{
   return operator<<=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the uniform right-shift of a dense vector.
// \ingroup dense_vector
//
// \param vec The dense vector for the uniform right-shift operation.
// \param count The number of bits to shift all vector elements.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid right-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline VT& operator>>=( DenseVector<VT,TF>& vec, int count )
{
   if( IsRestricted_v<VT> ) {
      if( !tryShift( *vec, 0UL, (*vec).size(), count ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-shift of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left >> count );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the uniform shift of a temporary dense vector.
// \ingroup dense_vector
//
// \param vec The temporary dense vector for the uniform right-shift operation.
// \param count The number of bits to shift all vector elements.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid right-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
inline VT& operator>>=( DenseVector<VT,TF>&& vec, int count )
{
   return operator>>=( *vec, count );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Right-shift assignment operator for the elementwise right-shift of a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector to be shifted.
// \param rhs The right-hand side dense vector of bits to shift.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid right-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator>>=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   if( IsRestricted_v<VT1> ) {
      if( !tryShiftAssign( *lhs, *rhs, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid right-shift of restricted vector" );
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
// \ingroup dense_vector
//
// \param lhs The left-hand side temporary dense vector to be shifted.
// \param rhs The right-hand side dense vector of bits to shift.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid right-shift of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator>>=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs )
{
   return operator>>=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise AND of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator&=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryBitand( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise AND of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left & scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a temporary dense vector and
//        a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the bitwise AND.
// \param scalar The right-hand side scalar value for the bitwise AND.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise AND of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator&=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator&=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the bitwise AND operation.
// \param rhs The right-hand side dense vector for the bitwise AND operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise AND of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator&=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   if( IsRestricted_v<VT1> ) {
      if( !tryBitandAssign( *lhs, *rhs, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise AND of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left & (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise AND assignment operator for the bitwise AND of a temporary dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side temporary dense vector for the bitwise AND operation.
// \param rhs The right-hand side dense vector for the bitwise AND operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise AND of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator&=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs )
{
   return operator&=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise OR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator|=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryBitor( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise OR of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left | scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a temporary dense vector and
//        a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the bitwise OR.
// \param scalar The right-hand side scalar value for the bitwise OR.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise OR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator|=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator|=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the bitwise OR operation.
// \param rhs The right-hand side dense vector for the bitwise OR operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise OR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator|=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   if( IsRestricted_v<VT1> ) {
      if( !tryBitorAssign( *lhs, *rhs, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise OR of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left | (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise OR assignment operator for the bitwise OR of a temporary dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side temporary dense vector for the bitwise OR operation.
// \param rhs The right-hand side dense vector for the bitwise OR operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise OR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator|=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs )
{
   return operator|=( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a dense vector and a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side dense vector for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise XOR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator^=( DenseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryBitxor( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise XOR of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *vec ) );

   smpAssign( left, left ^ scalar );

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a temporary dense vector and
//        a scalar value.
// \ingroup dense_vector
//
// \param vec The left-hand side temporary dense vector for the bitwise XOR.
// \param scalar The right-hand side scalar value for the bitwise XOR.
// \return Reference to the left-hand side dense vector.
// \exception std::invalid_argument Invalid bitwise XOR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side dense vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator^=( DenseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator^=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the bitwise XOR operation.
// \param rhs The right-hand side dense vector for the bitwise XOR operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise XOR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator^=( DenseVector<VT1,TF>& lhs, const DenseVector<VT2,TF>& rhs )
{
   if( IsRestricted_v<VT1> ) {
      if( !tryBitxorAssign( *lhs, *rhs, 0UL ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid bitwise XOR of restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *lhs ) );

   smpAssign( left, left ^ (*rhs) );

   BLAZE_INTERNAL_ASSERT( isIntact( *lhs ), "Invariant violation detected" );

   return *lhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Bitwise XOR assignment operator for the bitwise XOR of a temporary dense vector.
// \ingroup dense_vector
//
// \param lhs The left-hand side temporary dense vector for the bitwise XOR operation.
// \param rhs The right-hand side dense vector for the bitwise XOR operation.
// \return Reference to the dense vector.
// \exception std::invalid_argument Invalid bitwise XOR of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF >     // Transpose flag
inline VT1& operator^=( DenseVector<VT1,TF>&& lhs, const DenseVector<VT2,TF>& rhs )
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
/*!\name DenseVector functions */
//@{
template< typename VT, bool TF >
bool isnan( const DenseVector<VT,TF>& dv );

template< typename VT, bool TF >
bool isinf( const DenseVector<VT,TF>& dv );

template< typename VT, bool TF >
bool isfinite( const DenseVector<VT,TF>& dv );

template< typename VT, bool TF >
bool isDivisor( const DenseVector<VT,TF>& dv );

template< RelaxationFlag RF, typename VT, bool TF >
bool isUniform( const DenseVector<VT,TF>& dv );

template< RelaxationFlag RF, typename VT, bool TF >
bool isZero( const DenseVector<VT,TF>& dv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense vector for not-a-number elements.
// \ingroup dense_vector
//
// \param dv The dense vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
//
// This function checks the N-dimensional dense vector for not-a-number (NaN) elements. If at
// least one element of the vector is not-a-number, the function returns \a true, otherwise it
// returns \a false.

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
bool isnan( const DenseVector<VT,TF>& dv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return false;

   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   for( size_t i=0UL; i<a.size(); ++i ) {
      if( isnan( a[i] ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense vector for infinite elements.
// \ingroup dense_vector
//
// \param dv The dense vector to be checked for infinite elements.
// \return \a true if at least one element of the vector is infinite, \a false otherwise.
//
// This function checks the N-dimensional dense vector for infinite elements. If at least one
// element of the vector is infinite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isinf( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
bool isinf( const DenseVector<VT,TF>& dv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return false;

   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   for( size_t i=0UL; i<a.size(); ++i ) {
      if( isinf( a[i] ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given dense vector for finite elements.
// \ingroup dense_vector
//
// \param dv The dense vector to be checked for finite elements.
// \return \a true if all elements of the vector are finite, \a false otherwise.
//
// This function checks if all elements of the N-dimensional dense vector are finite elements
// (i.e. normal, subnormal or zero elements, but not infinite or NaN). If all elements of the
// vector are finite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::DynamicVector<double> a;
   // ... Resizing and initialization
   if( isfinite( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
bool isfinite( const DenseVector<VT,TF>& dv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return true;

   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   for( size_t i=0UL; i<a.size(); ++i ) {
      if( !isfinite( a[i] ) ) return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given dense vector is a valid divisor.
// \ingroup dense_vector
//
// \param dv The dense vector to be tested.
// \return \a true in case the given vector is a valid divisor, \a false otherwise.
//
// This function checks if the given dense vector is a valid divisor. If all elements of the
// vector are valid divisors the function returns \a true, if at least one element of the vector
// is not a valid divisor, the function returns \a false.

   \code
   StaticVector<int,3UL> a{ 1, -1, 2 };  // isDivisor( a ) returns true
   StaticVector<int,3UL> b{ 1, -1, 0 };  // isDivisor( b ) returns false
   \endcode
*/
template< typename VT  // Type of the dense vector
        , bool TF >    // Transpose flag
bool isDivisor( const DenseVector<VT,TF>& dv )
{
   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   for( size_t i=0UL; i<a.size(); ++i ) {
      if( !isDivisor( a[i] ) ) return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense vector is a uniform vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be checked.
// \return \a true if the vector is a uniform vector, \a false if not.
//
// This function checks if the given dense vector is a uniform vector. The vector is considered
// to be uniform if all its elements are identical. The following code example demonstrates the
// use of the function:

   \code
   blaze::DynamicVector<int,blaze::columnVector> a, b;
   // ... Initialization
   if( isUniform( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniform<relaxed>( a ) ) { ... }
   \endcode

// It is also possible to check if a vector expression results is a uniform vector:

   \code
   if( isUniform( a + b ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary vector.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename VT        // Type of the dense vector
        , bool TF >          // Transpose flag
bool isUniform( const DenseVector<VT,TF>& dv )
{
   if( IsUniform_v<VT> || (*dv).size() < 2UL )
      return true;

   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   const auto& cmp( a[0UL] );

   for( size_t i=1UL; i<a.size(); ++i ) {
      if( !equal<RF>( a[i], cmp ) )
         return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given dense vector is a zero vector.
// \ingroup dense_vector
//
// \param dv The dense vector to be checked.
// \return \a true if the vector is a zero vector, \a false if not.
//
// This function checks if the given dense vector is a zero vector. The vector is considered to
// be zero if all its elements are zero. The following code example demonstrates the use of the
// function:

   \code
   blaze::DynamicVector<int,blaze::columnVector> a, b;
   // ... Initialization
   if( isZero( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isZero<relaxed>( a ) ) { ... }
   \endcode

// It is also possible to check if a vector expression results is a zero vector:

   \code
   if( isZero( a + b ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary vector.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename VT        // Type of the dense vector
        , bool TF >          // Transpose flag
bool isZero( const DenseVector<VT,TF>& dv )
{
   if( IsZero_v<VT> || (*dv).size() == 0UL )
      return true;

   CompositeType_t<VT> a( *dv );  // Evaluation of the dense vector operand

   for( size_t i=0UL; i<a.size(); ++i ) {
      if( !isZero<RF>( a[i] ) )
         return false;
   }

   return true;
}
//*************************************************************************************************

} // namespace blaze

#endif
