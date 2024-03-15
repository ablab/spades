//=================================================================================================
/*!
//  \file blaze/math/sparse/SparseVector.h
//  \brief Header file for utility functions for sparse vectors
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

#ifndef _BLAZE_MATH_SPARSE_SPARSEVECTOR_H_
#define _BLAZE_MATH_SPARSE_SPARSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsFinite.h>
#include <blaze/math/shims/IsInf.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/util/constraints/Numeric.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseVector operators */
//@{
template< typename VT, bool TF, typename ST >
auto operator*=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator*=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator/=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;

template< typename VT, bool TF, typename ST >
auto operator/=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector and
//        a scalar value (\f$ \vec{a}*=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side sparse vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   if( IsRestricted_v<VT> ) {
      if( !tryMult( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   if( !IsResizable_v< ElementType_t<VT> > && isZero( scalar ) )
   {
      reset( *vec );
   }
   else
   {
      decltype(auto) left( derestrict( *vec ) );

      const auto last( left.end() );
      for( auto element=left.begin(); element!=last; ++element ) {
         element->value() *= scalar;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a temporary sparse vector
//        and a scalar (\f$ v*=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side temporary sparse vector for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator*=( *vec, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse vector by a scalar value
//        (\f$ \vec{a}/=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side sparse vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseVector<VT,TF>& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   BLAZE_USER_ASSERT( !isZero( scalar ), "Division by zero detected" );

   if( IsRestricted_v<VT> ) {
      if( !tryDiv( *vec, 0UL, (*vec).size(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted vector" );
      }
   }

   using ScalarType = If_t< IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > ||
                            IsFloatingPoint_v< UnderlyingBuiltin_t<ST> >
                          , If_t< IsComplex_v< UnderlyingScalar_t<VT> > && IsBuiltin_v<ST>
                                , DivTrait_t< UnderlyingBuiltin_t<VT>, ST >
                                , DivTrait_t< UnderlyingScalar_t<VT>, ST > >
                          , ST >;

   decltype(auto) left( derestrict( *vec ) );

   if( IsInvertible_v<ScalarType> ) {
      const ScalarType tmp( ScalarType(1)/static_cast<ScalarType>( scalar ) );
      for( auto element=left.begin(); element!=left.end(); ++element ) {
         element->value() *= tmp;
      }
   }
   else {
      for( auto element=left.begin(); element!=left.end(); ++element ) {
         element->value() /= scalar;
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *vec ), "Invariant violation detected" );

   return *vec;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a temporary sparse vector by a scalar
//        value (\f$ \vec{a}/=s \f$).
// \ingroup sparse_vector
//
// \param vec The left-hand side temporary sparse vector for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse vector.
// \exception std::invalid_argument Invalid scaling of restricted vector.
//
// In case the vector \a VT is restricted and the assignment would violate an invariant of the
// vector, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename VT    // Type of the left-hand side sparse vector
        , bool TF        // Transpose flag
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseVector<VT,TF>&& vec, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, VT& >
{
   return operator/=( *vec, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseVector functions */
//@{
template< typename VT, bool TF >
bool isnan( const SparseVector<VT,TF>& sv );

template< typename VT, bool TF >
bool isinf( const SparseVector<VT,TF>& sv );

template< typename VT, bool TF >
bool isfinite( const SparseVector<VT,TF>& sv );

template< RelaxationFlag RF, typename VT, bool TF >
bool isUniform( const SparseVector<VT,TF>& sv );

template< RelaxationFlag RF, typename VT, bool TF >
bool isZero( const SparseVector<VT,TF>& sv );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse vector for not-a-number elements.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked for not-a-number elements.
// \return \a true if at least one element of the vector is not-a-number, \a false otherwise.
//
// This function checks the N-dimensional sparse vector for not-a-number (NaN) elements. If
// at least one element of the vector is not-a-number, the function returns \a true, otherwise
// it returns \a false.

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isnan( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline bool isnan( const SparseVector<VT,TF>& sv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return false;

   using CT = CompositeType_t<VT>;

   CT a( *sv );  // Evaluation of the sparse vector operand

   const auto end( a.end() );
   for( auto element=a.begin(); element!=end; ++element ) {
      if( isnan( element->value() ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse vector for infinite elements.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked for infinite elements.
// \return \a true if at least one element of the vector is infinite, \a false otherwise.
//
// This function checks the N-dimensional sparse vector for infinite (inf) elements. If at
// least one element of the vector is infinite, the function returns \a true, otherwise
// it returns \a false.

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isinf( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline bool isinf( const SparseVector<VT,TF>& sv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return false;

   using CT = CompositeType_t<VT>;

   CT a( *sv );  // Evaluation of the sparse vector operand

   const auto end( a.end() );
   for( auto element=a.begin(); element!=end; ++element ) {
      if( isinf( element->value() ) ) return true;
   }
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse vector for finite elements.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked for finite elements.
// \return \a true if all elements of the vector are finite, \a false otherwise.
//
// This function checks if all elements of the N-dimensional sparse vector are finite elements
// (i.e. normal, subnormal or zero elements, but not infinite or NaN). If all elements of the
// vector are finite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isfinite( a ) ) { ... }
   \endcode
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
inline bool isfinite( const SparseVector<VT,TF>& sv )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<VT> > )
      return true;

   using CT = CompositeType_t<VT>;

   CT a( *sv );  // Evaluation of the sparse vector operand

   const auto end( a.end() );
   for( auto element=a.begin(); element!=end; ++element ) {
      if( !isfinite( element->value() ) ) return false;
   }
   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse vector is a uniform vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked.
// \return \a true if the vector is a uniform vector, \a false if not.
//
// This function checks if the given sparse vector is a uniform vector. The vector is considered
// to be uniform if all its elements are identical. The following code example demonstrates the
// use of the function:

   \code
   blaze::CompressedVector<int,blaze::columnVector> a, b;
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
        , typename VT        // Type of the sparse vector
        , bool TF >          // Transpose flag
bool isUniform( const SparseVector<VT,TF>& sv )
{
   using CT = CompositeType_t<VT>;

   if( IsUniform_v<VT> || (*sv).size() < 2UL )
      return true;

   CT a( *sv );  // Evaluation of the sparse vector operand

   if( a.nonZeros() != a.size() )
   {
      const auto end( a.end() );
      for( auto element=a.begin(); element!=end; ++element ) {
         if( !isDefault<RF>( element->value() ) )
            return false;
      }
   }
   else
   {
      const auto& cmp( a[0] );
      auto element( a.begin() );
      const auto end( a.end() );

      ++element;

      for( ; element!=end; ++element ) {
         if( !equal<RF>( element->value(), cmp ) )
            return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse vector is a zero vector.
// \ingroup sparse_vector
//
// \param sv The sparse vector to be checked.
// \return \a true if the vector is a zero vector, \a false if not.
//
// This function checks if the given sparse vector is a zero vector. The vector is considered
// to be zero if all its elements are zero. The following code example demonstrates the use
// of the function:

   \code
   blaze::CompressedVector<int,blaze::columnVector> a, b;
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
        , typename VT        // Type of the sparse vector
        , bool TF >          // Transpose flag
bool isZero( const SparseVector<VT,TF>& sv )
{
   if( IsZero_v<VT> || (*sv).nonZeros() == 0UL )
      return true;

   CompositeType_t<VT> a( *sv );  // Evaluation of the sparse vector operand

   const auto end( a.end() );
   for( auto element=a.begin(); element!=end; ++element ) {
      if( !isZero<RF>( element->value() ) )
         return false;
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a single element, a range or selection of elements from the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector.
// \param args The runtime arguments for the erase call.
// \return The result of the according erase member function.
//
// This function represents an abstract interface for erasing a single element, a range of
// elements or a selection of elements from the given sparse vector. It forwards the given
// arguments to the according \a erase() member function of the sparse vector and returns
// the result of the function call.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... Args >  // Type of the erase arguments
auto erase( SparseVector<VT,TF>& sv, Args&&... args )
   -> decltype( (*sv).erase( std::forward<Args>( args )... ) )
{
   return (*sv).erase( std::forward<Args>( args )... );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
