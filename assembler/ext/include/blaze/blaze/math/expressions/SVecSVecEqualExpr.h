//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecSVecEqualExpr.h
//  \brief Header file for the sparse vector/sparse vector equality comparison expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECSVECEQUALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECSVECEQUALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/system/MacroDisable.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL BINARY RELATIONAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two sparse vectors.
// \ingroup sparse_vector
//
// \param a The left-hand side sparse vector for the comparison.
// \param b The right-hand side sparse vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two sparse vectors. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point vectors with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename VT1       // Type of the left-hand side sparse vector
        , bool TF1           // Transpose flag of the left-hand side sparse vector
        , typename VT2       // Type of the right-hand side sparse vector
        , bool TF2 >         // Transpose flag of the right-hand side sparse vector
inline bool equal( const SparseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   using CT1 = CompositeType_t<VT1>;
   using CT2 = CompositeType_t<VT2>;

   // Early exit in case the vector sizes don't match
   if( (*lhs).size() != (*rhs).size() ) return false;

   // Evaluation of the two sparse vector operands
   CT1 a( *lhs );
   CT2 b( *rhs );

   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   const auto lend( a.end() );
   const auto rend( b.end() );

   auto lelem( a.begin() );
   auto relem( b.begin() );

   while( lelem != lend && relem != rend )
   {
      if( isDefault<RF>( lelem->value() ) ) { ++lelem; continue; }
      if( isDefault<RF>( relem->value() ) ) { ++relem; continue; }

      if( lelem->index() != relem->index() || !equal<RF>( lelem->value(), relem->value() ) ) {
         return false;
      }
      else {
         ++lelem;
         ++relem;
      }
   }

   while( lelem != lend ) {
      if( !isDefault<RF>( lelem->value() ) )
         return false;
      ++lelem;
   }

   while( relem != rend ) {
      if( !isDefault<RF>( relem->value() ) )
         return false;
      ++relem;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two sparse vectors.
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the comparison.
// \param rhs The right-hand side sparse vector for the comparison.
// \return \a true if the two sparse vectors are equal, \a false if not.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , bool TF1      // Transpose flag of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
inline bool operator==( const SparseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   return equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two sparse vectors.
// \ingroup sparse_vector
//
// \param lhs The left-hand side sparse vector for the comparison.
// \param rhs The right-hand side sparse vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename VT1  // Type of the left-hand side sparse vector
        , bool TF1      // Transpose flag of the left-hand side sparse vector
        , typename VT2  // Type of the right-hand side sparse vector
        , bool TF2 >    // Transpose flag of the right-hand side sparse vector
inline bool operator!=( const SparseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs )
{
   return !equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
