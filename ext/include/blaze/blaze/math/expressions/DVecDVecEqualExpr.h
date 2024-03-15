//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecDVecEqualExpr.h
//  \brief Header file for the dense vector/dense vector equality comparison expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECDVECEQUALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECDVECEQUALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/HasSIMDEqual.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Optimizations.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the dense vector/dense vector equality comparison.
// \ingroup dense_vector
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side dense vector
struct DVecDVecEqualExprHelper
{
   //**Type definitions****************************************************************************
   //! Composite type of the left-hand side dense vector expression.
   using CT1 = RemoveReference_t< CompositeType_t<VT1> >;

   //! Composite type of the right-hand side dense vector expression.
   using CT2 = RemoveReference_t< CompositeType_t<VT2> >;
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool value =
      ( useOptimizedKernels &&
        CT1::simdEnabled &&
        CT2::simdEnabled &&
        HasSIMDEqual_v< ElementType_t<CT1>, ElementType_t<CT2> > );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY RELATIONAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default equality check of two dense vectors.
// \ingroup dense_vector
//
// \param a The left-hand side dense vector for the comparison.
// \param b The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two dense vectors. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point vectors with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename VT1       // Type of the left-hand side dense vector
        , bool TF1           // Transpose flag of the left-hand side dense vector
        , typename VT2       // Type of the right-hand side dense vector
        , bool TF2 >         // Transpose flag of the right-hand side dense vector
inline auto equal( const DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
   -> DisableIf_t< DVecDVecEqualExprHelper<VT1,VT2>::value, bool >
{
   using CT1 = CompositeType_t<VT1>;
   using CT2 = CompositeType_t<VT2>;

   // Early exit in case the vector sizes don't match
   if( (*lhs).size() != (*rhs).size() ) return false;

   // Evaluation of the two dense vector operands
   CT1 a( *lhs );
   CT2 b( *rhs );

   // In order to compare the two vectors, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0UL; i<a.size(); ++i )
      if( !equal<RF>( a[i], b[i] ) ) return false;
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized equality check of two dense vectors.
// \ingroup dense_vector
//
// \param a The left-hand side dense vector for the comparison.
// \param b The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
//
// Equal function for the comparison of two dense vectors. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point vectors with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename VT1       // Type of the left-hand side dense vector
        , bool TF1           // Transpose flag of the left-hand side dense vector
        , typename VT2       // Type of the right-hand side dense vector
        , bool TF2 >         // Transpose flag of the right-hand side dense vector
inline auto equal( const DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
   -> EnableIf_t< DVecDVecEqualExprHelper<VT1,VT2>::value, bool >
{
   using CT1 = CompositeType_t<VT1>;
   using CT2 = CompositeType_t<VT2>;
   using XT1 = RemoveReference_t<CT1>;
   using XT2 = RemoveReference_t<CT2>;

   // Early exit in case the vector sizes don't match
   if( (*lhs).size() != (*rhs).size() ) return false;

   // Evaluation of the two dense vector operands
   CT1 a( *lhs );
   CT2 b( *rhs );

   constexpr size_t SIMDSIZE = SIMDTrait< ElementType_t<VT1> >::size;
   constexpr bool remainder( !IsPadded_v<XT1> || !IsPadded_v<XT2> );

   const size_t N( a.size() );

   const size_t ipos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
   BLAZE_INTERNAL_ASSERT( ipos <= N, "Invalid end calculation" );

   size_t i( 0UL );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      if( !equal<RF>( a.load(i             ), b.load(i             ) ) ) return false;
      if( !equal<RF>( a.load(i+SIMDSIZE    ), b.load(i+SIMDSIZE    ) ) ) return false;
      if( !equal<RF>( a.load(i+SIMDSIZE*2UL), b.load(i+SIMDSIZE*2UL) ) ) return false;
      if( !equal<RF>( a.load(i+SIMDSIZE*3UL), b.load(i+SIMDSIZE*3UL) ) ) return false;
   }
   for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL ) {
      if( !equal<RF>( a.load(i         ), b.load(i         ) ) ) return false;
      if( !equal<RF>( a.load(i+SIMDSIZE), b.load(i+SIMDSIZE) ) ) return false;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      if( !equal<RF>( a.load(i), b.load(i) ) ) return false;
   }
   for( ; remainder && i<N; ++i ) {
      if( !equal<RF>( a[i], b[i] ) ) return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense vectors.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are equal, \a false if not.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
inline bool operator==( const DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   return equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense vectors.
// \ingroup dense_vector
//
// \param lhs The left-hand side dense vector for the comparison.
// \param rhs The right-hand side dense vector for the comparison.
// \return \a true if the two vectors are not equal, \a false if they are equal.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side dense vector
        , bool TF2 >    // Transpose flag of the right-hand side dense vector
inline bool operator!=( const DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs )
{
   return !equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
