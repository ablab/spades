//=================================================================================================
/*!
//  \file blaze/math/expressions/DVecSVecInnerExpr.h
//  \brief Header file for the dense vector/sparse vector inner product expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DVECSVECINNEREXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DVECSVECINNEREXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/ColumnVector.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RowVector.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL BINARY ARITHMETIC OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (inner product) of a dense and a
//        sparse vector (\f$ s=\vec{a}*\vec{b} \f$).
// \ingroup sparse_vector
//
// \param lhs The left-hand side dense vector for the inner product.
// \param rhs The right-hand side sparse vector for the inner product.
// \return The scalar product.
// \exception std::invalid_argument Vector sizes do not match.
//
// This operator represents the scalar product (inner product) of a dense vector and a sparse
// vector:

   \code
   using blaze::columnVector;
   using blaze::rowVector;

   blaze::DynamicVector<double,rowVector> a;
   blaze::CompressedVector<double,columnVector> b;
   blaze::real res;
   // ... Resizing and initialization
   res = a * b;
   \endcode

// The operator returns a scalar value of the higher-order element type of the two involved
// vector element types \a VT1::ElementType and \a VT2::ElementType. Both vector types \a VT1
// and \a VT2 as well as the two element types \a VT1::ElementType and \a VT2::ElementType
// have to be supported by the MultTrait class template.\n
// In case the current sizes of the two given vectors don't match, a \a std::invalid_argument
// is thrown.
*/
template< typename VT1    // Type of the left-hand side dense vector
        , typename VT2 >  // Type of the right-hand side sparse vector
inline decltype(auto)
   operator*( const DenseVector<VT1,true>& lhs, const SparseVector<VT2,false>& rhs )
{
   BLAZE_FUNCTION_TRACE;

   using CT1      = CompositeType_t<VT1>;    // Composite type of the left-hand side dense vector expression
   using CT2      = CompositeType_t<VT2>;    // Composite type of the right-hand side sparse vector expression
   using XT1      = RemoveReference_t<CT1>;  // Auxiliary type for the left-hand side composite type
   using XT2      = RemoveReference_t<CT2>;  // Auxiliary type for the right-hand side composite type
   using ET1      = ElementType_t<XT1>;      // Element type of the left-hand side dense vector expression
   using ET2      = ElementType_t<XT2>;      // Element type of the right-hand side sparse vector expression
   using MultType = MultTrait_t<ET1,ET2>;    // Multiplication result type

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( VT2 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_VECTOR_TYPE   ( VT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_VECTOR_TYPE( VT2 );

   if( (*lhs).size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   CT1 left ( *lhs );
   CT2 right( *rhs );

   auto element( right.begin() );
   auto end    ( right.end()   );

   MultType sp{};

   if( element != end ) {
      sp = left[ element->index() ] * element->value();
      ++element;
      for( ; element!=end; ++element )
         sp += left[ element->index() ] * element->value();
   }

   return sp;
}
//*************************************************************************************************

} // namespace blaze

#endif
