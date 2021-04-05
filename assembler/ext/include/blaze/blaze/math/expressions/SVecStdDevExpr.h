//=================================================================================================
/*!
//  \file blaze/math/expressions/SVecStdDevExpr.h
//  \brief Header file for the sparse vector standard deviation expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_SVECSTDDEVEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_SVECSTDDEVEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/SVecMapExpr.h>
#include <blaze/math/expressions/SVecVarExpr.h>
#include <blaze/util/FunctionTrace.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Computes the standard deviation for the given sparse vector.
// \ingroup sparse_vector
//
// \param sv The given sparse vector for the standard deviation computation.
// \return The standard deviation of the given vector.
// \exception std::invalid_argument Invalid input vector.
//
// This function computes the
// <a href="https://en.wikipedia.org/wiki/Standard_deviation">standard deviation</a> for the
// given sparse vector \a sv. Both the non-zero and zero elements of the sparse vector are taken
// into account. Example:

   \code
   using blaze::CompressedVector;

   CompressedVector<int> v{ 1, 4, 3, 6, 7 };

   const double m = stddev( v );  // Results in 2.38747
   \endcode

// In case the size of the given vector is smaller than 2, a \a std::invalid_argument is thrown.
*/
template< typename VT  // Type of the sparse vector
        , bool TF >    // Transpose flag
decltype(auto) stddev( const SparseVector<VT,TF>& sv )
{
   BLAZE_FUNCTION_TRACE;

   return sqrt( var( *sv ) );
}
//*************************************************************************************************

} // namespace blaze

#endif
