//=================================================================================================
/*!
//  \file blaze/math/smp/ThreadMapping.h
//  \brief Header file for the SMP thread mapping functionality
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

#ifndef _BLAZE_MATH_SMP_THREADMAPPING_H_
#define _BLAZE_MATH_SMP_THREADMAPPING_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/expressions/Matrix.h>
#include <blaze/math/shims/Round.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  THREADMAPPING FUNCTIONALITY
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
using ThreadMapping = std::pair<size_t,size_t>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Creates a 2D mapping of threads.
// \ingroup smp
//
// \param threads The total number of threads to be mapped.
// \param A The matrix the mapping is created for.
// \return 2D mapping of the given number of threads.
//
// This function creates a 2D mapping of the given number of threads for the given matrix \a A.
// The mapping will depend on the ratio between rows and columns of the matrix and its storage
// order.
*/
template< typename MT  // Type of the matrix
        , bool SO >    // Storage order of the matrix
ThreadMapping createThreadMapping( size_t threads, const Matrix<MT,SO>& A )
{
   const size_t M( (*A).rows()    );
   const size_t N( (*A).columns() );

   if( M > N || ( M == N && !SO ) )
   {
      const double ratio( double(M)/double(N) );
      size_t m = min( threads, max( 1UL, static_cast<size_t>( round( sqrt( threads*ratio ) ) ) ) );
      size_t n = threads / m;

      while( m * n != threads ) {
         ++m;
         n = threads / m;
      }

      return ThreadMapping( m, n );
   }
   else
   {
      const double ratio( double(N)/double(M) );
      size_t n = min( threads, max( 1UL, static_cast<size_t>( round( sqrt( threads*ratio ) ) ) ) );
      size_t m = threads / n;

      while( m * n != threads ) {
         ++n;
         m = threads / n;
      }

      return ThreadMapping( m, n );
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
