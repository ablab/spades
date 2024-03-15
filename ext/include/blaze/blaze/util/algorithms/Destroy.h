//=================================================================================================
/*!
//  \file blaze/util/algorithms/Destroy.h
//  \brief Header file for the generic destroy algorithm
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

#ifndef _BLAZE_UTIL_ALGORITHMS_DESTROY_H_
#define _BLAZE_UTIL_ALGORITHMS_DESTROY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/util/algorithms/DestroyAt.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  DESTROY ALGORITHMS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Destroys the given range of objects .
// \ingroup algorithms
//
// \param first Iterator to the first element to be destroyed.
// \param last Iterator to the element one past the last element to be destroyed.
// \return void
//
// This function explicitly calls the destructor of all object in the given range.
*/
template< typename ForwardIt >
void destroy( ForwardIt first, ForwardIt last )
{
   for( ; first!=last; ++first ) {
      blaze::destroy_at( std::addressof( *first ) );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Destroys the given range of objects .
// \ingroup algorithms
//
// \param first Iterator to the first element to be destroyed.
// \param n The number of elements to be destroyed.
// \return void
//
// This function explicitly calls the destructor of all object in the given range.
*/
template< typename ForwardIt >
void destroy_n( ForwardIt first, size_t n )
{
   for( ; n > 0UL; (void) ++first, --n ) {
      blaze::destroy_at( std::addressof( *first ) );
   }
}
//*************************************************************************************************

} // namespace blaze

#endif
