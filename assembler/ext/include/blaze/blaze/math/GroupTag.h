//=================================================================================================
/*!
//  \file blaze/math/GroupTag.h
//  \brief Header file for the GroupTag class template
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

#ifndef _BLAZE_MATH_GROUPTAG_H_
#define _BLAZE_MATH_GROUPTAG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//*************************************************************************************************
/*!\brief Group tag for vectors and matrices.
// \ingroup math
//
// \section grouptag_general General
//
// Via the GroupTag class template it is possible to define distinct groups of vectors and matrices.
// Only vectors and matrices that belong to the same group can be used together. The attempt to
// combine vectors and matrices from different groups results in a compilation error:

   \code
   using blaze::columnVector;
   using blaze::Group0;
   using blaze::Group1;

   blaze::DynamicVector<int,columnVector,Group0> a0, b0;
   blaze::DynamicVector<int,columnVector,Group1> a1, b1;

   a0 + b0;  // Compiles, a0 and b0 are in the same group (Group0)
   a1 + b1;  // Compiles, a1 and b1 are in the same group (Group1)
   a0 + b1;  // Compilation error: a0 and b1 are not in the same group
   \endcode

// \section grouptag_custom_group Custom Group Tags
//
// The \b Blaze library provides 10 different groups based on the GroupTag class template (Group0
// through Group9). In order to create further groups, it is possible to instantiate the GroupTag
// class template with IDs beyond 9:

   \code
   using Group10 = blaze::GroupTag<10>;
   \endcode
*/
template< size_t ID >
struct GroupTag
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 0. This is the default tag for vectors and matrices.
// \ingroup math
*/
using Group0 = GroupTag<0UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 1.
// \ingroup math
*/
using Group1 = GroupTag<1UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 2.
// \ingroup math
*/
using Group2 = GroupTag<2UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 3.
// \ingroup math
*/
using Group3 = GroupTag<3UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 4.
// \ingroup math
*/
using Group4 = GroupTag<4UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 5.
// \ingroup math
*/
using Group5 = GroupTag<5UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 6.
// \ingroup math
*/
using Group6 = GroupTag<6UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 7.
// \ingroup math
*/
using Group7 = GroupTag<7UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 8.
// \ingroup math
*/
using Group8 = GroupTag<8UL>;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Tag for group 9.
// \ingroup math
*/
using Group9 = GroupTag<9UL>;
//*************************************************************************************************

} // namespace blaze

#endif
