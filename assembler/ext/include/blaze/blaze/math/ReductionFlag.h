//=================================================================================================
/*!
//  \file blaze/math/ReductionFlag.h
//  \brief Header file for the reduction flags
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

#ifndef _BLAZE_MATH_REDUCTIONFLAG_H_
#define _BLAZE_MATH_REDUCTIONFLAG_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  REDUCTION FLAGS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Type of the reduction flags.
*/
using ReductionFlag = size_t;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduction flag for row-wise reduction operations.
//
// This flag can be used to perform row-wise reduction operations on matrices. The following
// example shows the row-wise summation of a row-major matrix:

   \code
   using blaze::rowMajor;
   using blaze::columnVector;

   blaze::DynamicMatrix<int,rowMajor> A{ { 4, 1, 2 }, { -2, 0, 3 } };
   blaze::DynamicVector<int,columnVector> v;

   v = sum<rowwise>( A );  // Results in ( 7, 1 )
   \endcode
*/
constexpr ReductionFlag rowwise = 1UL;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reduction flag for column-wise reduction operations.
//
// This flag can be used to perform column-wise reduction operations on matrices. The following
// example shows the column-wise summation of a column-major matrix:

   \code
   using blaze::columnMajor;
   using blaze::rowVector;

   blaze::DynamicMatrix<int,columnMajor> A{ { 4, 1, 2 }, { -2, 0, 3 } };
   blaze::DynamicVector<int,rowVector> v;

   v = sum<columnwise>( A );  // Results in ( 2, 1, 5 )
   \endcode
*/
constexpr ReductionFlag columnwise = 0UL;
//*************************************************************************************************

} // namespace blaze

#endif
