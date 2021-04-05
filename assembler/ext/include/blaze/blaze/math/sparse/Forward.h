//=================================================================================================
/*!
//  \file blaze/math/sparse/Forward.h
//  \brief Header file for all forward declarations for sparse vectors and matrices
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

#ifndef _BLAZE_MATH_SPARSE_FORWARD_H_
#define _BLAZE_MATH_SPARSE_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/GroupTag.h>
#include <blaze/system/StorageOrder.h>
#include <blaze/system/TransposeFlag.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< typename Type                   // Data type of the vector
        , bool TF = defaultTransposeFlag  // Transpose flag
        , typename Tag = Group0 >         // Type tag
class CompressedVector;

template< typename Type                  // Data type of the matrix
        , bool SO = defaultStorageOrder  // Storage order
        , typename Tag = Group0 >        // Type tag
class CompressedMatrix;

template< typename Type                  // Data type of the matrix
        , bool SO = defaultStorageOrder  // Storage order
        , typename Tag = Group0 >        // Type tag
class IdentityMatrix;

template< typename Type                   // Data type of the vector
        , bool TF = defaultTransposeFlag  // Transpose flag
        , typename Tag = Group0 >         // Type tag
class ZeroVector;

template< typename Type                  // Data type of the matrix
        , bool SO = defaultStorageOrder  // Storage order
        , typename Tag = Group0 >        // Type tag
class ZeroMatrix;

} // namespace blaze

#endif
