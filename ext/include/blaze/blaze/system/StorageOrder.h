//=================================================================================================
/*!
//  \file blaze/system/StorageOrder.h
//  \brief Header file for the default storage order for all vectors of the Blaze library
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

#ifndef _BLAZE_SYSTEM_STORAGEORDER_H_
#define _BLAZE_SYSTEM_STORAGEORDER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/config/StorageOrder.h>
#include <blaze/math/StorageOrder.h>


namespace blaze {

//=================================================================================================
//
//  STORAGE ORDER
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default storage order for all matrices of the Blaze library.
// \ingroup system
//
// This value specifies the default storage order for all matrices of the Blaze library.
// In case no explicit storage order is specified with the according matrix type, this
// setting is used.

   \code
   // Explicit specification of the storage order => row-major matrix
   StaticMatrix<double,3UL,3UL,rowMajor> A;

   // No explicit specification of the storage order => use of the default storage order
   StaticMatrix<double,3UL,3UL> B;
   \endcode

// The default storage order is defined via the BLAZE_DEFAULT_STORAGE_ORDER compilation switch
// (see the \ref storage_order section). Valid settings for this value are blaze::rowMajor and
// blaze::columnMajor.
*/
constexpr bool defaultStorageOrder = BLAZE_DEFAULT_STORAGE_ORDER;
//*************************************************************************************************

} // namespace blaze

#endif
