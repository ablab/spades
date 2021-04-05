//=================================================================================================
/*!
//  \file blaze/math/typetraits/StorageOrder.h
//  \brief Header file for the StorageOrder type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_STORAGEORDER_H_
#define _BLAZE_MATH_TYPETRAITS_STORAGEORDER_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Evaluation of the storage order of a given matrix type.
// \ingroup math_type_traits
//
// Via this type trait it is possible to evaluate the storage order of a given matrix type.
// In case the given type is a row-major matrix type the nested boolean \a value is set to
// \a rowMajor, in case it is a column-major matrix type it is set to \a columnMajor. If the
// given type is not a matrix type a compilation error is created.

   \code
   using RowMajorMatrix    = blaze::DynamicMatrix<int,blaze::rowMajor>;
   using ColumnMajorMatrix = blaze::DynamicMatrix<int,blaze::columnMajor>;

   blaze::StorageOrder<RowMajorMatrix>::value     // Evaluates to blaze::rowMajor
   blaze::StorageOrder<ColumnMajorMatrix>::value  // Evaluates to blaze::columnMajor
   blaze::StorageOrder<int>::value                // Compilation error!
   \endcode
*/
template< typename T >
struct StorageOrder
   : public BoolConstant< T::storageOrder >
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the StorageOrder type trait.
// \ingroup math_type_traits
//
// The StorageOrder_v variable template provides a convenient shortcut to access the nested
// \a value of the StorageOrder class template. For instance, given the matrix type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::StorageOrder<T>::value;
   constexpr bool value2 = blaze::StorageOrder_v<T>;
   \endcode
*/
template< typename T >
constexpr bool StorageOrder_v = T::storageOrder;
//*************************************************************************************************

} // namespace blaze

#endif
