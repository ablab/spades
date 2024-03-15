//=================================================================================================
/*!
//  \file blaze/math/views/band/BaseTemplate.h
//  \brief Header file for the implementation of the Band base template
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

#ifndef _BLAZE_MATH_VIEWS_BAND_BASETEMPLATE_H_
#define _BLAZE_MATH_VIEWS_BAND_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsMatMatMultExpr.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base template of the Band class template.
// \ingroup band
*/
template< typename MT                       // Type of the matrix
        , bool TF = defaultTransposeFlag    // Transpose flag
        , bool DF = IsDenseMatrix_v<MT>     // Density flag
        , bool MF = IsMatMatMultExpr_v<MT>  // Multiplication flag
        , ptrdiff_t... CBAs >               // Compile time band arguments
class Band
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ALIAS DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary alias declaration for the Band class template.
// \ingroup band
//
// The Band_ alias declaration represents a convenient shortcut for the specification of the
// non-derived template arguments of the Band class template.
*/
template< typename MT          // Type of the matrix
        , ptrdiff_t... CBAs >  // Compile time band arguments
using Band_ = Band< MT
                  , defaultTransposeFlag
                  , IsDenseMatrix_v<MT>
                  , IsMatMatMultExpr_v<MT>
                  , CBAs... >;
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
