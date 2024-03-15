//=================================================================================================
/*!
//  \file blaze/math/views/subvector/BaseTemplate.h
//  \brief Header file for the implementation of the Subvector base template
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_BASETEMPLATE_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_BASETEMPLATE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Base template of the Subvector class template.
// \ingroup subvector
*/
template< typename VT                    // Type of the vector
        , AlignmentFlag AF = unaligned   // Alignment flag
        , bool TF = IsRowVector_v<VT>    // Transpose flag
        , bool DF = IsDenseVector_v<VT>  // Density flag
        , size_t... CSAs >               // Compile time subvector arguments
class Subvector
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
/*!\brief Auxiliary alias declaration for the Subvector class template.
// \ingroup subvector
//
// The Subvector_ alias declaration represents a convenient shortcut for the specification of the
// non-derived template arguments of the Subvector class template.
*/
template< typename VT                   // Type of the vector
        , AlignmentFlag AF = unaligned  // Alignment flag
        , size_t... CSAs >              // Compile time subvector arguments
using Subvector_ = Subvector< VT
                            , AF
                            , IsRowVector_v<VT>
                            , IsDenseVector_v<VT>
                            , CSAs... >;
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
