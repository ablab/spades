//=================================================================================================
/*!
//  \file blaze/math/views/Forward.h
//  \brief Header file for all forward declarations for views
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

#ifndef _BLAZE_MATH_VIEWS_FORWARD_H_
#define _BLAZE_MATH_VIEWS_FORWARD_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/views/band/BaseTemplate.h>
#include <blaze/math/views/column/BaseTemplate.h>
#include <blaze/math/views/columns/BaseTemplate.h>
#include <blaze/math/views/elements/BaseTemplate.h>
#include <blaze/math/views/row/BaseTemplate.h>
#include <blaze/math/views/rows/BaseTemplate.h>
#include <blaze/math/views/submatrix/BaseTemplate.h>
#include <blaze/math/views/subvector/BaseTemplate.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ::blaze NAMESPACE FORWARD DECLARATIONS
//
//=================================================================================================

template< size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&, RSAs... );

template< size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( const Vector<VT,TF>&, RSAs... );

template< size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( const Vector<VT,TF>&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t N, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&&, RSAs... );

template< typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&, size_t, size_t, RSAs... );

template< typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( const Vector<VT,TF>&, size_t, size_t, RSAs... );

template< typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&&, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( const Vector<VT,TF>&, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename VT, bool TF, typename... RSAs >
decltype(auto) subvector( Vector<VT,TF>&&, size_t, size_t, RSAs... );


template< size_t I, size_t... Is, typename VT, bool TF, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&, REAs... );

template< size_t I, size_t... Is, typename VT, bool TF, typename... REAs >
decltype(auto) elements( const Vector<VT,TF>&, REAs... );

template< size_t I, size_t... Is, typename VT, bool TF, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&&, REAs... );

template< typename VT, bool TF, typename T, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&, T*, size_t, REAs... );

template< typename VT, bool TF, typename T, typename... REAs >
decltype(auto) elements( const Vector<VT,TF>&, T*, size_t, REAs... );

template< typename VT, bool TF, typename T, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&&, T*, size_t, REAs... );

template< typename VT, bool TF, typename P, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&, P, size_t, REAs... );

template< typename VT, bool TF, typename P, typename... REAs >
decltype(auto) elements( const Vector<VT,TF>&, P, size_t, REAs... );

template< typename VT, bool TF, typename P, typename... REAs >
decltype(auto) elements( Vector<VT,TF>&&, P, size_t, REAs... );


template< size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&, RSAs... );

template< size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( const Matrix<MT,SO>&, RSAs... );

template< size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( const Matrix<MT,SO>&, RSAs... );

template< AlignmentFlag AF, size_t I, size_t J, size_t M, size_t N, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&&, RSAs... );

template< typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&, size_t, size_t, size_t, size_t, RSAs... );

template< typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( const Matrix<MT,SO>&, size_t, size_t, size_t, size_t, RSAs... );

template< typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&&, size_t, size_t, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&, size_t, size_t, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( const Matrix<MT,SO>&, size_t, size_t, size_t, size_t, RSAs... );

template< AlignmentFlag AF, typename MT, bool SO, typename... RSAs >
decltype(auto) submatrix( Matrix<MT,SO>&&, size_t, size_t, size_t, size_t, RSAs... );


template< size_t I, typename MT, bool SO, typename... RRAs >
decltype(auto) row( Matrix<MT,SO>&, RRAs... );

template< size_t I, typename MT, bool SO, typename... RRAs >
decltype(auto) row( const Matrix<MT,SO>&, RRAs... );

template< size_t I, typename MT, bool SO, typename... RRAs >
decltype(auto) row( Matrix<MT,SO>&&, RRAs... );

template< typename MT, bool SO, typename... RRAs >
decltype(auto) row( Matrix<MT,SO>&, size_t, RRAs... );

template< typename MT, bool SO, typename... RRAs >
decltype(auto) row( const Matrix<MT,SO>&, size_t, RRAs... );

template< typename MT, bool SO, typename... RRAs >
decltype(auto) row( Matrix<MT,SO>&&, size_t, RRAs... );


template< size_t I, size_t... Is, typename MT, bool SO, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&, RRAs... );

template< size_t I, size_t... Is, typename MT, bool SO, typename... RRAs >
decltype(auto) rows( const Matrix<MT,SO>&, RRAs... );

template< size_t I, size_t... Is, typename MT, bool SO, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&&, RRAs... );

template< typename MT, bool SO, typename T, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&, T*, size_t, RRAs... );

template< typename MT, bool SO, typename T, typename... RRAs >
decltype(auto) rows( const Matrix<MT,SO>&, T*, size_t, RRAs... );

template< typename MT, bool SO, typename T, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&&, T*, size_t, RRAs... );

template< typename MT, bool SO, typename P, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&, P, size_t, RRAs... );

template< typename MT, bool SO, typename P, typename... RRAs >
decltype(auto) rows( const Matrix<MT,SO>&, P, size_t, RRAs... );

template< typename MT, bool SO, typename P, typename... RRAs >
decltype(auto) rows( Matrix<MT,SO>&&, P, size_t, RRAs... );


template< size_t I, typename MT, bool SO, typename... RCAs >
decltype(auto) column( Matrix<MT,SO>&, RCAs... );

template< size_t I, typename MT, bool SO, typename... RCAs >
decltype(auto) column( const Matrix<MT,SO>&, RCAs... );

template< size_t I, typename MT, bool SO, typename... RCAs >
decltype(auto) column( Matrix<MT,SO>&&, RCAs... );

template< typename MT, bool SO, typename... RCAs >
decltype(auto) column( Matrix<MT,SO>&, size_t, RCAs... );

template< typename MT, bool SO, typename... RCAs >
decltype(auto) column( const Matrix<MT,SO>&, size_t, RCAs... );

template< typename MT, bool SO, typename... RCAs >
decltype(auto) column( Matrix<MT,SO>&&, size_t, RCAs... );


template< size_t I, size_t... Is, typename MT, bool SO, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&, RCAs... );

template< size_t I, size_t... Is, typename MT, bool SO, typename... RCAs >
decltype(auto) columns( const Matrix<MT,SO>&, RCAs... );

template< size_t I, size_t... Is, typename MT, bool SO, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&&, RCAs... );

template< typename MT, bool SO, typename T, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&, T*, size_t, RCAs... );

template< typename MT, bool SO, typename T, typename... RCAs >
decltype(auto) columns( const Matrix<MT,SO>&, T*, size_t, RCAs... );

template< typename MT, bool SO, typename T, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&&, T*, size_t, RCAs... );

template< typename MT, bool SO, typename P, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&, P, size_t, RCAs... );

template< typename MT, bool SO, typename P, typename... RCAs >
decltype(auto) columns( const Matrix<MT,SO>&, P, size_t, RCAs... );

template< typename MT, bool SO, typename P, typename... RCAs >
decltype(auto) columns( Matrix<MT,SO>&&, P, size_t, RCAs... );


template< ptrdiff_t I, typename MT, bool SO, typename... RBAs >
decltype(auto) band( Matrix<MT,SO>&, RBAs... );

template< ptrdiff_t I, typename MT, bool SO, typename... RBAs >
decltype(auto) band( const Matrix<MT,SO>&, RBAs... );

template< ptrdiff_t I, typename MT, bool SO, typename... RBAs >
decltype(auto) band( Matrix<MT,SO>&&, RBAs... );

template< typename MT, bool SO, typename... RBAs >
decltype(auto) band( Matrix<MT,SO>&, ptrdiff_t, RBAs... );

template< typename MT, bool SO, typename... RBAs >
decltype(auto) band( const Matrix<MT,SO>&, ptrdiff_t, RBAs... );

template< typename MT, bool SO, typename... RBAs >
decltype(auto) band( Matrix<MT,SO>&&, ptrdiff_t, RBAs... );

template< typename MT, bool SO, typename... RDAs >
decltype(auto) diagonal( Matrix<MT,SO>&, RDAs... );

template< typename MT, bool SO, typename... RDAs >
decltype(auto) diagonal( const Matrix<MT,SO>&, RDAs... );

template< typename MT, bool SO, typename... RDAs >
decltype(auto) diagonal( Matrix<MT,SO>&&, RDAs... );

} // namespace blaze

#endif
