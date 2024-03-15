//=================================================================================================
/*!
//  \file blaze/math/typetraits/IsSparseMatrix.h
//  \brief Header file for the IsSparseMatrix type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_ISSPARSEMATRIX_H_
#define _BLAZE_MATH_TYPETRAITS_ISSPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/expressions/Forward.h>
#include <blaze/util/IntegralConstant.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper functions for the IsSparseMatrix type trait.
// \ingroup math_type_traits
*/
template< typename MT, bool SO >
TrueType isSparseMatrix_backend( const volatile SparseMatrix<MT,SO>* );

FalseType isSparseMatrix_backend( ... );
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Compile time check for sparse matrix types.
// \ingroup math_type_traits
//
// This type trait tests whether or not the given template parameter is a sparse matrix type
// (i.e. whether \a T is derived from the SparseMatrix base class). In case the type is a sparse
// matrix type, the \a value member constant is set to \a true, the nested type definition
// \a Type is \a TrueType, and the class derives from \a TrueType. Otherwise \a value is set
// to \a false, \a Type is \a FalseType, and the class derives from \a FalseType.

   \code
   using namespace blaze;

   IsSparseMatrix< CompressedMatrix<double> >::value      // Evaluates to 1
   IsSparseMatrix< const CompressedMatrix<float> >::Type  // Results in TrueType
   IsSparseMatrix< volatile CompressedMatrix<int> >       // Is derived from TrueType
   IsSparseMatrix< DynamicVector<double> >::value         // Evaluates to 0
   IsSparseMatrix< const DynamicMatrix<double> >::Type    // Results in FalseType
   IsSparseMatrix< CompressedVector<double> >             // Is derived from FalseType
   \endcode
*/
template< typename T >
struct IsSparseMatrix
   : public decltype( isSparseMatrix_backend( std::declval<T*>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the IsSparseMatrix type trait for references.
// \ingroup math_type_traits
*/
template< typename T >
struct IsSparseMatrix<T&>
   : public FalseType
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary variable template for the IsSparseMatrix type trait.
// \ingroup math_type_traits
//
// The IsSparseMatrix_v variable template provides a convenient shortcut to access the nested
// \a value of the IsSparseMatrix class template. For instance, given the type \a T the
// following two statements are identical:

   \code
   constexpr bool value1 = blaze::IsSparseMatrix<T>::value;
   constexpr bool value2 = blaze::IsSparseMatrix_v<T>;
   \endcode
*/
template< typename T >
constexpr bool IsSparseMatrix_v = IsSparseMatrix<T>::value;
//*************************************************************************************************

} // namespace blaze

#endif
