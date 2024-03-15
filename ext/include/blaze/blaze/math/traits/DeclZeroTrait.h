//=================================================================================================
/*!
//  \file blaze/math/traits/DeclZeroTrait.h
//  \brief Header file for the declzero trait
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

#ifndef _BLAZE_MATH_TRAITS_DECLZEROTRAIT_H_
#define _BLAZE_MATH_TRAITS_DECLZEROTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/math/typetraits/StorageOrder.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename > struct DeclZeroTrait;
template< typename, typename = void > struct DeclZeroTraitEval;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T >
auto evalDeclZeroTrait( const volatile T& ) -> DeclZeroTraitEval<T>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the DeclZeroTrait class.
// \ingroup math_traits
//
// \section declzerotrait_general General
//
// The DeclZeroTrait class template offers the possibility to select the resulting data type
// of a generic declzero() operation on the given type \a T. In case the given type \a T is a
// fitting data type, DeclZeroTrait defines the nested type \a Type, which represents the
// resulting data type of the declzero() operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section declzerotrait_specializations Creating custom specializations
//
// Per default, DeclZeroTrait supports all vector and matrix types of the Blaze library (including
// views and adaptors). For all other data types it is possible to specialize the DeclZeroTrait
// template. The following example shows the according specialization for the SymmetricMatrix
// class template:

   \code
   template< typename MT, bool SO, bool DF, bool SF >
   struct DeclZeroTrait< SymmetricMatrix<MT,SO,DF,SF> >
   {
      using Type = ZeroMatrix< ElementType_t<MT>, SO >;
   };
   \endcode

// \n \section declzerotrait_examples Examples
//
// The following example demonstrates the use of the DeclZeroTrait template, where depending on
// the given vector or matrix type the resulting type is selected:

   \code
   using blaze::DynamicVector;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::DeclZeroTrait;
   using blaze::rowVector;
   using blaze::columnMajor;

   // Definition of the resulting type of a row-major dynamic matrix
   using VectorType    = DynamicVector<int,rowVector>;
   using DeclZeroType1 = typename DeclZeroTrait<VectorType>::Type;

   // Definition of the resulting type of a symmetric column-major static matrix
   using MatrixType    = SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> >;
   using DeclZeroType2 = typename DeclZeroTrait<MatrixType>::Type;
   \endcode
*/
template< typename T >  // Type of the vector or matrix
struct DeclZeroTrait
   : public decltype( evalDeclZeroTrait( std::declval<T&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DeclZeroTrait type trait.
// \ingroup math_traits
//
// The DeclZeroTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the DeclZeroTrait class template. For instance, given the vector or matrix type
// \a T the following two type definitions are identical:

   \code
   using Type1 = typename blaze::DeclZeroTrait<T>::Type;
   using Type2 = blaze::DeclZeroTrait_t<T>;
   \endcode
*/
template< typename T >  // Type of the matrix
using DeclZeroTrait_t = typename DeclZeroTrait<T>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the DeclZeroTrait type trait.
// \ingroup math_traits
*/
template< typename T  // Type of the vector or matrix
        , typename >  // Restricting condition
struct DeclZeroTraitEval
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DeclZeroTraitEval class template for vector types.
// \ingroup math_traits
*/
template< typename T >  // Type of the vector
struct DeclZeroTraitEval< T
                        , EnableIf_t< IsVector_v<T> > >
{
   using Type = ZeroVector< typename T::ElementType, TransposeFlag_v<T> >;
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DeclZeroTraitEval class template for matrix types.
// \ingroup math_traits
*/
template< typename T >  // Type of the matrix
struct DeclZeroTraitEval< T
                        , EnableIf_t< IsMatrix_v<T> > >
{
   using Type = ZeroMatrix< typename T::ElementType, StorageOrder_v<T> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
