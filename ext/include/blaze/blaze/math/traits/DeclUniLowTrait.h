//=================================================================================================
/*!
//  \file blaze/math/traits/DeclUniLowTrait.h
//  \brief Header file for the declunilow trait
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

#ifndef _BLAZE_MATH_TRAITS_DECLUNILOWTRAIT_H_
#define _BLAZE_MATH_TRAITS_DECLUNILOWTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/adaptors/unilowermatrix/BaseTemplate.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/util/EnableIf.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename > struct DeclUniLowTrait;
template< typename, typename = void > struct DeclUniLowTraitEval;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T >
auto evalDeclUniLowTrait( const volatile T& ) -> DeclUniLowTraitEval<T>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the DeclUniLowTrait class.
// \ingroup math_traits
//
// \section declunilowtrait_general General
//
// The DeclUniLowTrait class template offers the possibility to select the resulting data type
// of a generic declunilow() operation on the given type \a MT. In case the given type \a MT
// is a dense or sparse matrix type, DeclUniLowTrait defines the nested type \a Type, which
// represents the resulting data type of the declunilow() operation. Otherwise there is no
// nested type \a Type. Note that \a const and \a volatile qualifiers and reference modifiers
// are generally ignored.
//
//
// \section declunilowtrait_specializations Creating custom specializations
//
// Per default, DeclUniLowTrait supports all matrix types of the Blaze library (including views
// and adaptors). For all other data types it is possible to specialize the DeclUniLowTrait
// template. The following example shows the according specialization for the SymmetricMatrix
// class template:

   \code
   template< typename MT, bool SO, bool DF, bool SF >
   struct DeclUniLowTrait< SymmetricMatrix<MT,SO,DF,SF> >
   {
      using Type = IdentityMatrix<ElementType_t<MT>,SO>;
   };
   \endcode

// \n \section declunilowtrait_examples Examples
//
// The following example demonstrates the use of the DeclUniLowTrait template, where depending
// on the given matrix type the resulting type is selected:

   \code
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::DeclUniLowTrait;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the resulting type of a row-major dynamic matrix
   using MatrixType1     = DynamicMatrix<int,rowMajor>;
   using DeclUniLowType1 = typename DeclUniLowTrait<MatrixType1>::Type;

   // Definition of the resulting type of a symmetric column-major static matrix
   using MatrixType2     = SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> >;
   using DeclUniLowType2 = typename DeclUniLowTrait<MatrixType2>::Type;
   \endcode
*/
template< typename MT >  // Type of the matrix
struct DeclUniLowTrait
   : public decltype( evalDeclUniLowTrait( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DeclUniLowTrait type trait.
// \ingroup math_traits
//
// The DeclUniLowTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the DeclUniLowTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::DeclUniLowTrait<MT>::Type;
   using Type2 = blaze::DeclUniLowTrait_t<MT>;
   \endcode
*/
template< typename MT >  // Type of the matrix
using DeclUniLowTrait_t = typename DeclUniLowTrait<MT>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the DeclUniLowTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , typename >   // Restricting condition
struct DeclUniLowTraitEval
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DeclUniLowTraitEval class template for square matrix types.
// \ingroup math_traits
*/
template< typename MT >  // Type of the matrix
struct DeclUniLowTraitEval< MT
                          , EnableIf_t< IsMatrix_v<MT> &&
                                        ( Size_v<MT,0UL> == DefaultSize_v ||
                                          Size_v<MT,1UL> == DefaultSize_v ||
                                          Size_v<MT,0UL> == Size_v<MT,1UL> ) > >
{
   using Type = UniLowerMatrix<typename MT::ResultType>;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
