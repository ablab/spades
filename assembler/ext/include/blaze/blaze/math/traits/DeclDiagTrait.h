//=================================================================================================
/*!
//  \file blaze/math/traits/DeclDiagTrait.h
//  \brief Header file for the decldiag trait
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

#ifndef _BLAZE_MATH_TRAITS_DECLDIAGTRAIT_H_
#define _BLAZE_MATH_TRAITS_DECLDIAGTRAIT_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/adaptors/diagonalmatrix/BaseTemplate.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/Size.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename > struct DeclDiagTrait;
template< typename, typename = void > struct DeclDiagTraitEval;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T >
auto evalDeclDiagTrait( const volatile T& ) -> DeclDiagTraitEval<T>;
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the DeclDiagTrait class.
// \ingroup math_traits
//
// \section decldiagtrait_general General
//
// The DeclDiagTrait class template offers the possibility to select the resulting data type
// of a generic decldiag() operation on the given type \a MT. In case the given type \a MT is
// a dense or sparse matrix type, DeclDiagTrait defines the nested type \a Type, which represents
// the resulting data type of the decldiag() operation. Otherwise there is no nested type \a Type.
// Note that \a const and \a volatile qualifiers and reference modifiers are generally ignored.
//
//
// \section decldiagtrait_specializations Creating custom specializations
//
// Per default, DeclDiagTrait supports all matrix types of the Blaze library (including views and
// adaptors). For all other data types it is possible to specialize the DeclDiagTrait template.
// The following example shows the according specialization for the SymmetricMatrix class template:

   \code
   template< typename MT, bool SO, bool DF, bool SF >
   struct DeclDiagTrait< SymmetricMatrix<MT,SO,DF,SF> >
   {
      using Type = DiagonalMatrix<MT>;
   };
   \endcode

// \n \section decldiagtrait_examples Examples
//
// The following example demonstrates the use of the DeclDiagTrait template, where depending on
// the given matrix type the resulting type is selected:

   \code
   using blaze::DynamicMatrix;
   using blaze::StaticMatrix;
   using blaze::SymmetricMatrix;
   using blaze::DeclDiagTrait;
   using blaze::rowMajor;
   using blaze::columnMajor;

   // Definition of the resulting type of a row-major dynamic matrix
   using MatrixType1   = DynamicMatrix<int,rowMajor>;
   using DeclDiagType1 = typename DeclDiagTrait<MatrixType1>::Type;

   // Definition of the resulting type of a symmetric column-major static matrix
   using MatrixType2   = SymmetricMatrix< StaticMatrix<int,3UL,3UL,columnMajor> >;
   using DeclDiagType2 = typename DeclDiagTrait<MatrixType2>::Type;
   \endcode
*/
template< typename MT >  // Type of the matrix
struct DeclDiagTrait
   : public decltype( evalDeclDiagTrait( std::declval<MT&>() ) )
{};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the DeclDiagTrait type trait.
// \ingroup math_traits
//
// The DeclDiagTrait_t alias declaration provides a convenient shortcut to access the nested
// \a Type of the DeclDiagTrait class template. For instance, given the matrix type \a MT the
// following two type definitions are identical:

   \code
   using Type1 = typename blaze::DeclDiagTrait<MT>::Type;
   using Type2 = blaze::DeclDiagTrait_t<MT>;
   \endcode
*/
template< typename MT >  // Type of the matrix
using DeclDiagTrait_t = typename DeclDiagTrait<MT>::Type;
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the DeclDiagTrait type trait.
// \ingroup math_traits
*/
template< typename MT  // Type of the matrix
        , typename >   // Restricting condition
struct DeclDiagTraitEval
{};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the DeclDiagTraitEval class template for square matrix types.
// \ingroup math_traits
*/
template< typename MT >  // Type of the matrix
struct DeclDiagTraitEval< MT
                        , EnableIf_t< IsMatrix_v<MT> &&
                                      ( Size_v<MT,0UL> == DefaultSize_v ||
                                        Size_v<MT,1UL> == DefaultSize_v ||
                                        Size_v<MT,0UL> == Size_v<MT,1UL> ) > >
{
   using Type = DiagonalMatrix<typename MT::ResultType>;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
