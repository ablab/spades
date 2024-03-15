//=================================================================================================
/*!
//  \file blaze/math/typetraits/MakeComplex.h
//  \brief Header file for the MakeComplex type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_MAKECOMPLEX_H_
#define _BLAZE_MATH_TYPETRAITS_MAKECOMPLEX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Converting the given type to the matching 'complex' type.
// \ingroup type_traits
//
// The MakeComplex type trait converts the given floating point type \a T to the matching
// 'complex' type with the same cv-qualifiers. 'complex' types are preserved. For all other
// types, including integral types, no type conversion is performed.

   \code
   blaze::MakeComplex< float                    >::Type  // Results in 'complex<float>'
   blaze::MakeComplex< complex<double>          >::Type  // Results in 'complex<double>'
   blaze::MakeComplex< const long double        >::Type  // Results in 'const complex<long double>'
   blaze::MakeComplex< volatile complex<double> >::Type  // Results in 'volatile complex<double>'
   \endcode

// Note that it is possible to add support for other data types by specializing the MakeComplex
// class template.
*/
template< typename T >
struct MakeComplex
{};
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'float'.
// \ingroup type_traits
*/
template<>
struct MakeComplex<float>
{
   //**********************************************************************************************
   using Type = complex<float>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'double'.
// \ingroup type_traits
*/
template<>
struct MakeComplex<double>
{
   //**********************************************************************************************
   using Type = complex<double>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'long double'.
// \ingroup type_traits
*/
template<>
struct MakeComplex<long double>
{
   //**********************************************************************************************
   using Type = complex<long double>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for complex numbers.
// \ingroup type_traits
*/
template< typename T >
struct MakeComplex< complex<T> >
{
   //**********************************************************************************************
   using Type = complex<T>;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'const T'.
// \ingroup type_traits
*/
template< typename T >
struct MakeComplex<const T>
{
   //**********************************************************************************************
   using Type = const typename MakeComplex<T>::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'volatile T'.
// \ingroup type_traits
*/
template< typename T >
struct MakeComplex<volatile T>
{
   //**********************************************************************************************
   using Type = volatile typename MakeComplex<T>::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*! Specialization of the MakeComplex type trait for 'const volatile T'.
// \ingroup type_traits
*/
template< typename T >
struct MakeComplex<const volatile T>
{
   //**********************************************************************************************
   using Type = const volatile typename MakeComplex<T>::Type;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the MakeComplex type trait.
// \ingroup type_traits
//
// The MakeComplex_t alias declaration provides a convenient shortcut to access the nested \a Type
// of the MakeComplex class template. For instance, given the type \a T the following two type
// definitions are identical:

   \code
   using Type1 = typename blaze::MakeComplex<T>::Type;
   using Type2 = blaze::MakeComplex_t<T>;
   \endcode
*/
template< typename T >
using MakeComplex_t = typename MakeComplex<T>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
