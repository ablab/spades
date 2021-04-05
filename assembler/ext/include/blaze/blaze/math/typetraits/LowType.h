//=================================================================================================
/*!
//  \file blaze/math/typetraits/LowType.h
//  \brief Header file for the LowType type trait
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

#ifndef _BLAZE_MATH_TYPETRAITS_LOWTYPE_H_
#define _BLAZE_MATH_TYPETRAITS_LOWTYPE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Complex.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsSigned.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the LowType type trait.
// \ingroup math_traits
*/
template< typename T1        // First operand
        , typename T2        // Second operand
        , typename = void >  // Restricting condition
struct LowTypeHelper
{
 public:
   //**********************************************************************************************
   using Type = INVALID_TYPE;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for a small and a large integral type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsIntegral_v<T1> &&
                                   IsIntegral_v<T2> &&
                                   ( sizeof( T1 ) < sizeof( T2 ) ) > >
{
 public:
   //**********************************************************************************************
   using Type = T1;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for a large and a small integral type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsIntegral_v<T1> &&
                                   IsIntegral_v<T2> &&
                                   ( sizeof( T1 ) > sizeof( T2 ) ) > >
{
 public:
   //**********************************************************************************************
   using Type = T2;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for two integral types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsIntegral_v<T1> &&
                                   IsIntegral_v<T2> &&
                                   ( sizeof( T1 ) == sizeof( T2 ) ) > >
{
 public:
   //**********************************************************************************************
   using Type = If_t< IsSigned_v<T1>, T2, T1 >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for an integral and a floating point type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsIntegral_v<T1> && IsFloatingPoint_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = T1;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for a floating point and an integral type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsFloatingPoint_v<T1> && IsIntegral_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = T2;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for two floating point types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsFloatingPoint_v<T1> && IsFloatingPoint_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = If_t< ( sizeof( T1 ) < sizeof( T2 ) ), T1, T2 >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for a complex and another type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsComplex_v<T1> && !IsComplex_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = T2;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowTypeHelper class template for another and a complex type.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< !IsComplex_v<T1> && IsComplex_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = T1;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the LowType class template for two complex types.
// \ingroup math_traits
*/
template< typename T1, typename T2 >
struct LowTypeHelper< T1, T2
                     , EnableIf_t< IsComplex_v<T1> && IsComplex_v<T2> > >
{
 public:
   //**********************************************************************************************
   using Type = complex< typename LowTypeHelper< typename T1::value_type, typename T2::value_type >::Type >;
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Base template for the LowType type trait.
// \ingroup math_traits
//
// \section lowtype_general General
//
// The LowType class template determines the more significant, dominating data type of the two
// given data types \a T1 and \a T2. In case both \a T1 and \a T2 are built-in data types, the
// nested type \a Type is set to the smaller or unsigned data type. If case both built-in types
// have the same size and the same signedness, the selected type is implementation defined. In
// case no lower data type can be selected, \a Type is set to \a INVALID_TYPE.
//
// Per default, the LowType template provides support for the following built-in data types:
//
// <ul>
//    <li>Integral types</li>
//    <ul>
//       <li>unsigned char, signed char, char, wchar_t</li>
//       <li>char16_t, char32_t</li>
//       <li>unsigned short, short</li>
//       <li>unsigned int, int</li>
//       <li>unsigned long, long</li>
//       <li>std::size_t, std::ptrdiff_t (for certain 64-bit compilers)</li>
//    </ul>
//    <li>Floating point types</li>
//    <ul>
//       <li>float</li>
//       <li>double</li>
//       <li>long double</li>
//    </ul>
// </ul>
//
// Additionally, the Blaze library provides specializations for the following user-defined
// arithmetic types, wherever a less significant data type can be selected:
//
// <ul>
//    <li>std::complex</li>
//    <li>blaze::StaticVector</li>
//    <li>blaze::HybridVector</li>
//    <li>blaze::DynamicVector</li>
//    <li>blaze::CompressedVector</li>
//    <li>blaze::StaticMatrix</li>
//    <li>blaze::HybridMatrix</li>
//    <li>blaze::DynamicMatrix</li>
//    <li>blaze::CompressedMatrix</li>
//    <li>blaze::SymmetricMatrix</li>
//    <li>blaze::HermitianMatrix</li>
//    <li>blaze::LowerMatrix</li>
//    <li>blaze::UniLowerMatrix</li>
//    <li>blaze::StrictlyLowerMatrix</li>
//    <li>blaze::UpperMatrix</li>
//    <li>blaze::UniUpperMatrix</li>
//    <li>blaze::StrictlyUpperMatrix</li>
//    <li>blaze::DiagonalMatrix</li>
// </ul>
//
//
// \n \section lowtype_specializations Creating custom specializations
//
// It is possible to specialize the LowType template for additional user-defined data types.
// The following example shows the according specialization for two dynamic column vectors:

   \code
   template< typename T1, typename T2 >
   struct LowType< DynamicVector<T1,false>, DynamicVector<T2,false> >
   {
      using Type = DynamicVector< typename LowType<T1,T2>::Type, false >;
   };
   \endcode
*/
template< typename T1        // First operand
        , typename T2        // Second operand
        , typename = void >  // Restricting condition
struct LowType
{
 public:
   //**********************************************************************************************
   /*! \cond BLAZE_INTERNAL */
   using Type = typename LowTypeHelper<T1,T2>::Type;
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Auxiliary alias declaration for the LowType type trait.
// \ingroup math_type_traits
//
// The LowType_t alias declaration provides a convenient shortcut to access the nested \a Type of
// the LowType class template. For instance, given the types \a T1 and \a T2 the following two
// type definitions are identical:

   \code
   using Type1 = typename blaze::LowType<T1,T2>::Type;
   using Type2 = blaze::LowType_t<T1,T2>;
   \endcode
*/
template< typename T1    // First operand
        , typename T2 >  // Second operand
using LowType_t = typename LowType<T1,T2>::Type;
//*************************************************************************************************

} // namespace blaze

#endif
