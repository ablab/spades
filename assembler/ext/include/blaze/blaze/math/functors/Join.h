//=================================================================================================
/*!
//  \file blaze/math/functors/Join.h
//  \brief Header file for the Join functor
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

#ifndef _BLAZE_MATH_FUNCTORS_JOIN_H_
#define _BLAZE_MATH_FUNCTORS_JOIN_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/typetraits/IsPaddingEnabled.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/YieldsHermitian.h>
#include <blaze/math/typetraits/YieldsLower.h>
#include <blaze/math/typetraits/YieldsStrictlyLower.h>
#include <blaze/math/typetraits/YieldsStrictlyUpper.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/math/typetraits/YieldsUniLower.h>
#include <blaze/math/typetraits/YieldsUniUpper.h>
#include <blaze/math/typetraits/YieldsUpper.h>
#include <blaze/math/typetraits/YieldsZero.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for unpacking paired function arguments.
// \ingroup functors
*/
template< typename OP >  // Type of the operation
struct Join
{
 public:
   //**********************************************************************************************
   /*!\brief Constructor of the Join functor.
   //
   // \param op The wrapped operation.
   */
   explicit constexpr Join( const OP& op )
      : op_( op )  // The wrapped operation.
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given arguments.
   //
   // \param args The given arguments.
   // \return The result of the wrapped operation for the given arguments.
   */
   template< typename... Args >
   BLAZE_ALWAYS_INLINE constexpr decltype(auto) operator()( const Args&... args ) const
   {
      return op_( args... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given arguments.
   //
   // \param p The first two arguments in form of a pair.
   // \param args The remaining arguments.
   // \return The result of the wrapped operation for the given arguments.
   */
   template< typename T1, typename T2, typename... Args >
   BLAZE_ALWAYS_INLINE constexpr decltype(auto) operator()( const std::pair<T1,T2>& p, const Args&... args ) const
   {
      return (*this)( p.first, p.second, args... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a Args.
   //
   // \return \a true in case SIMD is enabled for the data types \a Args, \a false if not.
   */
   template< typename... Args >
   static constexpr bool simdEnabled() { return IsSIMDEnabled_v<OP,Args...>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return IsPaddingEnabled_v<OP>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given SIMD vectors.
   //
   // \param args The given SIMD vectors.
   // \return The result of the wrapped operation for the given SIMD vectors.
   */
   template< typename... Args >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const Args&... args ) const
   {
      return op_.load( args... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given SIMD vectors.
   //
   // \param p The first two SIMD vectors in form of a pair.
   // \param args The remaining SIMD vectors.
   // \return The result of the wrapped operation for the given SIMD vectors.
   */
   template< typename T1, typename T2, typename... Args >
   BLAZE_ALWAYS_INLINE decltype(auto) load( const std::pair<T1,T2>& p, const Args&... args ) const
   {
      return (*this).load( p.first, p.second, args... );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   OP op_;  //!< The wrapped operation.
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Wraps the given operation to unpack paired arguments.
// \ingroup functors
//
// \param op The operation to be wrapped.
// \return The wrapped operation.
//
// The \a join() function returns a function object that represents the given operation, but
// is able to unpack arguments in the form of pairs.
*/
template< typename OP >  // Type of the operation
constexpr Join<OP> join( const OP& op )
{
   return Join<OP>( op );
}
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename T, typename... Args >
struct YieldsUniform<Join<OP>,T,Args...>
   : public YieldsUniform<OP,T,Args...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSYMMETRIC SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsSymmetric<Join<OP>,MT,MTs...>
   : public YieldsSymmetric<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSHERMITIAN SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsHermitian<Join<OP>,MT,MTs...>
   : public YieldsHermitian<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsLower<Join<OP>,MT,MTs...>
   : public YieldsLower<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNILOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsUniLower<Join<OP>,MT,MTs...>
   : public YieldsUniLower<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYLOWER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsStrictlyLower<Join<OP>,MT,MTs...>
   : public YieldsStrictlyLower<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsUpper<Join<OP>,MT,MTs...>
   : public YieldsUpper<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsUniUpper<Join<OP>,MT,MTs...>
   : public YieldsUniUpper<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSSTRICTLYUPPER SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename MT, typename... MTs >
struct YieldsStrictlyUpper<Join<OP>,MT,MTs...>
   : public YieldsStrictlyUpper<OP,MT,MTs...>
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSZERO SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename T, typename... Args >
struct YieldsZero<Join<OP>,T,Args...>
   : public YieldsZero<OP,T,Args...>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
