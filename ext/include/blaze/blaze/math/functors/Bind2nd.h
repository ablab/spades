//=================================================================================================
/*!
//  \file blaze/math/functors/Bind2nd.h
//  \brief Header file for the Bind2nd functor
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

#ifndef _BLAZE_MATH_FUNCTORS_BIND2ND_H_
#define _BLAZE_MATH_FUNCTORS_BIND2ND_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/simd/Set.h>
#include <blaze/math/typetraits/IsSIMDEnabled.h>
#include <blaze/math/typetraits/YieldsSymmetric.h>
#include <blaze/math/typetraits/YieldsUniform.h>
#include <blaze/system/Inline.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Generic wrapper for an operation with fixed 2nd argument.
// \ingroup functors
*/
template< typename OP    // Type of the operation
        , typename A2 >  // Type of the bound argument
struct Bind2nd
{
 public:
   //**********************************************************************************************
   /*!\brief Constructor of the Bind2nd functor.
   //
   // \param op The wrapped operation.
   // \param a2 The 2nd argument.
   */
   constexpr Bind2nd( const OP& op, const A2& a2 )
      : op_( op )  // The wrapped operation.
      , a2_( a2 )  // The 2nd argument
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given arguments.
   //
   // \param a1 The 1st argument.
   // \param args The remaining arguments.
   // \return The result of the wrapped operation for the given arguments.
   */
   template< typename A1, typename... Args >
   BLAZE_ALWAYS_INLINE constexpr decltype(auto) operator()( A1&& a1, Args&&... args ) const
   {
      return op_( std::forward<A1>( a1 ), a2_, std::forward<Args>( args )... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a Ts.
   //
   // \return \a true in case SIMD is enabled for the data types \a Ts, \a false if not.
   */
   template< typename A1, typename... Args >
   static constexpr bool simdEnabled() { return IsSIMDEnabled_v<OP,A1,A2,Args...>; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the operation supports padding, i.e. whether it can deal with zeros.
   //
   // \return \a true in case padding is supported, \a false if not.
   */
   static constexpr bool paddingEnabled() { return false; }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given SIMD vectors.
   //
   // \param a1 The 1st given SIMD vector.
   // \param args The remaining SIMD vectors.
   // \return The result of the wrapped operation for the given SIMD vectors.
   */
   template< typename A1, typename... Args >
   BLAZE_ALWAYS_INLINE decltype(auto) load( A1&& a1, Args&&... args ) const
   {
      return op_.load( std::forward<A1>( a1 ), set( a2_ ), std::forward<Args>( args )... );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   OP op_;  //!< The wrapped operation.
   A2 a2_;  //!< The 2nd argument.
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Binds the given object/value to the 2nd parameter of the given operation.
// \ingroup functors
//
// \param op The operation to be wrapped.
// \param a2 The argument to be bound to the second parameter of the operation.
// \return The operation with bound 2nd argument.
//
// The \a bind2nd() function binds the given argument \a x to the 2nd parameter of the given
// operation \a op.
*/
template< typename OP    // Type of the operation
        , typename A2 >  // Type of the bound argument
constexpr Bind2nd<OP,A2> bind2nd( const OP& op, const A2& a2 )
{
   return Bind2nd<OP,A2>( op, a2 );
}
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename A2, typename T >
struct YieldsUniform<Bind2nd<OP,A2>,T>
   : public YieldsUniform<OP,T>
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
template< typename OP, typename A2, typename MT >
struct YieldsSymmetric<Bind2nd<OP,A2>,MT>
   : public YieldsSymmetric<OP,MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
