//=================================================================================================
/*!
//  \file blaze/math/functors/Bind1st.h
//  \brief Header file for the Bind1st functor
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

#ifndef _BLAZE_MATH_FUNCTORS_BIND1ST_H_
#define _BLAZE_MATH_FUNCTORS_BIND1ST_H_


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
/*!\brief Generic wrapper for an operation with fixed 1st argument.
// \ingroup functors
*/
template< typename OP    // Type of the operation
        , typename A1 >  // Type of the bound argument
struct Bind1st
{
 public:
   //**********************************************************************************************
   /*!\brief Constructor of the Bind1st functor.
   //
   // \param op The wrapped operation.
   // \param a1 The 1st argument.
   */
   constexpr Bind1st( const OP& op, const A1& a1 )
      : op_( op )  // The wrapped operation.
      , a1_( a1 )  // The 1st argument
   {}
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns the result of the wrapped operation for the given arguments.
   //
   // \param args The remaining arguments.
   // \return The result of the wrapped operation for the given arguments.
   */
   template< typename... Args >
   BLAZE_ALWAYS_INLINE constexpr decltype(auto) operator()( Args&&... args ) const
   {
      return op_( a1_, std::forward<Args>( args )... );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether SIMD is enabled for the specified data types \a Args.
   //
   // \return \a true in case SIMD is enabled for the data types \a Args, \a false if not.
   */
   template< typename... Args >
   static constexpr bool simdEnabled() { return IsSIMDEnabled_v<OP,A1,Args...>; }
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
   // \param args The remaining SIMD vectors.
   // \return The result of the wrapped operation for the given SIMD vectors.
   */
   template< typename... Args >
   BLAZE_ALWAYS_INLINE decltype(auto) load( Args&&... args ) const
   {
      return op_.load( set( a1_ ), std::forward<Args>( args )... );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   OP op_;  //!< The wrapped operation.
   A1 a1_;  //!< The 1st argument.
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Binds the given object/value to the 1st parameter of the given operation.
// \ingroup functors
//
// \param op The operation to be wrapped.
// \param a1 The argument to be bound to the second parameter of the operation.
// \return The operation with bound 1st argument.
//
// The \a bind1st() function binds the given argument \a x to the 1st parameter of the given
// operation \a op.
*/
template< typename OP    // Type of the operation
        , typename A1 >  // Type of the bound argument
constexpr Bind1st<OP,A1> bind1st( const OP& op, const A1& a1 )
{
   return Bind1st<OP,A1>( op, a1 );
}
//*************************************************************************************************




//=================================================================================================
//
//  YIELDSUNIFORM SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename OP, typename A1, typename T >
struct YieldsUniform<Bind1st<OP,A1>,T>
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
template< typename OP, typename A1, typename MT >
struct YieldsSymmetric<Bind1st<OP,A1>,MT>
   : public YieldsSymmetric<OP,MT>
{};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
