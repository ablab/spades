//=================================================================================================
/*!
//  \file blaze/math/smp/threads/DenseVector.h
//  \brief Header file for the C++11/Boost thread-based dense vector SMP implementation
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

#ifndef _BLAZE_MATH_SMP_THREADS_DENSEVECTOR_H_
#define _BLAZE_MATH_SMP_THREADS_DENSEVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/SMPAssignable.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/smp/ParallelSection.h>
#include <blaze/math/smp/SerialSection.h>
#include <blaze/math/smp/threads/ThreadBackend.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/views/Subvector.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/SMP.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/FunctionTrace.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  THREAD-BASED ASSIGNMENT KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP (compound) assignment of a dense vector to
//        a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side dense vector to be assigned.
// \param op The (compound) assignment operation.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP assignment
// of a dense vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1   // Type of the left-hand side dense vector
        , bool TF1       // Transpose flag of the left-hand side dense vector
        , typename VT2   // Type of the right-hand side dense vector
        , bool TF2       // Transpose flag of the right-hand side dense vector
        , typename OP >  // Type of the assignment operation
void threadAssign( DenseVector<VT1,TF1>& lhs, const DenseVector<VT2,TF2>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   using ET1 = ElementType_t<VT1>;
   using ET2 = ElementType_t<VT2>;

   constexpr bool simdEnabled( VT1::simdEnabled && VT2::simdEnabled && IsSIMDCombinable_v<ET1,ET2> );
   constexpr size_t SIMDSIZE( SIMDTrait< ElementType_t<VT1> >::size );

   const bool lhsAligned( (*lhs).isAligned() );
   const bool rhsAligned( (*rhs).isAligned() );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (*lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t equalShare   ( (*lhs).size() / threads + addon );
   const size_t rest         ( equalShare & ( SIMDSIZE - 1UL ) );
   const size_t sizePerThread( ( simdEnabled && rest )?( equalShare - rest + SIMDSIZE ):( equalShare ) );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (*lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (*lhs).size() - index ) );

      if( simdEnabled && lhsAligned && rhsAligned ) {
         auto       target( subvector<aligned>( *lhs, index, size, unchecked ) );
         const auto source( subvector<aligned>( *rhs, index, size, unchecked ) );
         TheThreadBackend::schedule( target, source, op );
      }
      else if( simdEnabled && lhsAligned ) {
         auto       target( subvector<aligned>( *lhs, index, size, unchecked ) );
         const auto source( subvector<unaligned>( *rhs, index, size, unchecked ) );
         TheThreadBackend::schedule( target, source, op );
      }
      else if( simdEnabled && rhsAligned ) {
         auto       target( subvector<unaligned>( *lhs, index, size, unchecked ) );
         const auto source( subvector<aligned>( *rhs, index, size, unchecked ) );
         TheThreadBackend::schedule( target, source, op );
      }
      else {
         auto       target( subvector<unaligned>( *lhs, index, size, unchecked ) );
         const auto source( subvector<unaligned>( *rhs, index, size, unchecked ) );
         TheThreadBackend::schedule( target, source, op );
      }
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the C++11/Boost thread-based SMP assignment of a sparse vector to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be assigned.
// \param op The (compound) assignment operation.
// \return void
//
// This function is the backend implementation of the C++11/Boost thread-based SMP assignment
// of a sparse vector to a dense vector.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1   // Type of the left-hand side dense vector
        , bool TF1       // Transpose flag of the left-hand side dense vector
        , typename VT2   // Type of the right-hand side sparse vector
        , bool TF2       // Transpose flag of the right-hand side sparse vector
        , typename OP >  // Type of the assignment operation
void threadAssign( DenseVector<VT1,TF1>& lhs, const SparseVector<VT2,TF2>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   const size_t threads      ( TheThreadBackend::size() );
   const size_t addon        ( ( ( (*lhs).size() % threads ) != 0UL )? 1UL : 0UL );
   const size_t sizePerThread( (*lhs).size() / threads + addon );

   for( size_t i=0UL; i<threads; ++i )
   {
      const size_t index( i*sizePerThread );

      if( index >= (*lhs).size() )
         continue;

      const size_t size( min( sizePerThread, (*lhs).size() - index ) );
      auto       target( subvector<unaligned>( *lhs, index, size, unchecked ) );
      const auto source( subvector<unaligned>( *rhs, index, size, unchecked ) );
      TheThreadBackend::schedule( target, source, op );
   }

   TheThreadBackend::wait();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  PLAIN ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP assignment to a dense vector.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && ( !IsSMPAssignable_v<VT1> || !IsSMPAssignable_v<VT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   assign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function performs the C++11/Boost thread-based SMP assignment to a dense vector. Due to
// the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands are
// not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && IsSMPAssignable_v<VT1> && IsSMPAssignable_v<VT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         assign( *lhs, *rhs );
      }
      else {
         threadAssign( *lhs, *rhs, []( auto& a, const auto& b ){ assign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDITION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP addition assignment to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP addition assignment to a
// dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && ( !IsSMPAssignable_v<VT1> || !IsSMPAssignable_v<VT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   addAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP addition assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function performs the C++11/Boost thread-based SMP addition assignment to a dense vector.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpAddAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && IsSMPAssignable_v<VT1> && IsSMPAssignable_v<VT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         addAssign( *lhs, *rhs );
      }
      else {
         threadAssign( *lhs, *rhs, []( auto& a, const auto& b ){ addAssign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRACTION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP subtraction assignment to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP subtraction assignment to a
// dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && ( !IsSMPAssignable_v<VT1> || !IsSMPAssignable_v<VT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   subAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP subtraction assignment to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function performs the C++11/Boost thread-based SMP subtraction assignment to a dense
// vector. Due to the explicit application of the SFINAE principle, this function can only be
// selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpSubAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && IsSMPAssignable_v<VT1> && IsSMPAssignable_v<VT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         subAssign( *lhs, *rhs );
      }
      else {
         threadAssign( *lhs, *rhs, []( auto& a, const auto& b ){ subAssign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTIPLICATION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP multiplication assignment
//        to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP multiplication assignment to
// a dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && ( !IsSMPAssignable_v<VT1> || !IsSMPAssignable_v<VT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   multAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP multiplication assignment to a dense
//        vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function performs the C++11/Boost thread-based SMP multiplication assignment to a dense
// vector. Due to the explicit application of the SFINAE principle, this function can only be
// selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpMultAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && IsSMPAssignable_v<VT1> && IsSMPAssignable_v<VT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         multAssign( *lhs, *rhs );
      }
      else {
         threadAssign( *lhs, *rhs, []( auto& a, const auto& b ){ multAssign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVISION ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the C++11/Boost thread-based SMP division assignment to a
//        dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side vector to be assigned.
// \return void
//
// This function implements the default C++11/Boost thread-based SMP division assignment to a
// dense vector. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpDivAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && ( !IsSMPAssignable_v<VT1> || !IsSMPAssignable_v<VT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   divAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the C++11/Boost thread-based SMP division assignment to a dense vector.
// \ingroup smp
//
// \param lhs The target left-hand side dense vector.
// \param rhs The right-hand side sparse vector divisor.
// \return void
//
// This function performs the C++11/Boost thread-based SMP division assignment to a dense vector.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT1  // Type of the left-hand side dense vector
        , bool TF1      // Transpose flag of the left-hand side dense vector
        , typename VT2  // Type of the right-hand side vector
        , bool TF2 >    // Transpose flag of the right-hand side vector
inline auto smpDivAssign( Vector<VT1,TF1>& lhs, const Vector<VT2,TF2>& rhs )
   -> EnableIf_t< IsDenseVector_v<VT1> && IsSMPAssignable_v<VT1> && IsSMPAssignable_v<VT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<VT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).size() == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         divAssign( *lhs, *rhs );
      }
      else {
         threadAssign( *lhs, *rhs, []( auto& a, const auto& b ){ divAssign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( BLAZE_CPP_THREADS_PARALLEL_MODE || BLAZE_BOOST_THREADS_PARALLEL_MODE );

}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
