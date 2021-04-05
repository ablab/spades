//=================================================================================================
/*!
//  \file blaze/math/smp/openmp/DenseMatrix.h
//  \brief Header file for the OpenMP-based dense matrix SMP implementation
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

#ifndef _BLAZE_MATH_SMP_OPENMP_DENSEMATRIX_H_
#define _BLAZE_MATH_SMP_OPENMP_DENSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <omp.h>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/SMPAssignable.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/simd/SIMDTrait.h>
#include <blaze/math/smp/ParallelSection.h>
#include <blaze/math/smp/SerialSection.h>
#include <blaze/math/smp/ThreadMapping.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/views/Submatrix.h>
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
//  OPENMP-BASED ASSIGNMENT KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the OpenMP-based SMP (compound) assignment of a dense matrix to a dense matrix.
// \ingroup math
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side dense matrix to be assigned.
// \param op The (compound) assignment operation.
// \return void
//
// This function is the backend implementation of the OpenMP-based SMP assignment of a dense
// matrix to a dense matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1   // Type of the left-hand side dense matrix
        , bool SO1       // Storage order of the left-hand side dense matrix
        , typename MT2   // Type of the right-hand side dense matrix
        , bool SO2       // Storage order of the right-hand side dense matrix
        , typename OP >  // Type of the assignment operation
void openmpAssign( DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;

   constexpr bool simdEnabled( MT1::simdEnabled && MT2::simdEnabled && IsSIMDCombinable_v<ET1,ET2> );
   constexpr size_t SIMDSIZE( SIMDTrait< ElementType_t<MT1> >::size );

   const bool lhsAligned( (*lhs).isAligned() );
   const bool rhsAligned( (*rhs).isAligned() );

   const int threads( omp_get_num_threads() );
   const ThreadMapping threadmap( createThreadMapping( threads, *rhs ) );

   const size_t addon1     ( ( ( (*rhs).rows() % threadmap.first ) != 0UL )? 1UL : 0UL );
   const size_t equalShare1( (*rhs).rows() / threadmap.first + addon1 );
   const size_t rest1      ( equalShare1 & ( SIMDSIZE - 1UL ) );
   const size_t rowsPerThread( ( simdEnabled && rest1 )?( equalShare1 - rest1 + SIMDSIZE ):( equalShare1 ) );

   const size_t addon2     ( ( ( (*rhs).columns() % threadmap.second ) != 0UL )? 1UL : 0UL );
   const size_t equalShare2( (*rhs).columns() / threadmap.second + addon2 );
   const size_t rest2      ( equalShare2 & ( SIMDSIZE - 1UL ) );
   const size_t colsPerThread( ( simdEnabled && rest2 )?( equalShare2 - rest2 + SIMDSIZE ):( equalShare2 ) );

#pragma omp for schedule(dynamic,1) nowait
   for( int i=0; i<threads; ++i )
   {
      const size_t row   ( ( i / threadmap.second ) * rowsPerThread );
      const size_t column( ( i % threadmap.second ) * colsPerThread );

      if( row >= (*rhs).rows() || column >= (*rhs).columns() )
         continue;

      const size_t m( min( rowsPerThread, (*rhs).rows()    - row    ) );
      const size_t n( min( colsPerThread, (*rhs).columns() - column ) );

      if( simdEnabled && lhsAligned && rhsAligned ) {
         auto       target( submatrix<aligned>( *lhs, row, column, m, n ) );
         const auto source( submatrix<aligned>( *rhs, row, column, m, n ) );
         op( target, source );
      }
      else if( simdEnabled && lhsAligned ) {
         auto       target( submatrix<aligned>( *lhs, row, column, m, n ) );
         const auto source( submatrix<unaligned>( *rhs, row, column, m, n ) );
         op( target, source );
      }
      else if( simdEnabled && rhsAligned ) {
         auto       target( submatrix<unaligned>( *lhs, row, column, m, n ) );
         const auto source( submatrix<aligned>( *rhs, row, column, m, n ) );
         op( target, source );
      }
      else {
         auto       target( submatrix<unaligned>( *lhs, row, column, m, n ) );
         const auto source( submatrix<unaligned>( *rhs, row, column, m, n ) );
         op( target, source );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Backend of the OpenMP-based SMP (compound) assignment of a sparse matrix to a dense matrix.
// \ingroup math
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side sparse matrix to be assigned.
// \param op The (compound) assignment operation.
// \return void
//
// This function is the backend implementation of the OpenMP-based SMP assignment of a sparse
// matrix to a dense matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1   // Type of the left-hand side dense matrix
        , bool SO1       // Storage order of the left-hand side dense matrix
        , typename MT2   // Type of the right-hand side sparse matrix
        , bool SO2       // Storage order of the right-hand side sparse matrix
        , typename OP >  // Type of the assignment operation
void openmpAssign( DenseMatrix<MT1,SO1>& lhs, const SparseMatrix<MT2,SO2>& rhs, OP op )
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( isParallelSectionActive(), "Invalid call outside a parallel section" );

   const size_t threads( omp_get_num_threads() );
   const ThreadMapping threadmap( createThreadMapping( threads, *rhs ) );

   const size_t addon1       ( ( ( (*rhs).rows() % threadmap.first ) != 0UL )? 1UL : 0UL );
   const size_t rowsPerThread( (*rhs).rows() / threadmap.first + addon1 );

   const size_t addon2       ( ( ( (*rhs).columns() % threadmap.second ) != 0UL )? 1UL : 0UL );
   const size_t colsPerThread( (*rhs).columns() / threadmap.second + addon2 );

#pragma omp for schedule(dynamic,1) nowait
   for( size_t i=0; i<threads; ++i )
   {
      const size_t row   ( ( i / threadmap.second ) * rowsPerThread );
      const size_t column( ( i % threadmap.second ) * colsPerThread );

      if( row >= (*rhs).rows() || column >= (*rhs).columns() )
         continue;

      const size_t m( min( rowsPerThread, (*lhs).rows()    - row    ) );
      const size_t n( min( colsPerThread, (*lhs).columns() - column ) );

      auto       target( submatrix<unaligned>( *lhs, row, column, m, n ) );
      const auto source( submatrix<unaligned>( *rhs, row, column, m, n ) );
      op( target, source );
   }
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
/*!\brief Default implementation of the OpenMP-based SMP assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be assigned.
// \return void
//
// This function implements the default OpenMP-based SMP assignment to a dense matrix. Due to
// the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands are
// not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && ( !IsSMPAssignable_v<MT1> || !IsSMPAssignable_v<MT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   assign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the OpenMP-based SMP assignment to a dense matrix.
// \ingroup math
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be assigned.
// \return void
//
// This function implements the OpenMP-based SMP assignment to a dense matrix. Due to the
// explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && IsSMPAssignable_v<MT1> && IsSMPAssignable_v<MT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         assign( *lhs, *rhs );
      }
      else {
#pragma omp parallel shared( lhs, rhs )
         openmpAssign( *lhs, *rhs, []( auto& a, const auto& b ){ assign( a, b ); } );
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
/*!\brief Default implementation of the OpenMP-based SMP addition assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function implements the default OpenMP-based SMP addition assignment to a dense matrix.
// Due to the explicit application of the SFINAE principle, this function can only be selected
// by the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && ( !IsSMPAssignable_v<MT1> || !IsSMPAssignable_v<MT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   addAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the OpenMP-based SMP addition assignment to a dense matrix.
// \ingroup math
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be added.
// \return void
//
// This function implements the OpenMP-based SMP addition assignment to a dense matrix. Due to
// the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands are
// not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpAddAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && IsSMPAssignable_v<MT1> && IsSMPAssignable_v<MT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         addAssign( *lhs, *rhs );
      }
      else {
#pragma omp parallel shared( lhs, rhs )
         openmpAssign( *lhs, *rhs, []( auto& a, const auto& b ){ addAssign( a, b ); } );
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
/*!\brief Default implementation of the OpenMP-based SMP subtracction assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function implements the default OpenMP-based SMP subtraction assignment to a dense matrix.
// Due to the explicit application of the SFINAE principle, this function can only be selected by
// the compiler in case both operands are SMP-assignable and the element types of both operands
// are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && ( !IsSMPAssignable_v<MT1> || !IsSMPAssignable_v<MT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   subAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the OpenMP-based SMP subtracction assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be subtracted.
// \return void
//
// This function implements the default OpenMP-based SMP subtraction assignment of a matrix to a
// dense matrix. Due to the explicit application of the SFINAE principle, this function can only
// be selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSubAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && IsSMPAssignable_v<MT1> && IsSMPAssignable_v<MT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         subAssign( *lhs, *rhs );
      }
      else {
#pragma omp parallel shared( lhs, rhs )
         openmpAssign( *lhs, *rhs, []( auto& a, const auto& b ){ subAssign( a, b ); } );
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SCHUR PRODUCT ASSIGNMENT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the OpenMP-based SMP Schur product assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix for the Schur product.
// \return void
//
// This function implements the default OpenMP-based SMP Schur product assignment to a dense
// matrix. Due to the explicit application of the SFINAE principle, this function can only be
// selected by the compiler in case both operands are SMP-assignable and the element types of
// both operands are not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && ( !IsSMPAssignable_v<MT1> || !IsSMPAssignable_v<MT2> ) >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   schurAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Implementation of the OpenMP-based SMP Schur product assignment to a dense matrix.
// \ingroup math
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix for the Schur product.
// \return void
//
// This function implements the OpenMP-based SMP Schur product assignment to a dense matrix. Due
// to the explicit application of the SFINAE principle, this function can only be selected by the
// compiler in case both operands are SMP-assignable and the element types of both operands are
// not SMP-assignable.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpSchurAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> && IsSMPAssignable_v<MT1> && IsSMPAssignable_v<MT2> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT1> );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SMP_ASSIGNABLE( ElementType_t<MT2> );

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   BLAZE_PARALLEL_SECTION
   {
      if( isSerialSectionActive() || !(*rhs).canSMPAssign() ) {
         schurAssign( *lhs, *rhs );
      }
      else {
#pragma omp parallel shared( lhs, rhs )
         openmpAssign( *lhs, *rhs, []( auto& a, const auto& b ){ schurAssign( a, b ); } );
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
/*!\brief Default implementation of the OpenMP-based SMP multiplication assignment to a dense matrix.
// \ingroup smp
//
// \param lhs The target left-hand side dense matrix.
// \param rhs The right-hand side matrix to be multiplied.
// \return void
//
// This function implements the default OpenMP-based SMP multiplication assignment to a dense
// matrix.\n
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side matrix
        , typename MT2  // Type of the right-hand side matrix
        , bool SO2 >    // Storage order of the right-hand side matrix
inline auto smpMultAssign( Matrix<MT1,SO1>& lhs, const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< IsDenseMatrix_v<MT1> >
{
   BLAZE_FUNCTION_TRACE;

   BLAZE_INTERNAL_ASSERT( (*lhs).rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( (*lhs).columns() == (*rhs).columns(), "Invalid number of columns" );

   multAssign( *lhs, *rhs );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COMPILE TIME CONSTRAINT
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
namespace {

BLAZE_STATIC_ASSERT( BLAZE_OPENMP_PARALLEL_MODE );

}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
