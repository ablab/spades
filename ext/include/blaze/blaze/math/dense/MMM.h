//=================================================================================================
/*!
//  \file blaze/math/dense/MMM.h
//  \brief Header file for the dense matrix multiplication kernels
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

#ifndef _BLAZE_MATH_DENSE_MMM_H_
#define _BLAZE_MATH_DENSE_MMM_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/SIMDCombinable.h>
#include <blaze/math/constraints/StrictlyLower.h>
#include <blaze/math/constraints/StrictlyUpper.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/UniLower.h>
#include <blaze/math/constraints/UniUpper.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/dense/DynamicMatrix.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/Submatrix.h>
#include <blaze/system/Blocking.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/StaticAssert.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blaze {

//=================================================================================================
//
//  GENERAL DENSE MATRIX MULTIPLICATION KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a general dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side row-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a general dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// row-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void mmm( DenseMatrix<MT1,false>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t JBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( JBLOCK >= SIMDSIZE && JBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( M, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, JBLOCK );

   if( isDefault( beta ) ) {
      reset( *C );
   }
   else if( !isOne( beta ) ) {
      (*C) *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t iend  ( IsUpper_v<MT2> ? kk+kblock : M );
      const size_t isize ( iend - ibegin );

      A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ibegin, kk, isize, kblock, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( ( IsLower_v<MT3> && kk+kblock <= jj ) ||
             ( IsUpper_v<MT3> && jj+jblock <= kk ) ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jj, kblock, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               size_t j( 0UL );

               for( ; (j+2UL) <= jblock; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a2 * b1;
                     xmm4  += a2 * b2;
                     xmm5  += a3 * b1;
                     xmm6  += a3 * b2;
                     xmm7  += a4 * b1;
                     xmm8  += a4 * b2;
                     xmm9  += a5 * b1;
                     xmm10 += a5 * b2;
                  }

                  (*C)(ibegin+i    ,jj+j    ) += sum( xmm1  ) * alpha;
                  (*C)(ibegin+i    ,jj+j+1UL) += sum( xmm2  ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j    ) += sum( xmm3  ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j+1UL) += sum( xmm4  ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j    ) += sum( xmm5  ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j+1UL) += sum( xmm6  ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j    ) += sum( xmm7  ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j+1UL) += sum( xmm8  ) * alpha;
                  (*C)(ibegin+i+4UL,jj+j    ) += sum( xmm9  ) * alpha;
                  (*C)(ibegin+i+4UL,jj+j+1UL) += sum( xmm10 ) * alpha;
               }

               if( j<jblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                     xmm5 += a5 * b1;
                  }

                  (*C)(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
                  (*C)(ibegin+i+4UL,jj+j) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               size_t j( 0UL );

               for( ; (j+2UL) <= jblock; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a2 * b1;
                     xmm4 += a2 * b2;
                     xmm5 += a3 * b1;
                     xmm6 += a3 * b2;
                     xmm7 += a4 * b1;
                     xmm8 += a4 * b2;
                  }

                  (*C)(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
                  (*C)(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j    ) += sum( xmm5 ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j+1UL) += sum( xmm6 ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j    ) += sum( xmm7 ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j+1UL) += sum( xmm8 ) * alpha;
               }

               if( j<jblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                  }

                  (*C)(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  (*C)(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  (*C)(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            size_t j( 0UL );

            for( ; (j+4UL) <= jblock; j+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );
                  const SIMDType b3( B2.load(k,j+2UL) );
                  const SIMDType b4( B2.load(k,j+3UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a1 * b3;
                  xmm4 += a1 * b4;
                  xmm5 += a2 * b1;
                  xmm6 += a2 * b2;
                  xmm7 += a2 * b3;
                  xmm8 += a2 * b4;
               }

               (*C)(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               (*C)(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               (*C)(ibegin+i    ,jj+j+2UL) += sum( xmm3 ) * alpha;
               (*C)(ibegin+i    ,jj+j+3UL) += sum( xmm4 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j    ) += sum( xmm5 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j+1UL) += sum( xmm6 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j+2UL) += sum( xmm7 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j+3UL) += sum( xmm8 ) * alpha;
            }

            for( ; (j+2UL) <= jblock; j+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               (*C)(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               (*C)(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( j<jblock )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
               }

               (*C)(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
               (*C)(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
            }
         }

         if( i<isize )
         {
            size_t j( 0UL );

            for( ; (j+2UL) <= jblock; j+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               (*C)(ibegin+i,jj+j    ) += sum( xmm1 ) * alpha;
               (*C)(ibegin+i,jj+j+1UL) += sum( xmm2 ) * alpha;
            }

            if( j<jblock )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j);
               }

               (*C)(ibegin+i,jj+j) += sum( xmm1 ) * alpha;
            }
         }

         jj += jblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t isize ( M - ibegin );

      A2 = serial( submatrix( A, ibegin, kk, isize, ksize, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( IsUpper_v<MT3> && jj+jblock <= kk ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix( B, kk, jj, ksize, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               size_t j( 0UL );

               for( ; (j+2UL) <= jblock; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+4UL,jj+j    ) += A2(i+4UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+4UL,jj+j+1UL) += A2(i+4UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+4UL,jj+j) += A2(i+4UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               size_t j( 0UL );

               for( ; (j+2UL) <= jblock; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     (*C)(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     (*C)(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            size_t j( 0UL );

            for( ; (j+2UL) <= jblock; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  (*C)(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  (*C)(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  (*C)(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                  (*C)(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }
         }

         if( i<isize )
         {
            size_t j( 0UL );

            for( ; (j+2UL) <= jblock; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ibegin+i,jj+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  (*C)(ibegin+i,jj+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ibegin+i,jj+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         jj += jblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a general dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a general dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// column-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void mmm( DenseMatrix<MT1,true>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE        ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t IBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( IBLOCK >= SIMDSIZE && IBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( IBLOCK, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, N );

   if( isDefault( beta ) ) {
      reset( *C );
   }
   else if( !isOne( beta ) ) {
      (*C) *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jend  ( IsLower_v<MT3> ? kk+kblock : N );
      const size_t jsize ( jend - jbegin );

      B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jbegin, kblock, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( ( IsLower_v<MT2> && ii+iblock <= kk ) ||
             ( IsUpper_v<MT2> && kk+kblock <= ii ) ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ii, kk, iblock, kblock, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET3> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               size_t i( 0UL );

               for( ; (i+2UL) <= iblock; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );
                     const SIMDType b5( B2.load(k,j+4UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a1 * b3;
                     xmm4  += a1 * b4;
                     xmm5  += a1 * b5;
                     xmm6  += a2 * b1;
                     xmm7  += a2 * b2;
                     xmm8  += a2 * b3;
                     xmm9  += a2 * b4;
                     xmm10 += a2 * b5;
                  }

                  (*C)(ii+i    ,jbegin+j    ) += sum( xmm1  ) * alpha;
                  (*C)(ii+i    ,jbegin+j+1UL) += sum( xmm2  ) * alpha;
                  (*C)(ii+i    ,jbegin+j+2UL) += sum( xmm3  ) * alpha;
                  (*C)(ii+i    ,jbegin+j+3UL) += sum( xmm4  ) * alpha;
                  (*C)(ii+i    ,jbegin+j+4UL) += sum( xmm5  ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j    ) += sum( xmm6  ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+1UL) += sum( xmm7  ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+2UL) += sum( xmm8  ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+3UL) += sum( xmm9  ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+4UL) += sum( xmm10 ) * alpha;
               }

               if( i<iblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                     xmm5 += a1 * B2.load(k,j+4UL);
                  }

                  (*C)(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  (*C)(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  (*C)(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  (*C)(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  (*C)(ii+i,jbegin+j+4UL) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               size_t i( 0UL );

               for( ; (i+2UL) <= iblock; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a1 * b3;
                     xmm4 += a1 * b4;
                     xmm5 += a2 * b1;
                     xmm6 += a2 * b2;
                     xmm7 += a2 * b3;
                     xmm8 += a2 * b4;
                  }

                  (*C)(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
                  (*C)(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  (*C)(ii+i    ,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  (*C)(ii+i    ,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j    ) += sum( xmm5 ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+2UL) += sum( xmm7 ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+3UL) += sum( xmm8 ) * alpha;
               }

               if( i<iblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                  }

                  (*C)(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  (*C)(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  (*C)(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  (*C)(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            size_t i( 0UL );

            for( ; (i+4UL) <= iblock; i+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );
                  const SIMDType a3( A2.load(i+2UL,k) );
                  const SIMDType a4( A2.load(i+3UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
                  xmm5 += a3 * b1;
                  xmm6 += a3 * b2;
                  xmm7 += a4 * b1;
                  xmm8 += a4 * b2;
               }

               (*C)(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               (*C)(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               (*C)(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               (*C)(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
               (*C)(ii+i+2UL,jbegin+j    ) += sum( xmm5 ) * alpha;
               (*C)(ii+i+2UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
               (*C)(ii+i+3UL,jbegin+j    ) += sum( xmm7 ) * alpha;
               (*C)(ii+i+3UL,jbegin+j+1UL) += sum( xmm8 ) * alpha;
            }

            for( ; (i+2UL) <= iblock; i+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               (*C)(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               (*C)(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               (*C)(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               (*C)(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( i<iblock )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               (*C)(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
               (*C)(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
            }
         }

         if( j<jsize )
         {
            size_t i( 0UL );

            for( ; (i+2UL) <= iblock; i+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += A2.load(i    ,k) * b1;
                  xmm2 += A2.load(i+1UL,k) * b1;
               }

               (*C)(ii+i    ,jbegin+j) += sum( xmm1 ) * alpha;
               (*C)(ii+i+1UL,jbegin+j) += sum( xmm2 ) * alpha;
            }

            if( i<iblock )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  xmm1 += A2.load(i,k) * B2.load(k,j);
               }

               (*C)(ii+i,jbegin+j) += sum( xmm1 ) * alpha;
            }
         }

         ii += iblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jsize ( N - jbegin );

      B2 = serial( submatrix( B, kk, jbegin, ksize, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( IsLower_v<MT2> && ii+iblock <= kk ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix( A, ii, kk, iblock, ksize, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               size_t i( 0UL );

               for( ; (i+2UL) <= iblock; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     (*C)(ii+i    ,jbegin+j+4UL) += A2(i    ,k) * B2(k,j+4UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+4UL) += A2(i+1UL,k) * B2(k,j+4UL) * alpha;
                  }
               }

               if( i<iblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                     (*C)(ii+i,jbegin+j+4UL) += A2(i,k) * B2(k,j+4UL) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               size_t i( 0UL );

               for( ; (i+2UL) <= iblock; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                  }
               }

               if( i<iblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     (*C)(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     (*C)(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     (*C)(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     (*C)(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                  }
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            size_t i( 0UL );

            for( ; (i+2UL) <= iblock; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  (*C)(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  (*C)(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  (*C)(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( i<iblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  (*C)(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }
         }

         if( j<jsize )
         {
            size_t i( 0UL );

            for( ; (i+2UL) <= iblock; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ii+i    ,jbegin+j) += A2(i    ,k) * B2(k,j) * alpha;
                  (*C)(ii+i+1UL,jbegin+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }

            if( i<iblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  (*C)(ii+i,jbegin+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         ii += iblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a general dense matrix/dense matrix multiplication (\f$ C=A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \return void
//
// This function implements the compute kernel for a general dense matrix/dense matrix
// multiplication of the form \f$ C=A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must
// provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3 >
inline void mmm( MT1& C, const MT2& A, const MT3& B )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   mmm( C, A, B, ET1(1), ET1(0) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOWER DENSE MATRIX MULTIPLICATION KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a lower dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side row-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a lower dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// row-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void lmmm( DenseMatrix<MT1,false>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE             ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE          ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t JBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( JBLOCK >= SIMDSIZE && JBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( M, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, JBLOCK );

   decltype(auto) c( derestrict( *C ) );

   if( isDefault( beta ) ) {
      reset( c );
   }
   else if( !isOne( beta ) ) {
      c *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t iend  ( IsUpper_v<MT2> ? kk+kblock : M );
      const size_t isize ( iend - ibegin );

      A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ibegin, kk, isize, kblock, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( ( IsLower_v<MT3> && kk+kblock <= jj ) ||
             ( IsUpper_v<MT3> && jj+jblock <= kk ) ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jj, kblock, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               if( jj > ibegin+i+4UL ) continue;

               const size_t jend( min( ibegin+i-jj+5UL, jblock ) );
               size_t j( 0UL );

               for( ; (j+2UL) <= jend; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a2 * b1;
                     xmm4  += a2 * b2;
                     xmm5  += a3 * b1;
                     xmm6  += a3 * b2;
                     xmm7  += a4 * b1;
                     xmm8  += a4 * b2;
                     xmm9  += a5 * b1;
                     xmm10 += a5 * b2;
                  }

                  c(ibegin+i    ,jj+j    ) += sum( xmm1  ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += sum( xmm2  ) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += sum( xmm3  ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4  ) * alpha;
                  c(ibegin+i+2UL,jj+j    ) += sum( xmm5  ) * alpha;
                  c(ibegin+i+2UL,jj+j+1UL) += sum( xmm6  ) * alpha;
                  c(ibegin+i+3UL,jj+j    ) += sum( xmm7  ) * alpha;
                  c(ibegin+i+3UL,jj+j+1UL) += sum( xmm8  ) * alpha;
                  c(ibegin+i+4UL,jj+j    ) += sum( xmm9  ) * alpha;
                  c(ibegin+i+4UL,jj+j+1UL) += sum( xmm10 ) * alpha;
               }

               if( j<jend )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                     xmm5 += a5 * b1;
                  }

                  c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  c(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  c(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
                  c(ibegin+i+4UL,jj+j) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               if( jj > ibegin+i+3UL ) continue;

               const size_t jend( min( ibegin+i-jj+4UL, jblock ) );
               size_t j( 0UL );

               for( ; (j+2UL) <= jend; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a2 * b1;
                     xmm4 += a2 * b2;
                     xmm5 += a3 * b1;
                     xmm6 += a3 * b2;
                     xmm7 += a4 * b1;
                     xmm8 += a4 * b2;
                  }

                  c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
                  c(ibegin+i+2UL,jj+j    ) += sum( xmm5 ) * alpha;
                  c(ibegin+i+2UL,jj+j+1UL) += sum( xmm6 ) * alpha;
                  c(ibegin+i+3UL,jj+j    ) += sum( xmm7 ) * alpha;
                  c(ibegin+i+3UL,jj+j+1UL) += sum( xmm8 ) * alpha;
               }

               if( j<jend )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                  }

                  c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  c(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  c(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            if( jj > ibegin+i+1UL ) continue;

            const size_t jend( min( ibegin+i-jj+2UL, jblock ) );
            size_t j( 0UL );

            for( ; (j+4UL) <= jend; j+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );
                  const SIMDType b3( B2.load(k,j+2UL) );
                  const SIMDType b4( B2.load(k,j+3UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a1 * b3;
                  xmm4 += a1 * b4;
                  xmm5 += a2 * b1;
                  xmm6 += a2 * b2;
                  xmm7 += a2 * b3;
                  xmm8 += a2 * b4;
               }

               c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               c(ibegin+i    ,jj+j+2UL) += sum( xmm3 ) * alpha;
               c(ibegin+i    ,jj+j+3UL) += sum( xmm4 ) * alpha;
               c(ibegin+i+1UL,jj+j    ) += sum( xmm5 ) * alpha;
               c(ibegin+i+1UL,jj+j+1UL) += sum( xmm6 ) * alpha;
               c(ibegin+i+1UL,jj+j+2UL) += sum( xmm7 ) * alpha;
               c(ibegin+i+1UL,jj+j+3UL) += sum( xmm8 ) * alpha;
            }

            for( ; (j+2UL) <= jend; j+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               c(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
               c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( j<jend )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
               }

               c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
               c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
            }
         }

         if( i<isize && jj <= ibegin+i )
         {
            const size_t jend( min( ibegin+i-jj+2UL, jblock ) );
            size_t j( 0UL );

            for( ; (j+2UL) <= jend; j+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               c(ibegin+i,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i,jj+j+1UL) += sum( xmm2 ) * alpha;
            }

            if( j<jend )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j);
               }

               c(ibegin+i,jj+j) += sum( xmm1 ) * alpha;
            }
         }

         jj += jblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t isize ( M - ibegin );

      A2 = serial( submatrix( A, ibegin, kk, isize, ksize, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( IsUpper_v<MT3> && jj+jblock <= kk ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix( B, kk, jj, ksize, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               if( jj > ibegin+i+4UL ) continue;

               const size_t jend( min( ibegin+i-jj+5UL, jblock ) );
               size_t j( 0UL );

               for( ; (j+2UL) <= jend; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+4UL,jj+j    ) += A2(i+4UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+4UL,jj+j+1UL) += A2(i+4UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jend ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+4UL,jj+j) += A2(i+4UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               if( jj > ibegin+i+3UL ) continue;

               const size_t jend( min( ibegin+i-jj+4UL, jblock ) );
               size_t j( 0UL );

               for( ; (j+2UL) <= jend; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jend ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            if( jj > ibegin+i+1UL ) continue;

            const size_t jend( min( ibegin+i-jj+2UL, jblock ) );
            size_t j( 0UL );

            for( ; (j+2UL) <= jend; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jend ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                  c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }
         }

         if( i<isize && jj <= ibegin+i )
         {
            const size_t jend( min( ibegin+i-jj+2UL, jblock ) );
            size_t j( 0UL );

            for( ; (j+2UL) <= jend; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i,jj+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  c(ibegin+i,jj+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jend ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i,jj+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         jj += jblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a lower dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a lower dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// column-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void lmmm( DenseMatrix<MT1,true>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE             ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNILOWER_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_LOWER_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE          ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t IBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( IBLOCK >= SIMDSIZE && IBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( IBLOCK, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, N );

   decltype(auto) c( derestrict( *C ) );

   if( isDefault( beta ) ) {
      reset( c );
   }
   else if( !isOne( beta ) ) {
      c *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jend  ( IsLower_v<MT3> ? kk+kblock : N );
      const size_t jsize ( jend - jbegin );

      B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jbegin, kblock, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( ( IsLower_v<MT2> && ii+iblock <= kk ) ||
             ( IsUpper_v<MT2> && kk+kblock <= ii ) ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ii, kk, iblock, kblock, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET3> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               if( ii+iblock < jbegin ) continue;

               size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

               for( ; (i+2UL) <= iblock; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );
                     const SIMDType b5( B2.load(k,j+4UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a1 * b3;
                     xmm4  += a1 * b4;
                     xmm5  += a1 * b5;
                     xmm6  += a2 * b1;
                     xmm7  += a2 * b2;
                     xmm8  += a2 * b3;
                     xmm9  += a2 * b4;
                     xmm10 += a2 * b5;
                  }

                  c(ii+i    ,jbegin+j    ) += sum( xmm1  ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += sum( xmm2  ) * alpha;
                  c(ii+i    ,jbegin+j+2UL) += sum( xmm3  ) * alpha;
                  c(ii+i    ,jbegin+j+3UL) += sum( xmm4  ) * alpha;
                  c(ii+i    ,jbegin+j+4UL) += sum( xmm5  ) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += sum( xmm6  ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += sum( xmm7  ) * alpha;
                  c(ii+i+1UL,jbegin+j+2UL) += sum( xmm8  ) * alpha;
                  c(ii+i+1UL,jbegin+j+3UL) += sum( xmm9  ) * alpha;
                  c(ii+i+1UL,jbegin+j+4UL) += sum( xmm10 ) * alpha;
               }

               if( i<iblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                     xmm5 += a1 * B2.load(k,j+4UL);
                  }

                  c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  c(ii+i,jbegin+j+4UL) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               if( ii+iblock < jbegin ) continue;

               size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

               for( ; (i+2UL) <= iblock; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a1 * b3;
                     xmm4 += a1 * b4;
                     xmm5 += a2 * b1;
                     xmm6 += a2 * b2;
                     xmm7 += a2 * b3;
                     xmm8 += a2 * b4;
                  }

                  c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i    ,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i    ,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += sum( xmm5 ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
                  c(ii+i+1UL,jbegin+j+2UL) += sum( xmm7 ) * alpha;
                  c(ii+i+1UL,jbegin+j+3UL) += sum( xmm8 ) * alpha;
               }

               if( i<iblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                  }

                  c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            if( ii+iblock < jbegin ) continue;

            size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

            for( ; (i+4UL) <= iblock; i+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );
                  const SIMDType a3( A2.load(i+2UL,k) );
                  const SIMDType a4( A2.load(i+3UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
                  xmm5 += a3 * b1;
                  xmm6 += a3 * b2;
                  xmm7 += a4 * b1;
                  xmm8 += a4 * b2;
               }

               c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               c(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               c(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
               c(ii+i+2UL,jbegin+j    ) += sum( xmm5 ) * alpha;
               c(ii+i+2UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
               c(ii+i+3UL,jbegin+j    ) += sum( xmm7 ) * alpha;
               c(ii+i+3UL,jbegin+j+1UL) += sum( xmm8 ) * alpha;
            }

            for( ; (i+2UL) <= iblock; i+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               c(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               c(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( i<iblock )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
            }
         }

         if( j<jsize && ii+iblock >= jbegin )
         {
            size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

            for( ; (i+2UL) <= iblock; i+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += A2.load(i    ,k) * b1;
                  xmm2 += A2.load(i+1UL,k) * b1;
               }

               c(ii+i    ,jbegin+j) += sum( xmm1 ) * alpha;
               c(ii+i+1UL,jbegin+j) += sum( xmm2 ) * alpha;
            }

            if( i<iblock )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  xmm1 += A2.load(i,k) * B2.load(k,j);
               }

               c(ii+i,jbegin+j) += sum( xmm1 ) * alpha;
            }
         }

         ii += iblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jsize ( N - jbegin );

      B2 = serial( submatrix( B, kk, jbegin, ksize, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( IsLower_v<MT2> && ii+iblock <= kk ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix( A, ii, kk, iblock, ksize, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               if( ii+iblock < jbegin ) continue;

               size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

               for( ; (i+2UL) <= iblock; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     c(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     c(ii+i    ,jbegin+j+4UL) += A2(i    ,k) * B2(k,j+4UL) * alpha;
                     c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     c(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                     c(ii+i+1UL,jbegin+j+4UL) += A2(i+1UL,k) * B2(k,j+4UL) * alpha;
                  }
               }

               if( i<iblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     c(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     c(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                     c(ii+i,jbegin+j+4UL) += A2(i,k) * B2(k,j+4UL) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               if( ii+iblock < jbegin ) continue;

               size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

               for( ; (i+2UL) <= iblock; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     c(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     c(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                  }
               }

               if( i<iblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     c(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     c(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                  }
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            if( ii+iblock < jbegin ) continue;

            size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

            for( ; (i+2UL) <= iblock; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( i<iblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }
         }

         if( j<jsize )
         {
            if( ii+iblock < jbegin ) continue;

            size_t i( ( ii > jbegin+j )?( 0UL ):( jbegin+j-ii ) );

            for( ; (i+2UL) <= iblock; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i    ,jbegin+j) += A2(i    ,k) * B2(k,j) * alpha;
                  c(ii+i+1UL,jbegin+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }

            if( i<iblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i,jbegin+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         ii += iblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a lower dense matrix/dense matrix multiplication (\f$ C=A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \return void
//
// This function implements the compute kernel for a lower dense matrix/dense matrix
// multiplication of the form \f$ C=A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must
// provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3 >
inline void lmmm( MT1& C, const MT2& A, const MT3& B )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   lmmm( C, A, B, ET1(1), ET1(0) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UPPER DENSE MATRIX MULTIPLICATION KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a upper dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side row-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a upper dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// row-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void ummm( DenseMatrix<MT1,false>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE             ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE          ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t JBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( JBLOCK >= SIMDSIZE && JBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( M, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, JBLOCK );

   decltype(auto) c( derestrict( *C ) );

   if( isDefault( beta ) ) {
      reset( c );
   }
   else if( !isOne( beta ) ) {
      c *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t iend  ( IsUpper_v<MT2> ? kk+kblock : M );
      const size_t isize ( iend - ibegin );

      A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ibegin, kk, isize, kblock, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( ( IsLower_v<MT3> && kk+kblock <= jj ) ||
             ( IsUpper_v<MT3> && jj+jblock <= kk ) ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jj, kblock, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               if( jj+jblock < ibegin ) continue;

               size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

               for( ; (j+2UL) <= jblock; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a2 * b1;
                     xmm4  += a2 * b2;
                     xmm5  += a3 * b1;
                     xmm6  += a3 * b2;
                     xmm7  += a4 * b1;
                     xmm8  += a4 * b2;
                     xmm9  += a5 * b1;
                     xmm10 += a5 * b2;
                  }

                  c(ibegin+i    ,jj+j    ) += sum( xmm1  ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += sum( xmm2  ) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += sum( xmm3  ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4  ) * alpha;
                  c(ibegin+i+2UL,jj+j    ) += sum( xmm5  ) * alpha;
                  c(ibegin+i+2UL,jj+j+1UL) += sum( xmm6  ) * alpha;
                  c(ibegin+i+3UL,jj+j    ) += sum( xmm7  ) * alpha;
                  c(ibegin+i+3UL,jj+j+1UL) += sum( xmm8  ) * alpha;
                  c(ibegin+i+4UL,jj+j    ) += sum( xmm9  ) * alpha;
                  c(ibegin+i+4UL,jj+j+1UL) += sum( xmm10 ) * alpha;
               }

               if( j<jblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );
                     const SIMDType a5( A2.load(i+4UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                     xmm5 += a5 * b1;
                  }

                  c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  c(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  c(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
                  c(ibegin+i+4UL,jj+j) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               if( jj+jblock < ibegin ) continue;

               size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

               for( ; (j+2UL) <= jblock; j+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a2 * b1;
                     xmm4 += a2 * b2;
                     xmm5 += a3 * b1;
                     xmm6 += a3 * b2;
                     xmm7 += a4 * b1;
                     xmm8 += a4 * b2;
                  }

                  c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
                  c(ibegin+i+2UL,jj+j    ) += sum( xmm5 ) * alpha;
                  c(ibegin+i+2UL,jj+j+1UL) += sum( xmm6 ) * alpha;
                  c(ibegin+i+3UL,jj+j    ) += sum( xmm7 ) * alpha;
                  c(ibegin+i+3UL,jj+j+1UL) += sum( xmm8 ) * alpha;
               }

               if( j<jblock )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );
                     const SIMDType a3( A2.load(i+2UL,k) );
                     const SIMDType a4( A2.load(i+3UL,k) );

                     const SIMDType b1( B2.load(k,j) );

                     xmm1 += a1 * b1;
                     xmm2 += a2 * b1;
                     xmm3 += a3 * b1;
                     xmm4 += a4 * b1;
                  }

                  c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
                  c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
                  c(ibegin+i+2UL,jj+j) += sum( xmm3 ) * alpha;
                  c(ibegin+i+3UL,jj+j) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            if( jj+jblock < ibegin ) continue;

            size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

            for( ; (j+4UL) <= jblock; j+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );
                  const SIMDType b3( B2.load(k,j+2UL) );
                  const SIMDType b4( B2.load(k,j+3UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a1 * b3;
                  xmm4 += a1 * b4;
                  xmm5 += a2 * b1;
                  xmm6 += a2 * b2;
                  xmm7 += a2 * b3;
                  xmm8 += a2 * b4;
               }

               c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               c(ibegin+i    ,jj+j+2UL) += sum( xmm3 ) * alpha;
               c(ibegin+i    ,jj+j+3UL) += sum( xmm4 ) * alpha;
               c(ibegin+i+1UL,jj+j    ) += sum( xmm5 ) * alpha;
               c(ibegin+i+1UL,jj+j+1UL) += sum( xmm6 ) * alpha;
               c(ibegin+i+1UL,jj+j+2UL) += sum( xmm7 ) * alpha;
               c(ibegin+i+1UL,jj+j+3UL) += sum( xmm8 ) * alpha;
            }

            for( ; (j+2UL) <= jblock; j+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               c(ibegin+i    ,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i    ,jj+j+1UL) += sum( xmm2 ) * alpha;
               c(ibegin+i+1UL,jj+j    ) += sum( xmm3 ) * alpha;
               c(ibegin+i+1UL,jj+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( j<jblock )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += a1 * b1;
                  xmm2 += a2 * b1;
               }

               c(ibegin+i    ,jj+j) += sum( xmm1 ) * alpha;
               c(ibegin+i+1UL,jj+j) += sum( xmm2 ) * alpha;
            }
         }

         if( i<isize && jj+jblock >= ibegin )
         {
            size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

            for( ; (j+2UL) <= jblock; j+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               c(ibegin+i,jj+j    ) += sum( xmm1 ) * alpha;
               c(ibegin+i,jj+j+1UL) += sum( xmm2 ) * alpha;
            }

            if( j<jblock )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j);
               }

               c(ibegin+i,jj+j) += sum( xmm1 ) * alpha;
            }
         }

         jj += jblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t ibegin( IsLower_v<MT2> ? kk : 0UL );
      const size_t isize ( M - ibegin );

      A2 = serial( submatrix( A, ibegin, kk, isize, ksize, unchecked ) );

      size_t jj( 0UL );
      size_t jblock( 0UL );

      while( jj < N )
      {
         jblock = ( ( jj+JBLOCK <= N )?( JBLOCK ):( N - jj ) );

         if( IsUpper_v<MT3> && jj+jblock <= kk ) {
            jj += jblock;
            continue;
         }

         B2 = serial( submatrix( B, kk, jj, ksize, jblock, unchecked ) );

         size_t i( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (i+5UL) <= isize; i+=5UL )
            {
               if( jj+jblock < ibegin ) continue;

               size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

               for( ; (j+2UL) <= jblock; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+4UL,jj+j    ) += A2(i+4UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+4UL,jj+j+1UL) += A2(i+4UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+4UL,jj+j) += A2(i+4UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (i+4UL) <= isize; i+=4UL )
            {
               if( jj+jblock < ibegin ) continue;

               size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

               for( ; (j+2UL) <= jblock; j+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+2UL,jj+j    ) += A2(i+2UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+2UL,jj+j+1UL) += A2(i+2UL,k) * B2(k,j+1UL) * alpha;
                     c(ibegin+i+3UL,jj+j    ) += A2(i+3UL,k) * B2(k,j    ) * alpha;
                     c(ibegin+i+3UL,jj+j+1UL) += A2(i+3UL,k) * B2(k,j+1UL) * alpha;
                  }
               }

               if( j<jblock ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                     c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+2UL,jj+j) += A2(i+2UL,k) * B2(k,j) * alpha;
                     c(ibegin+i+3UL,jj+j) += A2(i+3UL,k) * B2(k,j) * alpha;
                  }
               }
            }
         }

         for( ; (i+2UL) <= isize; i+=2UL )
         {
            if( jj+jblock < ibegin ) continue;

            size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

            for( ; (j+2UL) <= jblock; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i    ,jj+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  c(ibegin+i    ,jj+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  c(ibegin+i+1UL,jj+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  c(ibegin+i+1UL,jj+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i    ,jj+j) += A2(i    ,k) * B2(k,j) * alpha;
                  c(ibegin+i+1UL,jj+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }
         }

         if( i<isize && jj+jblock >= ibegin )
         {
            size_t j( ( jj > ibegin+i )?( 0UL ):( ibegin+i-jj ) );

            for( ; (j+2UL) <= jblock; j+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i,jj+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  c(ibegin+i,jj+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( j<jblock ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ibegin+i,jj+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         jj += jblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a upper dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B+\beta*C \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \param beta The scaling factor for \f$ C \f$.
// \return void
//
// This function implements the compute kernel for a upper dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B+\beta*C \f$. Both \a A and \a B must
// be non-expression dense matrix types, \a C must be a non-expression, non-adaptor,
// column-major dense matrix type. The element types of all three matrices must be SIMD
// combinable, i.e. must provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void ummm( DenseMatrix<MT1,true>& C, const MT2& A, const MT3& B, ST alpha, ST beta )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;
   using SIMDType = SIMDTrait_t<ET1>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE             ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE         ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIUPPER_MATRIX_TYPE      ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_UPPER_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE          ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   constexpr size_t SIMDSIZE( SIMDTrait<ET1>::size );

   constexpr bool remainder( !IsPadded_v<MT2> || !IsPadded_v<MT3> );

   constexpr size_t KBLOCK( MMM_OUTER_BLOCK_SIZE * ( 16UL/sizeof(ET1) ) );
   constexpr size_t IBLOCK( MMM_INNER_BLOCK_SIZE );

   BLAZE_STATIC_ASSERT( KBLOCK >= SIMDSIZE && KBLOCK % SIMDSIZE == 0UL );
   BLAZE_STATIC_ASSERT( IBLOCK >= SIMDSIZE && IBLOCK % SIMDSIZE == 0UL );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );
   const size_t K( A.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   DynamicMatrix<ET2,false> A2( IBLOCK, KBLOCK );
   DynamicMatrix<ET3,true>  B2( KBLOCK, N );

   decltype(auto) c( derestrict( *C ) );

   if( isDefault( beta ) ) {
      reset( c );
   }
   else if( !isOne( beta ) ) {
      c *= beta;
   }

   size_t kk( 0UL );
   size_t kblock( 0UL );

   while( kk + ( remainder ? SIMDSIZE-1UL : 0UL ) < K )
   {
      if( remainder ) {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( prevMultiple( K - kk, SIMDSIZE ) ) );
      }
      else {
         kblock = ( ( kk+KBLOCK <= K )?( KBLOCK ):( K - kk ) );
      }

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jend  ( IsLower_v<MT3> ? kk+kblock : N );
      const size_t jsize ( jend - jbegin );

      B2 = serial( submatrix< remainder ? unaligned : aligned >( B, kk, jbegin, kblock, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( ( IsLower_v<MT2> && ii+iblock <= kk ) ||
             ( IsUpper_v<MT2> && kk+kblock <= ii ) ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix< remainder ? unaligned : aligned >( A, ii, kk, iblock, kblock, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET3> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               if( ii > jbegin+j+4UL ) continue;

               const size_t iend( min( iblock, jbegin+j-ii+5UL ) );
               size_t i( 0UL );

               for( ; (i+2UL) <= iend; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8, xmm9, xmm10;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );
                     const SIMDType b5( B2.load(k,j+4UL) );

                     xmm1  += a1 * b1;
                     xmm2  += a1 * b2;
                     xmm3  += a1 * b3;
                     xmm4  += a1 * b4;
                     xmm5  += a1 * b5;
                     xmm6  += a2 * b1;
                     xmm7  += a2 * b2;
                     xmm8  += a2 * b3;
                     xmm9  += a2 * b4;
                     xmm10 += a2 * b5;
                  }

                  c(ii+i    ,jbegin+j    ) += sum( xmm1  ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += sum( xmm2  ) * alpha;
                  c(ii+i    ,jbegin+j+2UL) += sum( xmm3  ) * alpha;
                  c(ii+i    ,jbegin+j+3UL) += sum( xmm4  ) * alpha;
                  c(ii+i    ,jbegin+j+4UL) += sum( xmm5  ) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += sum( xmm6  ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += sum( xmm7  ) * alpha;
                  c(ii+i+1UL,jbegin+j+2UL) += sum( xmm8  ) * alpha;
                  c(ii+i+1UL,jbegin+j+3UL) += sum( xmm9  ) * alpha;
                  c(ii+i+1UL,jbegin+j+4UL) += sum( xmm10 ) * alpha;
               }

               if( i<iend )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                     xmm5 += a1 * B2.load(k,j+4UL);
                  }

                  c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  c(ii+i,jbegin+j+4UL) += sum( xmm5 ) * alpha;
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               if( ii > jbegin+j+3UL ) continue;

               const size_t iend( min( iblock, jbegin+j-ii+4UL ) );
               size_t i( 0UL );

               for( ; (i+2UL) <= iend; i+=2UL )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i    ,k) );
                     const SIMDType a2( A2.load(i+1UL,k) );

                     const SIMDType b1( B2.load(k,j    ) );
                     const SIMDType b2( B2.load(k,j+1UL) );
                     const SIMDType b3( B2.load(k,j+2UL) );
                     const SIMDType b4( B2.load(k,j+3UL) );

                     xmm1 += a1 * b1;
                     xmm2 += a1 * b2;
                     xmm3 += a1 * b3;
                     xmm4 += a1 * b4;
                     xmm5 += a2 * b1;
                     xmm6 += a2 * b2;
                     xmm7 += a2 * b3;
                     xmm8 += a2 * b4;
                  }

                  c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i    ,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i    ,jbegin+j+3UL) += sum( xmm4 ) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += sum( xmm5 ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
                  c(ii+i+1UL,jbegin+j+2UL) += sum( xmm7 ) * alpha;
                  c(ii+i+1UL,jbegin+j+3UL) += sum( xmm8 ) * alpha;
               }

               if( i<iend )
               {
                  SIMDType xmm1, xmm2, xmm3, xmm4;

                  for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
                  {
                     const SIMDType a1( A2.load(i,k) );

                     xmm1 += a1 * B2.load(k,j    );
                     xmm2 += a1 * B2.load(k,j+1UL);
                     xmm3 += a1 * B2.load(k,j+2UL);
                     xmm4 += a1 * B2.load(k,j+3UL);
                  }

                  c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
                  c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
                  c(ii+i,jbegin+j+2UL) += sum( xmm3 ) * alpha;
                  c(ii+i,jbegin+j+3UL) += sum( xmm4 ) * alpha;
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            if( ii > jbegin+j+1UL ) continue;

            const size_t iend( min( iblock, jbegin+j-ii+2UL ) );
            size_t i( 0UL );

            for( ; (i+4UL) <= iend; i+=4UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4, xmm5, xmm6, xmm7, xmm8;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );
                  const SIMDType a3( A2.load(i+2UL,k) );
                  const SIMDType a4( A2.load(i+3UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
                  xmm5 += a3 * b1;
                  xmm6 += a3 * b2;
                  xmm7 += a4 * b1;
                  xmm8 += a4 * b2;
               }

               c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               c(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               c(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
               c(ii+i+2UL,jbegin+j    ) += sum( xmm5 ) * alpha;
               c(ii+i+2UL,jbegin+j+1UL) += sum( xmm6 ) * alpha;
               c(ii+i+3UL,jbegin+j    ) += sum( xmm7 ) * alpha;
               c(ii+i+3UL,jbegin+j+1UL) += sum( xmm8 ) * alpha;
            }

            for( ; (i+2UL) <= iend; i+=2UL )
            {
               SIMDType xmm1, xmm2, xmm3, xmm4;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i    ,k) );
                  const SIMDType a2( A2.load(i+1UL,k) );

                  const SIMDType b1( B2.load(k,j    ) );
                  const SIMDType b2( B2.load(k,j+1UL) );

                  xmm1 += a1 * b1;
                  xmm2 += a1 * b2;
                  xmm3 += a2 * b1;
                  xmm4 += a2 * b2;
               }

               c(ii+i    ,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i    ,jbegin+j+1UL) += sum( xmm2 ) * alpha;
               c(ii+i+1UL,jbegin+j    ) += sum( xmm3 ) * alpha;
               c(ii+i+1UL,jbegin+j+1UL) += sum( xmm4 ) * alpha;
            }

            if( i<iend )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType a1( A2.load(i,k) );

                  xmm1 += a1 * B2.load(k,j    );
                  xmm2 += a1 * B2.load(k,j+1UL);
               }

               c(ii+i,jbegin+j    ) += sum( xmm1 ) * alpha;
               c(ii+i,jbegin+j+1UL) += sum( xmm2 ) * alpha;
            }
         }

         if( j<jsize && ii <= jbegin+j )
         {
            const size_t iend( min( iblock, jbegin+j-ii+2UL ) );
            size_t i( 0UL );

            for( ; (i+2UL) <= iend; i+=2UL )
            {
               SIMDType xmm1, xmm2;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  const SIMDType b1( B2.load(k,j) );

                  xmm1 += A2.load(i    ,k) * b1;
                  xmm2 += A2.load(i+1UL,k) * b1;
               }

               c(ii+i    ,jbegin+j) += sum( xmm1 ) * alpha;
               c(ii+i+1UL,jbegin+j) += sum( xmm2 ) * alpha;
            }

            if( i<iend )
            {
               SIMDType xmm1;

               for( size_t k=0UL; k<kblock; k+=SIMDSIZE )
               {
                  xmm1 += A2.load(i,k) * B2.load(k,j);
               }

               c(ii+i,jbegin+j) += sum( xmm1 ) * alpha;
            }
         }

         ii += iblock;
      }

      kk += kblock;
   }

   if( remainder && kk < K )
   {
      const size_t ksize( K - kk );

      const size_t jbegin( IsUpper_v<MT3> ? kk : 0UL );
      const size_t jsize ( N - jbegin );

      B2 = serial( submatrix( B, kk, jbegin, ksize, jsize, unchecked ) );

      size_t ii( 0UL );
      size_t iblock( 0UL );

      while( ii < M )
      {
         iblock = ( ( ii+IBLOCK <= M )?( IBLOCK ):( M - ii ) );

         if( IsLower_v<MT2> && ii+iblock <= kk ) {
            ii += iblock;
            continue;
         }

         A2 = serial( submatrix( A, ii, kk, iblock, ksize, unchecked ) );

         size_t j( 0UL );

         if( IsFloatingPoint_v<ET1> )
         {
            for( ; (j+5UL) <= jsize; j+=5UL )
            {
               if( ii > jbegin+j+4UL ) continue;

               const size_t iend( min( iblock, jbegin+j-ii+5UL ) );
               size_t i( 0UL );

               for( ; (i+2UL) <= iend; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     c(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     c(ii+i    ,jbegin+j+4UL) += A2(i    ,k) * B2(k,j+4UL) * alpha;
                     c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     c(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                     c(ii+i+1UL,jbegin+j+4UL) += A2(i+1UL,k) * B2(k,j+4UL) * alpha;
                  }
               }

               if( i<iend ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     c(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     c(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                     c(ii+i,jbegin+j+4UL) += A2(i,k) * B2(k,j+4UL) * alpha;
                  }
               }
            }
         }
         else
         {
            for( ; (j+4UL) <= jsize; j+=4UL )
            {
               if( ii > jbegin+j+3UL ) continue;

               const size_t iend( min( iblock, jbegin+j-ii+4UL ) );
               size_t i( 0UL );

               for( ; (i+2UL) <= iend; i+=2UL ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                     c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                     c(ii+i    ,jbegin+j+2UL) += A2(i    ,k) * B2(k,j+2UL) * alpha;
                     c(ii+i    ,jbegin+j+3UL) += A2(i    ,k) * B2(k,j+3UL) * alpha;
                     c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                     c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
                     c(ii+i+1UL,jbegin+j+2UL) += A2(i+1UL,k) * B2(k,j+2UL) * alpha;
                     c(ii+i+1UL,jbegin+j+3UL) += A2(i+1UL,k) * B2(k,j+3UL) * alpha;
                  }
               }

               if( i<iend ) {
                  for( size_t k=0UL; k<ksize; ++k ) {
                     c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                     c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
                     c(ii+i,jbegin+j+2UL) += A2(i,k) * B2(k,j+2UL) * alpha;
                     c(ii+i,jbegin+j+3UL) += A2(i,k) * B2(k,j+3UL) * alpha;
                  }
               }
            }
         }

         for( ; (j+2UL) <= jsize; j+=2UL )
         {
            if( ii > jbegin+j+1UL ) continue;

            const size_t iend( min( iblock, jbegin+j-ii+2UL ) );
            size_t i( 0UL );

            for( ; (i+2UL) <= iend; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i    ,jbegin+j    ) += A2(i    ,k) * B2(k,j    ) * alpha;
                  c(ii+i    ,jbegin+j+1UL) += A2(i    ,k) * B2(k,j+1UL) * alpha;
                  c(ii+i+1UL,jbegin+j    ) += A2(i+1UL,k) * B2(k,j    ) * alpha;
                  c(ii+i+1UL,jbegin+j+1UL) += A2(i+1UL,k) * B2(k,j+1UL) * alpha;
               }
            }

            if( i<iend ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i,jbegin+j    ) += A2(i,k) * B2(k,j    ) * alpha;
                  c(ii+i,jbegin+j+1UL) += A2(i,k) * B2(k,j+1UL) * alpha;
               }
            }
         }

         if( j<jsize && ii <= jbegin+j )
         {
            const size_t iend( min( iblock, jbegin+j-ii+2UL ) );
            size_t i( 0UL );

            for( ; (i+2UL) <= iend; i+=2UL ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i    ,jbegin+j) += A2(i    ,k) * B2(k,j) * alpha;
                  c(ii+i+1UL,jbegin+j) += A2(i+1UL,k) * B2(k,j) * alpha;
               }
            }

            if( i<iend ) {
               for( size_t k=0UL; k<ksize; ++k ) {
                  c(ii+i,jbegin+j) += A2(i,k) * B2(k,j) * alpha;
               }
            }
         }

         ii += iblock;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a upper dense matrix/dense matrix multiplication (\f$ C=A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \return void
//
// This function implements the compute kernel for a upper dense matrix/dense matrix
// multiplication of the form \f$ C=A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must
// provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3 >
inline void ummm( MT1& C, const MT2& A, const MT3& B )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   ummm( C, A, B, ET1(1), ET1(0) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SYMMETRIC DENSE MATRIX MULTIPLICATION KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a symmetric dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side row-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \return void
//
// This function implements the compute kernel for a symmetric dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix type.
// The element types of all three matrices must be SIMD combinable, i.e. must provide a common
// SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void smmm( DenseMatrix<MT1,false>& C, const MT2& A, const MT3& B, ST alpha )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   lmmm( C, A, B, alpha, ST(0) );

   for( size_t ii=0UL; ii<M; ii+=BLOCK_SIZE )
   {
      const size_t iend( min( M, ii+BLOCK_SIZE ) );

      for( size_t i=ii; i<iend; ++i ) {
         for( size_t j=i+1UL; j<iend; ++j ) {
            (*C)(i,j) = (*C)(j,i);
         }
      }

      for( size_t jj=ii+BLOCK_SIZE; jj<N; jj+=BLOCK_SIZE ) {
         const size_t jend( min( N, jj+BLOCK_SIZE ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               (*C)(i,j) = (*C)(j,i);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a symmetric dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \return void
//
// This function implements the compute kernel for a symmetric dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, column-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must provide
// a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void smmm( DenseMatrix<MT1,true>& C, const MT2& A, const MT3& B, ST alpha )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE        ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   ummm( C, A, B, alpha, ST(0) );

   for( size_t jj=0UL; jj<N; jj+=BLOCK_SIZE )
   {
      const size_t jend( min( N, jj+BLOCK_SIZE ) );

      for( size_t j=jj; j<jend; ++j ) {
         for( size_t i=jj+1UL; i<jend; ++i ) {
            (*C)(i,j) = (*C)(j,i);
         }
      }

      for( size_t ii=jj+BLOCK_SIZE; ii<M; ii+=BLOCK_SIZE ) {
         const size_t iend( min( M, ii+BLOCK_SIZE ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               (*C)(i,j) = (*C)(j,i);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a symmetric dense matrix/dense matrix multiplication (\f$ C=A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \return void
//
// This function implements the compute kernel for a symmetric dense matrix/dense matrix
// multiplication of the form \f$ C=A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must
// provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3 >
inline void smmm( MT1& C, const MT2& A, const MT3& B )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   smmm( C, A, B, ET1(1) );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HERMITIAN DENSE MATRIX MULTIPLICATION KERNELS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a Hermitian dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side row-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \return void
//
// This function implements the compute kernel for a Hermitian dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix type.
// The element types of all three matrices must be SIMD combinable, i.e. must provide a common
// SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void hmmm( DenseMatrix<MT1,false>& C, const MT2& A, const MT3& B, ST alpha )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE     ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   lmmm( C, A, B, alpha, ST(0) );

   for( size_t ii=0UL; ii<M; ii+=BLOCK_SIZE )
   {
      const size_t iend( min( M, ii+BLOCK_SIZE ) );

      for( size_t i=ii; i<iend; ++i ) {
         for( size_t j=i+1UL; j<iend; ++j ) {
            (*C)(i,j) = conj( (*C)(j,i) );
         }
      }

      for( size_t jj=ii+BLOCK_SIZE; jj<N; jj+=BLOCK_SIZE ) {
         const size_t jend( min( N, jj+BLOCK_SIZE ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               (*C)(i,j) = conj( (*C)(j,i) );
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a Hermitian dense matrix/dense matrix multiplication
//        (\f$ C=\alpha*A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \param alpha The scaling factor for \f$ A*B \f$.
// \return void
//
// This function implements the compute kernel for a Hermitian dense matrix/dense matrix
// multiplication of the form \f$ C=\alpha*A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, column-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must provide
// a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3, typename ST >
void hmmm( DenseMatrix<MT1,true>& C, const MT2& A, const MT3& B, ST alpha )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT1 );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE        ( MT1 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT1 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT2 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT2 );

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE   ( MT3 );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT3 );

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   const size_t M( A.rows()    );
   const size_t N( B.columns() );

   BLAZE_INTERNAL_ASSERT( A.columns() == B.rows(), "Invalid matrix sizes detected" );

   ummm( C, A, B, alpha, ST(0) );

   for( size_t jj=0UL; jj<N; jj+=BLOCK_SIZE )
   {
      const size_t jend( min( N, jj+BLOCK_SIZE ) );

      for( size_t j=jj; j<jend; ++j ) {
         for( size_t i=jj+1UL; i<jend; ++i ) {
            (*C)(i,j) = conj( (*C)(j,i) );
         }
      }

      for( size_t ii=jj+BLOCK_SIZE; ii<M; ii+=BLOCK_SIZE ) {
         const size_t iend( min( M, ii+BLOCK_SIZE ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               (*C)(i,j) = conj( (*C)(j,i) );
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compute kernel for a Hermitian dense matrix/dense matrix multiplication (\f$ C=A*B \f$).
// \ingroup dense_matrix
//
// \param C The target left-hand side column-major dense matrix.
// \param A The left-hand side multiplication operand.
// \param B The right-hand side multiplication operand.
// \return void
//
// This function implements the compute kernel for a Hermitian dense matrix/dense matrix
// multiplication of the form \f$ C=A*B \f$. Both \a A and \a B must be non-expression
// dense matrix types, \a C must be a non-expression, non-adaptor, row-major dense matrix
// type. The element types of all three matrices must be SIMD combinable, i.e. must
// provide a common SIMD interface.
*/
template< typename MT1, typename MT2, typename MT3 >
inline void hmmm( MT1& C, const MT2& A, const MT3& B )
{
   using ET1 = ElementType_t<MT1>;
   using ET2 = ElementType_t<MT2>;
   using ET3 = ElementType_t<MT3>;

   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET2 );
   BLAZE_CONSTRAINT_MUST_BE_SIMD_COMBINABLE_TYPES( ET1, ET3 );

   hmmm( C, A, B, ET1(1) );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
