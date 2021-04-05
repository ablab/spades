//=================================================================================================
/*!
//  \file blaze/math/expressions/DMatDMatEqualExpr.h
//  \brief Header file for the dense matrix/dense matrix equality comparison expression
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

#ifndef _BLAZE_MATH_EXPRESSIONS_DMATDMATEQUALEXPR_H_
#define _BLAZE_MATH_EXPRESSIONS_DMATDMATEQUALEXPR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/typetraits/HasSIMDEqual.h>
#include <blaze/math/typetraits/IsPadded.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/Optimizations.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary helper struct for the dense matrix/dense matrix equality comparison.
// \ingroup dense_matrix
*/
template< typename MT1    // Type of the left-hand side dense matrix
        , typename MT2 >  // Type of the right-hand side dense matrix
struct DMatDMatEqualExprHelper
{
   //**Type definitions****************************************************************************
   //! Composite type of the left-hand side dense matrix expression.
   using CT1 = RemoveReference_t< CompositeType_t<MT1> >;

   //! Composite type of the right-hand side dense matrix expression.
   using CT2 = RemoveReference_t< CompositeType_t<MT2> >;
   //**********************************************************************************************

   //**********************************************************************************************
   static constexpr bool value =
      ( useOptimizedKernels &&
        CT1::simdEnabled &&
        CT2::simdEnabled &&
        HasSIMDEqual_v< ElementType_t<CT1>, ElementType_t<CT2> > );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL BINARY RELATIONAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default equality check of two row-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side dense matrix
inline auto equal( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< DMatDMatEqualExprHelper<MT1,MT2>::value, bool >
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( *lhs );
   CT2 B( *rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t i=0UL; i<A.rows(); ++i ) {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         if( !equal<RF>( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized equality check of two row-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side dense matrix
inline auto equal( const DenseMatrix<MT1,false>& lhs, const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< DMatDMatEqualExprHelper<MT1,MT2>::value, bool >
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;
   using XT1 = RemoveReference_t<CT1>;
   using XT2 = RemoveReference_t<CT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( *lhs );
   CT2 B( *rhs );

   constexpr size_t SIMDSIZE = SIMDTrait< ElementType_t<MT1> >::size;
   constexpr bool remainder( !IsPadded_v<XT1> || !IsPadded_v<XT2> );

   const size_t M( A.rows()    );
   const size_t N( A.columns() );

   const size_t jpos( remainder ? prevMultiple( N, SIMDSIZE ) : N );
   BLAZE_INTERNAL_ASSERT( jpos <= N, "Invalid end calculation" );

   for( size_t i=0UL; i<M; ++i )
   {
      size_t j( 0UL );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         if( !equal<RF>( A.load(i,j             ), B.load(i,j             ) ) ) return false;
         if( !equal<RF>( A.load(i,j+SIMDSIZE    ), B.load(i,j+SIMDSIZE    ) ) ) return false;
         if( !equal<RF>( A.load(i,j+SIMDSIZE*2UL), B.load(i,j+SIMDSIZE*2UL) ) ) return false;
         if( !equal<RF>( A.load(i,j+SIMDSIZE*3UL), B.load(i,j+SIMDSIZE*3UL) ) ) return false;
      }
      for( ; (j+SIMDSIZE) < jpos; j+=SIMDSIZE*2UL ) {
         if( !equal<RF>( A.load(i,j         ), B.load(i,j         ) ) ) return false;
         if( !equal<RF>( A.load(i,j+SIMDSIZE), B.load(i,j+SIMDSIZE) ) ) return false;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         if( !equal<RF>( A.load(i,j), B.load(i,j) ) ) return false;
      }
      for( ; remainder && j<A.columns(); ++j ) {
         if( !equal<RF>( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default equality check of two column-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side dense matrix
inline auto equal( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< DMatDMatEqualExprHelper<MT1,MT2>::value, bool >
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( *lhs );
   CT2 B( *rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   for( size_t j=0UL; j<A.columns(); ++j ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         if( !equal<RF>( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized equality check of two column-major dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , typename MT2 >     // Type of the right-hand side dense matrix
inline auto equal( const DenseMatrix<MT1,true>& lhs, const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< DMatDMatEqualExprHelper<MT1,MT2>::value, bool >
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;
   using XT1 = RemoveReference_t<CT1>;
   using XT2 = RemoveReference_t<CT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( *lhs );
   CT2 B( *rhs );

   constexpr size_t SIMDSIZE = SIMDTrait< ElementType_t<MT1> >::size;
   constexpr bool remainder( !IsPadded_v<XT1> || !IsPadded_v<XT2> );

   const size_t M( A.rows()    );
   const size_t N( A.columns() );

   const size_t ipos( remainder ? prevMultiple( M, SIMDSIZE ) : M );
   BLAZE_INTERNAL_ASSERT( ipos <= M, "Invalid end calculation" );

   for( size_t j=0UL; j<N; ++j )
   {
      size_t i( 0UL );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         if( !equal<RF>( A.load(i             ,j), B.load(i             ,j) ) ) return false;
         if( !equal<RF>( A.load(i+SIMDSIZE    ,j), B.load(i+SIMDSIZE    ,j) ) ) return false;
         if( !equal<RF>( A.load(i+SIMDSIZE*2UL,j), B.load(i+SIMDSIZE*2UL,j) ) ) return false;
         if( !equal<RF>( A.load(i+SIMDSIZE*3UL,j), B.load(i+SIMDSIZE*3UL,j) ) ) return false;
      }
      for( ; (i+SIMDSIZE) < ipos; i+=SIMDSIZE*2UL ) {
         if( !equal<RF>( A.load(i         ,j), B.load(i         ,j) ) ) return false;
         if( !equal<RF>( A.load(i+SIMDSIZE,j), B.load(i+SIMDSIZE,j) ) ) return false;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         if( !equal<RF>( A.load(i,j), B.load(i,j) ) ) return false;
      }
      for( ; remainder && i<M; ++i ) {
         if( !equal<RF>( A(i,j), B(i,j) ) ) return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Equality check of two dense matrices with different storage order.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
//
// Equal function for the comparison of two dense matrices. Due to the limited machine accuracy,
// a direct comparison of two floating point numbers should be avoided. This function offers the
// possibility to compare two floating-point matrices with a certain accuracy margin.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT1       // Type of the left-hand side dense matrix
        , typename MT2       // Type of the right-hand side dense matrix
        , bool SO >          // Storage order
inline bool equal( const DenseMatrix<MT1,SO>& lhs, const DenseMatrix<MT2,!SO>& rhs )
{
   using CT1 = CompositeType_t<MT1>;
   using CT2 = CompositeType_t<MT2>;

   // Early exit in case the matrix sizes don't match
   if( (*lhs).rows() != (*rhs).rows() || (*lhs).columns() != (*rhs).columns() )
      return false;

   // Evaluation of the two dense matrix operands
   CT1 A( *lhs );
   CT2 B( *rhs );

   // In order to compare the two matrices, the data values of the lower-order data
   // type are converted to the higher-order data type within the equal function.
   const size_t rows   ( A.rows() );
   const size_t columns( A.columns() );
   const size_t block  ( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows; ii+=block ) {
      const size_t iend( ( rows < ii+block )?( rows ):( ii+block ) );
      for( size_t jj=0UL; jj<columns; jj+=block ) {
         const size_t jend( ( columns < jj+block )?( columns ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               if( !equal<RF>( A(i,j), B(i,j) ) ) return false;
            }
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality operator for the comparison of two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side matrix for the comparison.
// \param rhs The right-hand side matrix for the comparison.
// \return \a true if the two matrices are equal, \a false if not.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator==( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality operator for the comparison of two dense matrices.
// \ingroup dense_matrix
//
// \param lhs The left-hand side dense matrix for the comparison.
// \param rhs The right-hand side dense matrix for the comparison.
// \return \a true if the two matrices are not equal, \a false if they are equal.
*/
template< typename MT1  // Type of the left-hand side dense matrix
        , bool SO1      // Storage order of the left-hand side dense matrix
        , typename MT2  // Type of the right-hand side dense matrix
        , bool SO2 >    // Storage order of the right-hand side dense matrix
inline bool operator!=( const DenseMatrix<MT1,SO1>& lhs, const DenseMatrix<MT2,SO2>& rhs )
{
   return !equal<relaxed>( lhs, rhs );
}
//*************************************************************************************************

} // namespace blaze

#endif
