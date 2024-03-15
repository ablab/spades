//=================================================================================================
/*!
//  \file blaze/math/lapack/orgl2.h
//  \brief Header file for the LAPACK functions to reconstruct Q from a LQ decomposition (orgl2)
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

#ifndef _BLAZE_MATH_LAPACK_ORGL2_H_
#define _BLAZE_MATH_LAPACK_ORGL2_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Contiguous.h>
#include <blaze/math/constraints/MutableDataAccess.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/lapack/clapack/orgl2.h>
#include <blaze/math/lapack/clapack/org2r.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Builtin.h>
#include <blaze/util/NumericCast.h>


namespace blaze {

//=================================================================================================
//
//  LAPACK FUNCTIONS TO RECONSTRUCT Q FROM A LQ DECOMPOSITION (ORGL2)
//
//=================================================================================================

//*************************************************************************************************
/*!\name LAPACK functions to reconstruct Q from a LQ decomposition (orgl2) */
//@{
template< typename MT, bool SO >
void orgl2( DenseMatrix<MT,SO>& A, const ElementType_t<MT>* tau );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief LAPACK kernel for the reconstruction of the orthogonal matrix Q from a LQ decomposition.
// \ingroup lapack_decomposition
//
// \param A The decomposed matrix.
// \param tau Array for the scalar factors of the elementary reflectors; size >= min( \a m, \a n ).
// \return void
//
// This function reconstructs the orthogonal matrix \a Q of a LQ decomposition based on the LAPACK
// orgl2() functions from matrices that have already been LQ factorized by the gelqf() functions.
// Note that this function can only be used for general, non-adapted matrices with \c float or
// \c double element type. The attempt to call the function with any adapted matrix or matrices
// of any other element type results in a compile time error!\n
//
// The min(\a m,\a n)-by-\a n \a Q matrix is stored within the given matrix \a A:

   \code
   using blaze::DynamicMatrix;
   using blaze::columnMajor;

   DynamicMatrix<double,columnMajor> A;
   DynamicVector<double> tau;
   // ... Resizing and initialization

   gelqf( A, tau.data() );  // Performing the LQ decomposition
   orgl2( A, tau.data() );  // Reconstructing the Q matrix

   const size_t m( A.rows() );
   const size_t n( A.columns() );

   DynamicMatrix<double,columnMajor> Q( submatrix( A, 0UL, 0UL, min(m,n), n ) );
   \endcode

// For more information on the orgl2() functions (i.e. sorgl2() and dorgl2()) see the LAPACK
// online documentation browser:
//
//        http://www.netlib.org/lapack/explore-html/
//
// \note This function can only be used if a fitting LAPACK library, which supports this function,
// is available and linked to the executable. Otherwise a call to this function will result in a
// linker error.
*/
template< typename MT, bool SO >
inline void orgl2( DenseMatrix<MT,SO>& A, const ElementType_t<MT>* tau )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_HAVE_MUTABLE_DATA_ACCESS( MT );
   BLAZE_CONSTRAINT_MUST_BE_CONTIGUOUS_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );
   BLAZE_CONSTRAINT_MUST_BE_BUILTIN_TYPE( ElementType_t<MT> );

   using ET = ElementType_t<MT>;

   blas_int_t m   ( numeric_cast<blas_int_t>( SO ? (*A).rows() : (*A).columns() ) );
   blas_int_t n   ( numeric_cast<blas_int_t>( SO ? (*A).columns() : (*A).rows() ) );
   blas_int_t k   ( min( m, n ) );
   blas_int_t lda ( numeric_cast<blas_int_t>( (*A).spacing() ) );
   blas_int_t info( 0 );

   if( k == 0 ) {
      return;
   }

   const std::unique_ptr<ET[]> work( new ET[m] );

   if( SO ) {
      orgl2( k, n, k, (*A).data(), lda, tau, work.get(), &info );
   }
   else {
      org2r( m, k, k, (*A).data(), lda, tau, work.get(), &info );
   }

   BLAZE_INTERNAL_ASSERT( info == 0, "Invalid argument for Q reconstruction" );
}
//*************************************************************************************************

} // namespace blaze

#endif
