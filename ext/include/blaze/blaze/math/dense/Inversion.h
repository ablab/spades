//=================================================================================================
/*!
//  \file blaze/math/dense/Inversion.h
//  \brief Header file for the dense matrix in-place inversion kernels
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

#ifndef _BLAZE_MATH_DENSE_INVERSION_H_
#define _BLAZE_MATH_DENSE_INVERSION_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <memory>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Adaptor.h>
#include <blaze/math/constraints/BLASCompatible.h>
#include <blaze/math/constraints/StrictlyTriangular.h>
#include <blaze/math/constraints/Uniform.h>
#include <blaze/math/dense/StaticMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/lapack/getrf.h>
#include <blaze/math/lapack/getri.h>
#include <blaze/math/lapack/hetrf.h>
#include <blaze/math/lapack/hetri.h>
#include <blaze/math/lapack/potrf.h>
#include <blaze/math/lapack/potri.h>
#include <blaze/math/lapack/sytrf.h>
#include <blaze/math/lapack/sytri.h>
#include <blaze/math/lapack/trtri.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/IsDivisor.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/Real.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  AUXILIARY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the proper inversion flag for the given matrix type \a MT.
//
// \return The proper inversion flag for the matrix type \a MT.
*/
template< typename MT >
constexpr InversionFlag getInversionFlag() noexcept
{
   return ( IsDiagonal_v<MT>  ? asDiagonal
          : IsUniUpper_v<MT>  ? asUniUpper
          : IsUpper_v<MT>     ? asUpper
          : IsUniLower_v<MT>  ? asUniLower
          : IsLower_v<MT>     ? asLower
          : IsSymmetric_v<MT> ? asSymmetric
          : IsHermitian_v<MT> ? asHermitian
          :                     asGeneral );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR 2x2 MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertGeneral2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET det( A(0,0)*A(1,1) - A(0,1)*A(1,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   A(0,0) =  A(1,1) * idet;
   A(1,0) = -A(1,0) * idet;
   A(0,1) = -A(0,1) * idet;
   A(1,1) =  a11;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The symmetric dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense \f$ 2 \times 2 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertSymmetric2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const MT& A( *dm );
   MT& B( *dm );

   const ET det( A(0,0)*A(1,1) - A(0,1)*A(1,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   B(0,0) =  A(1,1) * idet;
   B(1,0) = -A(1,0) * idet;
   B(0,1) =  B(1,0);
   B(1,1) =  a11;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given Hermitian dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The Hermitian dense matrix to be inverted.
// \return void
//
// This function inverts the given Hermitian dense \f$ 2 \times 2 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertHermitian2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const MT& A( *dm );
   MT& B( *dm );

   const ET det( real( A(0,0)*A(1,1) - A(0,1)*A(1,0) ) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   B(0,0) =  ET( A(1,1) * idet );
   B(1,0) = -A(1,0) * idet;
   B(0,1) =  conj( B(1,0) );
   B(1,1) =  ET( a11 );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The lower dense matrix to be inverted.
// \return void
//
// This function inverts the given lower dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLower2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET det( A(0,0) * A(1,1) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   A(0,0) =  A(1,1) * idet;
   A(1,0) = -A(1,0) * idet;
   A(1,1) =  a11;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given unilower dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The unilower dense matrix to be inverted.
// \return void
//
// This function inverts the given unilower dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLower2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   MT& A( *dm );

   A(1,0) = -A(1,0);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The upper dense matrix to be inverted.
// \return void
//
// This function inverts the given upper dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpper2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET det( A(0,0) * A(1,1) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   A(0,0) =  A(1,1) * idet;
   A(0,1) = -A(0,1) * idet;
   A(1,1) =  a11;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The uniupper dense matrix to be inverted.
// \return void
//
// This function inverts the given uniupper dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpper2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   MT& A( *dm );

   A(0,1) = -A(0,1);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense \f$ 2 \times 2 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The diagonal dense matrix to be inverted.
// \return void
//
// This function inverts the given diagonal dense \f$ 2 \times 2 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonal2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 2UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 2UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET det( A(0,0) * A(1,1) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );
   const ET a11( A(0,0) * idet );

   A(0,0) =  A(1,1) * idet;
   A(1,1) =  a11;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense 2x2 matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense 2x2 matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert2x2( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertGeneral2x2  ( *dm ); break;
      case byLDLT     : invertSymmetric2x2( *dm ); break;
      case byLDLH     : invertHermitian2x2( *dm ); break;
      case byLLH      : invertHermitian2x2( *dm ); break;
      case asGeneral  : invertGeneral2x2  ( *dm ); break;
      case asSymmetric: invertSymmetric2x2( *dm ); break;
      case asHermitian: invertHermitian2x2( *dm ); break;
      case asLower    : invertLower2x2    ( *dm ); break;
      case asUniLower : invertUniLower2x2 ( *dm ); break;
      case asUpper    : invertUpper2x2    ( *dm ); break;
      case asUniUpper : invertUniUpper2x2 ( *dm ); break;
      case asDiagonal : invertDiagonal2x2 ( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR 3x3 MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertGeneral3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   B(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
   B(1,0) = A(1,2)*A(2,0) - A(1,0)*A(2,2);
   B(2,0) = A(1,0)*A(2,1) - A(1,1)*A(2,0);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,1) = A(0,2)*A(2,1) - A(0,1)*A(2,2);
   B(1,1) = A(0,0)*A(2,2) - A(0,2)*A(2,0);
   B(2,1) = A(0,1)*A(2,0) - A(0,0)*A(2,1);
   B(0,2) = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   B(1,2) = A(0,2)*A(1,0) - A(0,0)*A(1,2);
   B(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The symmetric dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense \f$ 3 \times 3 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertSymmetric3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   B(0,0) = A(1,1)*A(2,2) - A(1,2)*A(2,1);
   B(1,0) = A(1,2)*A(2,0) - A(1,0)*A(2,2);
   B(2,0) = A(1,0)*A(2,1) - A(1,1)*A(2,0);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,1) = B(1,0);
   B(1,1) = A(0,0)*A(2,2) - A(0,2)*A(2,0);
   B(2,1) = A(0,1)*A(2,0) - A(0,0)*A(2,1);
   B(0,2) = B(2,0);
   B(1,2) = B(2,1);
   B(2,2) = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given Hermitian dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The Hermitian dense matrix to be inverted.
// \return void
//
// This function inverts the given Hermitian dense \f$ 3 \times 3 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertHermitian3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   B(0,0) = ET( real( A(1,1)*A(2,2) - A(1,2)*A(2,1) ) );
   B(1,0) = A(1,2)*A(2,0) - A(1,0)*A(2,2);
   B(2,0) = A(1,0)*A(2,1) - A(1,1)*A(2,0);

   const ET det( real( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) ) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,1) = conj( B(1,0) );
   B(1,1) = ET( real( A(0,0)*A(2,2) - A(0,2)*A(2,0) ) );
   B(2,1) = A(0,1)*A(2,0) - A(0,0)*A(2,1);
   B(0,2) = conj( B(2,0) );
   B(1,2) = conj( B(2,1) );
   B(2,2) = ET( real( A(0,0)*A(1,1) - A(0,1)*A(1,0) ) );

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The lower dense matrix to be inverted.
// \return void
//
// This function inverts the given lower dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLower3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp( A(1,1)*A(2,2) );
   const ET det( A(0,0)*tmp );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,0) =   tmp;
   B(1,0) = - A(1,0)*A(2,2);
   B(2,0) =   A(1,0)*A(2,1) - A(1,1)*A(2,0);
   B(1,1) =   A(0,0)*A(2,2);
   B(2,1) = - A(0,0)*A(2,1);
   B(2,2) =   A(0,0)*A(1,1);

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given unilower dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The unilower dense matrix to be inverted.
// \return void
//
// This function inverts the given unilower dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLower3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   B(1,0) = - A(1,0);
   B(2,0) =   A(1,0)*A(2,1) - A(2,0);
   B(2,1) = - A(2,1);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The upper dense matrix to be inverted.
// \return void
//
// This function inverts the given upper dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpper3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp( A(1,1)*A(2,2) );
   const ET det( A(0,0)*tmp );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,0) = tmp;
   B(0,1) = - A(0,1)*A(2,2);
   B(1,1) =   A(0,0)*A(2,2);
   B(0,2) =   A(0,1)*A(1,2) - A(0,2)*A(1,1);
   B(1,2) = - A(0,0)*A(1,2);
   B(2,2) =   A(0,0)*A(1,1);

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The uniupper dense matrix to be inverted.
// \return void
//
// This function inverts the given uniupper dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpper3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,3UL,3UL,SO> A( *dm );
   MT& B( *dm );

   B(0,1) = - A(0,1);
   B(0,2) =   A(0,1)*A(1,2) - A(0,2);
   B(1,2) = - A(1,2);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense \f$ 3 \times 3 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The diagonal dense matrix to be inverted.
// \return void
//
// This function inverts the given diagonal dense \f$ 3 \times 3 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonal3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 3UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 3UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET tmp1( A(0,0)*A(1,1) );
   const ET tmp2( A(0,0)*A(2,2) );

   const ET det( tmp1*A(2,2) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );

   A(0,0) = A(1,1)*A(2,2)*idet;
   A(1,1) = tmp2*idet;
   A(2,2) = tmp1*idet;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense 3x3 matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense 3x3 matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert3x3( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertGeneral3x3  ( *dm ); break;
      case byLDLT     : invertSymmetric3x3( *dm ); break;
      case byLDLH     : invertHermitian3x3( *dm ); break;
      case byLLH      : invertHermitian3x3( *dm ); break;
      case asGeneral  : invertGeneral3x3  ( *dm ); break;
      case asSymmetric: invertSymmetric3x3( *dm ); break;
      case asHermitian: invertHermitian3x3( *dm ); break;
      case asLower    : invertLower3x3    ( *dm ); break;
      case asUniLower : invertUniLower3x3 ( *dm ); break;
      case asUpper    : invertUpper3x3    ( *dm ); break;
      case asUniUpper : invertUniUpper3x3 ( *dm ); break;
      case asDiagonal : invertDiagonal3x3 ( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR 4x4 MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertGeneral4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1( A(2,2)*A(3,3) - A(2,3)*A(3,2) );
   ET tmp2( A(2,1)*A(3,3) - A(2,3)*A(3,1) );
   ET tmp3( A(2,1)*A(3,2) - A(2,2)*A(3,1) );

   B(0,0) = A(1,1)*tmp1 - A(1,2)*tmp2 + A(1,3)*tmp3;
   B(0,1) = A(0,2)*tmp2 - A(0,1)*tmp1 - A(0,3)*tmp3;

   ET tmp4( A(2,0)*A(3,3) - A(2,3)*A(3,0) );
   ET tmp5( A(2,0)*A(3,2) - A(2,2)*A(3,0) );

   B(1,0) = A(1,2)*tmp4 - A(1,0)*tmp1 - A(1,3)*tmp5;
   B(1,1) = A(0,0)*tmp1 - A(0,2)*tmp4 + A(0,3)*tmp5;

   tmp1 = A(2,0)*A(3,1) - A(2,1)*A(3,0);

   B(2,0) = A(1,0)*tmp2 - A(1,1)*tmp4 + A(1,3)*tmp1;
   B(2,1) = A(0,1)*tmp4 - A(0,0)*tmp2 - A(0,3)*tmp1;
   B(3,0) = A(1,1)*tmp5 - A(1,0)*tmp3 - A(1,2)*tmp1;
   B(3,1) = A(0,0)*tmp3 - A(0,1)*tmp5 + A(0,2)*tmp1;

   tmp1 = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2 = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3 = A(0,1)*A(1,2) - A(0,2)*A(1,1);

   B(0,2) = A(3,1)*tmp1 - A(3,2)*tmp2 + A(3,3)*tmp3;
   B(0,3) = A(2,2)*tmp2 - A(2,1)*tmp1 - A(2,3)*tmp3;

   tmp4 = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5 = A(0,0)*A(1,2) - A(0,2)*A(1,0);

   B(1,2) = A(3,2)*tmp4 - A(3,0)*tmp1 - A(3,3)*tmp5;
   B(1,3) = A(2,0)*tmp1 - A(2,2)*tmp4 + A(2,3)*tmp5;

   tmp1 = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B(2,2) = A(3,0)*tmp2 - A(3,1)*tmp4 + A(3,3)*tmp1;
   B(2,3) = A(2,1)*tmp4 - A(2,0)*tmp2 - A(2,3)*tmp1;
   B(3,2) = A(3,1)*tmp5 - A(3,0)*tmp3 - A(3,2)*tmp1;
   B(3,3) = A(2,0)*tmp3 - A(2,1)*tmp5 + A(2,2)*tmp1;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The symmetric dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense \f$ 4 \times 4 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertSymmetric4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1( A(2,2)*A(3,3) - A(2,3)*A(3,2) );
   ET tmp2( A(2,1)*A(3,3) - A(2,3)*A(3,1) );
   ET tmp3( A(2,1)*A(3,2) - A(2,2)*A(3,1) );

   B(0,0) = A(1,1)*tmp1 - A(1,2)*tmp2 + A(1,3)*tmp3;
   B(0,1) = A(0,2)*tmp2 - A(0,1)*tmp1 - A(0,3)*tmp3;

   ET tmp4( A(2,0)*A(3,3) - A(2,3)*A(3,0) );
   ET tmp5( A(2,0)*A(3,2) - A(2,2)*A(3,0) );

   B(1,1) = A(0,0)*tmp1 - A(0,2)*tmp4 + A(0,3)*tmp5;

   tmp1 = A(2,0)*A(3,1) - A(2,1)*A(3,0);

   B(2,0) = A(1,0)*tmp2 - A(1,1)*tmp4 + A(1,3)*tmp1;
   B(2,1) = A(0,1)*tmp4 - A(0,0)*tmp2 - A(0,3)*tmp1;
   B(3,0) = A(1,1)*tmp5 - A(1,0)*tmp3 - A(1,2)*tmp1;
   B(3,1) = A(0,0)*tmp3 - A(0,1)*tmp5 + A(0,2)*tmp1;

   tmp1 = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp2 = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp3 = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp4 = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp5 = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B(2,2) = A(3,0)*tmp1 - A(3,1)*tmp3 + A(3,3)*tmp5;
   B(2,3) = A(2,1)*tmp3 - A(2,0)*tmp1 - A(2,3)*tmp5;
   B(3,3) = A(2,0)*tmp2 - A(2,1)*tmp4 + A(2,2)*tmp5;

   B(0,2) = B(2,0);
   B(0,3) = B(3,0);
   B(1,0) = B(0,1);
   B(1,2) = B(2,1);
   B(1,3) = B(3,1);
   B(3,2) = B(2,3);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given Hermitian dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The Hermitian dense matrix to be inverted.
// \return void
//
// This function inverts the given Hermitian dense \f$ 4 \times 4 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertHermitian4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1( A(2,2)*A(3,3) - A(2,3)*A(3,2) );
   ET tmp2( A(2,1)*A(3,3) - A(2,3)*A(3,1) );
   ET tmp3( A(2,1)*A(3,2) - A(2,2)*A(3,1) );

   B(0,0) = ET( real( A(1,1)*tmp1 - A(1,2)*tmp2 + A(1,3)*tmp3 ) );
   B(0,1) = A(0,2)*tmp2 - A(0,1)*tmp1 - A(0,3)*tmp3;

   ET tmp4( A(2,0)*A(3,3) - A(2,3)*A(3,0) );
   ET tmp5( A(2,0)*A(3,2) - A(2,2)*A(3,0) );

   B(1,1) = ET( real( A(0,0)*tmp1 - A(0,2)*tmp4 + A(0,3)*tmp5 ) );

   tmp1 = A(2,0)*A(3,1) - A(2,1)*A(3,0);

   B(2,0) = A(1,0)*tmp2 - A(1,1)*tmp4 + A(1,3)*tmp1;
   B(2,1) = A(0,1)*tmp4 - A(0,0)*tmp2 - A(0,3)*tmp1;
   B(3,0) = A(1,1)*tmp5 - A(1,0)*tmp3 - A(1,2)*tmp1;
   B(3,1) = A(0,0)*tmp3 - A(0,1)*tmp5 + A(0,2)*tmp1;

   tmp1 = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp2 = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp3 = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp4 = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp5 = A(0,0)*A(1,1) - A(0,1)*A(1,0);

   B(2,2) = ET( real( A(3,0)*tmp1 - A(3,1)*tmp3 + A(3,3)*tmp5 ) );
   B(2,3) = A(2,1)*tmp3 - A(2,0)*tmp1 - A(2,3)*tmp5;
   B(3,3) = ET( real( A(2,0)*tmp2 - A(2,1)*tmp4 + A(2,2)*tmp5 ) );

   B(0,2) = conj( B(2,0) );
   B(0,3) = conj( B(3,0) );
   B(1,0) = conj( B(0,1) );
   B(1,2) = conj( B(2,1) );
   B(1,3) = conj( B(3,1) );
   B(3,2) = conj( B(2,3) );

   const ET det( real( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) ) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The lower dense matrix to be inverted.
// \return void
//
// This function inverts the given lower dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLower4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(2,2)*A(3,3) );
   const ET tmp2( A(2,1)*A(3,3) );
   const ET tmp3( A(2,1)*A(3,2) - A(2,2)*A(3,1) );
   const ET tmp4( A(0,0)*A(1,1) );

   const ET det( tmp4 * A(2,2) * A(3,3) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,0) =   A(1,1)*tmp1;
   B(1,0) = - A(1,0)*tmp1;
   B(2,0) =   A(1,0)*tmp2 - A(1,1)*A(2,0)*A(3,3);
   B(3,0) =   A(1,1)*( A(2,0)*A(3,2) - A(2,2)*A(3,0) ) - A(1,0)*tmp3;
   B(1,1) =   A(0,0)*tmp1;
   B(2,1) = - A(0,0)*tmp2;
   B(3,1) =   A(0,0)*tmp3;
   B(2,2) =   A(3,3)*tmp4;
   B(3,2) = - A(3,2)*tmp4;
   B(3,3) =   A(2,2)*tmp4;

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given unilower dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The unilower dense matrix to be inverted.
// \return void
//
// This function inverts the given unilower dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLower4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp( A(2,1)*A(3,2) - A(3,1) );

   B(1,0) = - A(1,0);
   B(2,0) =   A(1,0)*A(2,1) - A(2,0);
   B(3,0) =   A(2,0)*A(3,2) - A(3,0) - A(1,0)*tmp;
   B(2,1) = - A(2,1);
   B(3,1) =   tmp;
   B(3,2) = - A(3,2);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The upper dense matrix to be inverted.
// \return void
//
// This function inverts the given upper dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpper4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1( A(2,2)*A(3,3) );
   ET tmp2( A(0,1)*A(1,2) - A(0,2)*A(1,1) );
   ET tmp3( A(0,0)*A(1,2) );
   ET tmp4( A(0,0)*A(1,1) );

   const ET det( A(0,0)*A(1,1)*tmp1 );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B(0,0) =   A(1,1)*tmp1;
   B(0,1) = - A(0,1)*tmp1;
   B(1,1) =   A(0,0)*tmp1;
   B(0,2) =   A(3,3)*tmp2;
   B(1,2) = - A(3,3)*tmp3;
   B(2,2) =   A(3,3)*tmp4;
   B(0,3) =   A(2,2)*( A(0,1)*A(1,3) - A(0,3)*A(1,1) ) - A(2,3)*tmp2;
   B(1,3) =   A(2,3)*tmp3 - A(2,2)*A(0,0)*A(1,3);
   B(2,3) = - A(2,3)*tmp4;
   B(3,3) =   A(2,2)*tmp4;

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The uniupper dense matrix to be inverted.
// \return void
//
// This function inverts the given uniupper dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpper4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,4UL,4UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp( A(0,1)*A(1,2) - A(0,2) );

   B(0,1) = - A(0,1);
   B(0,2) =   tmp;
   B(1,2) = - A(1,2);
   B(0,3) =   A(0,1)*A(1,3) - A(0,3) - A(2,3)*tmp;
   B(1,3) =   A(2,3)*A(1,2) - A(1,3);
   B(2,3) = - A(2,3);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense \f$ 4 \times 4 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The diagonal dense matrix to be inverted.
// \return void
//
// This function inverts the given diagonal dense \f$ 4 \times 4 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonal4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 4UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 4UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET tmp1( A(2,2)*A(3,3) );
   const ET tmp2( A(0,0)*A(1,1) );
   const ET tmp3( A(0,0)*tmp1 );
   const ET tmp4( A(2,2)*tmp2 );

   const ET det( tmp1 * tmp2 );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );

   A(0,0) = A(1,1)*tmp1*idet;
   A(1,1) = tmp3*idet;
   A(2,2) = A(3,3)*tmp2*idet;
   A(3,3) = tmp4*idet;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense 4x4 matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense 4x4 matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert4x4( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertGeneral4x4  ( *dm ); break;
      case byLDLT     : invertSymmetric4x4( *dm ); break;
      case byLDLH     : invertHermitian4x4( *dm ); break;
      case byLLH      : invertHermitian4x4( *dm ); break;
      case asGeneral  : invertGeneral4x4  ( *dm ); break;
      case asSymmetric: invertSymmetric4x4( *dm ); break;
      case asHermitian: invertHermitian4x4( *dm ); break;
      case asLower    : invertLower4x4    ( *dm ); break;
      case asUniLower : invertUniLower4x4 ( *dm ); break;
      case asUpper    : invertUpper4x4    ( *dm ); break;
      case asUniUpper : invertUniUpper4x4 ( *dm ); break;
      case asDiagonal : invertDiagonal4x4 ( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR 5x5 MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertGeneral5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(3,3)*A(4,4) - A(3,4)*A(4,3) );
   ET tmp2 ( A(3,2)*A(4,4) - A(3,4)*A(4,2) );
   ET tmp3 ( A(3,2)*A(4,3) - A(3,3)*A(4,2) );
   ET tmp4 ( A(3,1)*A(4,4) - A(3,4)*A(4,1) );
   ET tmp5 ( A(3,1)*A(4,3) - A(3,3)*A(4,1) );
   ET tmp6 ( A(3,1)*A(4,2) - A(3,2)*A(4,1) );
   ET tmp7 ( A(3,0)*A(4,4) - A(3,4)*A(4,0) );
   ET tmp8 ( A(3,0)*A(4,3) - A(3,3)*A(4,0) );
   ET tmp9 ( A(3,0)*A(4,2) - A(3,2)*A(4,0) );
   ET tmp10( A(3,0)*A(4,1) - A(3,1)*A(4,0) );

   ET tmp11( A(2,2)*tmp1 - A(2,3)*tmp2 + A(2,4)*tmp3  );
   ET tmp12( A(2,1)*tmp1 - A(2,3)*tmp4 + A(2,4)*tmp5  );
   ET tmp13( A(2,1)*tmp2 - A(2,2)*tmp4 + A(2,4)*tmp6  );
   ET tmp14( A(2,1)*tmp3 - A(2,2)*tmp5 + A(2,3)*tmp6  );
   ET tmp15( A(2,0)*tmp1 - A(2,3)*tmp7 + A(2,4)*tmp8  );
   ET tmp16( A(2,0)*tmp2 - A(2,2)*tmp7 + A(2,4)*tmp9  );
   ET tmp17( A(2,0)*tmp3 - A(2,2)*tmp8 + A(2,3)*tmp9  );

   B(0,0) =   A(1,1)*tmp11 - A(1,2)*tmp12 + A(1,3)*tmp13 - A(1,4)*tmp14;
   B(0,1) = - A(0,1)*tmp11 + A(0,2)*tmp12 - A(0,3)*tmp13 + A(0,4)*tmp14;
   B(1,0) = - A(1,0)*tmp11 + A(1,2)*tmp15 - A(1,3)*tmp16 + A(1,4)*tmp17;
   B(1,1) =   A(0,0)*tmp11 - A(0,2)*tmp15 + A(0,3)*tmp16 - A(0,4)*tmp17;

   ET tmp18( A(2,0)*tmp4 - A(2,1)*tmp7 + A(2,4)*tmp10 );
   ET tmp19( A(2,0)*tmp5 - A(2,1)*tmp8 + A(2,3)*tmp10 );
   ET tmp20( A(2,0)*tmp6 - A(2,1)*tmp9 + A(2,2)*tmp10 );

   B(2,0) =   A(1,0)*tmp12 - A(1,1)*tmp15 + A(1,3)*tmp18 - A(1,4)*tmp19;
   B(2,1) = - A(0,0)*tmp12 + A(0,1)*tmp15 - A(0,3)*tmp18 + A(0,4)*tmp19;
   B(3,0) = - A(1,0)*tmp13 + A(1,1)*tmp16 - A(1,2)*tmp18 + A(1,4)*tmp20;
   B(3,1) =   A(0,0)*tmp13 - A(0,1)*tmp16 + A(0,2)*tmp18 - A(0,4)*tmp20;
   B(4,0) =   A(1,0)*tmp14 - A(1,1)*tmp17 + A(1,2)*tmp19 - A(1,3)*tmp20;
   B(4,1) = - A(0,0)*tmp14 + A(0,1)*tmp17 - A(0,2)*tmp19 + A(0,3)*tmp20;

   tmp11 = A(1,2)*tmp1 - A(1,3)*tmp2 + A(1,4)*tmp3;
   tmp12 = A(1,1)*tmp1 - A(1,3)*tmp4 + A(1,4)*tmp5;
   tmp13 = A(1,1)*tmp2 - A(1,2)*tmp4 + A(1,4)*tmp6;
   tmp14 = A(1,1)*tmp3 - A(1,2)*tmp5 + A(1,3)*tmp6;
   tmp15 = A(1,0)*tmp1 - A(1,3)*tmp7 + A(1,4)*tmp8;
   tmp16 = A(1,0)*tmp2 - A(1,2)*tmp7 + A(1,4)*tmp9;
   tmp17 = A(1,0)*tmp3 - A(1,2)*tmp8 + A(1,3)*tmp9;
   tmp18 = A(1,0)*tmp4 - A(1,1)*tmp7 + A(1,4)*tmp10;
   tmp19 = A(1,0)*tmp5 - A(1,1)*tmp8 + A(1,3)*tmp10;

   B(0,2) =   A(0,1)*tmp11 - A(0,2)*tmp12 + A(0,3)*tmp13 - A(0,4)*tmp14;
   B(1,2) = - A(0,0)*tmp11 + A(0,2)*tmp15 - A(0,3)*tmp16 + A(0,4)*tmp17;
   B(2,2) =   A(0,0)*tmp12 - A(0,1)*tmp15 + A(0,3)*tmp18 - A(0,4)*tmp19;

   tmp1  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp4  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp6  = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp7  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp8  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp9  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp10 = A(0,3)*A(1,4) - A(0,4)*A(1,3);

   tmp11 = A(2,2)*tmp10 - A(2,3)*tmp7 + A(2,4)*tmp1;
   tmp12 = A(2,1)*tmp10 - A(2,3)*tmp8 + A(2,4)*tmp2;
   tmp13 = A(2,1)*tmp7  - A(2,2)*tmp8 + A(2,4)*tmp3;
   tmp14 = A(2,1)*tmp1  - A(2,2)*tmp2 + A(2,3)*tmp3;
   tmp15 = A(2,0)*tmp10 - A(2,3)*tmp9 + A(2,4)*tmp4;
   tmp16 = A(2,0)*tmp7  - A(2,2)*tmp9 + A(2,4)*tmp5;
   tmp17 = A(2,0)*tmp1  - A(2,2)*tmp4 + A(2,3)*tmp5;

   B(0,3) =   A(4,1)*tmp11 - A(4,2)*tmp12 + A(4,3)*tmp13 - A(4,4)*tmp14;
   B(0,4) = - A(3,1)*tmp11 + A(3,2)*tmp12 - A(3,3)*tmp13 + A(3,4)*tmp14;
   B(1,3) = - A(4,0)*tmp11 + A(4,2)*tmp15 - A(4,3)*tmp16 + A(4,4)*tmp17;
   B(1,4) =   A(3,0)*tmp11 - A(3,2)*tmp15 + A(3,3)*tmp16 - A(3,4)*tmp17;

   tmp18 = A(2,0)*tmp8  - A(2,1)*tmp9 + A(2,4)*tmp6;
   tmp19 = A(2,0)*tmp2  - A(2,1)*tmp4 + A(2,3)*tmp6;
   tmp20 = A(2,0)*tmp3  - A(2,1)*tmp5 + A(2,2)*tmp6;

   B(2,3) =   A(4,0)*tmp12 - A(4,1)*tmp15 + A(4,3)*tmp18 - A(4,4)*tmp19;
   B(2,4) = - A(3,0)*tmp12 + A(3,1)*tmp15 - A(3,3)*tmp18 + A(3,4)*tmp19;
   B(3,3) = - A(4,0)*tmp13 + A(4,1)*tmp16 - A(4,2)*tmp18 + A(4,4)*tmp20;
   B(3,4) =   A(3,0)*tmp13 - A(3,1)*tmp16 + A(3,2)*tmp18 - A(3,4)*tmp20;
   B(4,3) =   A(4,0)*tmp14 - A(4,1)*tmp17 + A(4,2)*tmp19 - A(4,3)*tmp20;
   B(4,4) = - A(3,0)*tmp14 + A(3,1)*tmp17 - A(3,2)*tmp19 + A(3,3)*tmp20;

   tmp11 = A(3,1)*tmp7  - A(3,2)*tmp8 + A(3,4)*tmp3;
   tmp12 = A(3,0)*tmp7  - A(3,2)*tmp9 + A(3,4)*tmp5;
   tmp13 = A(3,0)*tmp8  - A(3,1)*tmp9 + A(3,4)*tmp6;
   tmp14 = A(3,0)*tmp3  - A(3,1)*tmp5 + A(3,2)*tmp6;

   tmp15 = A(3,1)*tmp1  - A(3,2)*tmp2 + A(3,3)*tmp3;
   tmp16 = A(3,0)*tmp1  - A(3,2)*tmp4 + A(3,3)*tmp5;
   tmp17 = A(3,0)*tmp2  - A(3,1)*tmp4 + A(3,3)*tmp6;

   B(3,2) =   A(4,0)*tmp11 - A(4,1)*tmp12 + A(4,2)*tmp13 - A(4,4)*tmp14;
   B(4,2) = - A(4,0)*tmp15 + A(4,1)*tmp16 - A(4,2)*tmp17 + A(4,3)*tmp14;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) + A(0,4)*B(4,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric/Hermitian dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The symmetric/Hermitian dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense \f$ 5 \times 5 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertSymmetric5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(3,3)*A(4,4) - A(3,4)*A(4,3) );
   ET tmp2 ( A(3,2)*A(4,4) - A(3,4)*A(4,2) );
   ET tmp3 ( A(3,2)*A(4,3) - A(3,3)*A(4,2) );
   ET tmp4 ( A(3,1)*A(4,4) - A(3,4)*A(4,1) );
   ET tmp5 ( A(3,1)*A(4,3) - A(3,3)*A(4,1) );
   ET tmp6 ( A(3,1)*A(4,2) - A(3,2)*A(4,1) );
   ET tmp7 ( A(3,0)*A(4,4) - A(3,4)*A(4,0) );
   ET tmp8 ( A(3,0)*A(4,3) - A(3,3)*A(4,0) );
   ET tmp9 ( A(3,0)*A(4,2) - A(3,2)*A(4,0) );
   ET tmp10( A(3,0)*A(4,1) - A(3,1)*A(4,0) );

   ET tmp11( A(2,2)*tmp1 - A(2,3)*tmp2 + A(2,4)*tmp3  );
   ET tmp12( A(2,1)*tmp1 - A(2,3)*tmp4 + A(2,4)*tmp5  );
   ET tmp13( A(2,1)*tmp2 - A(2,2)*tmp4 + A(2,4)*tmp6  );
   ET tmp14( A(2,1)*tmp3 - A(2,2)*tmp5 + A(2,3)*tmp6  );
   ET tmp15( A(2,0)*tmp1 - A(2,3)*tmp7 + A(2,4)*tmp8  );
   ET tmp16( A(2,0)*tmp2 - A(2,2)*tmp7 + A(2,4)*tmp9  );
   ET tmp17( A(2,0)*tmp3 - A(2,2)*tmp8 + A(2,3)*tmp9  );

   B(0,0) =   A(1,1)*tmp11 - A(1,2)*tmp12 + A(1,3)*tmp13 - A(1,4)*tmp14;
   B(0,1) = - A(0,1)*tmp11 + A(0,2)*tmp12 - A(0,3)*tmp13 + A(0,4)*tmp14;
   B(1,1) =   A(0,0)*tmp11 - A(0,2)*tmp15 + A(0,3)*tmp16 - A(0,4)*tmp17;

   ET tmp18( A(2,0)*tmp4 - A(2,1)*tmp7 + A(2,4)*tmp10 );
   ET tmp19( A(2,0)*tmp5 - A(2,1)*tmp8 + A(2,3)*tmp10 );
   ET tmp20( A(2,0)*tmp6 - A(2,1)*tmp9 + A(2,2)*tmp10 );

   B(2,0) =   A(1,0)*tmp12 - A(1,1)*tmp15 + A(1,3)*tmp18 - A(1,4)*tmp19;
   B(2,1) = - A(0,0)*tmp12 + A(0,1)*tmp15 - A(0,3)*tmp18 + A(0,4)*tmp19;
   B(3,0) = - A(1,0)*tmp13 + A(1,1)*tmp16 - A(1,2)*tmp18 + A(1,4)*tmp20;
   B(3,1) =   A(0,0)*tmp13 - A(0,1)*tmp16 + A(0,2)*tmp18 - A(0,4)*tmp20;
   B(4,0) =   A(1,0)*tmp14 - A(1,1)*tmp17 + A(1,2)*tmp19 - A(1,3)*tmp20;
   B(4,1) = - A(0,0)*tmp14 + A(0,1)*tmp17 - A(0,2)*tmp19 + A(0,3)*tmp20;

   tmp11 = A(1,1)*tmp1 - A(1,3)*tmp4 + A(1,4)*tmp5;
   tmp12 = A(1,0)*tmp1 - A(1,3)*tmp7 + A(1,4)*tmp8;
   tmp13 = A(1,0)*tmp4 - A(1,1)*tmp7 + A(1,4)*tmp10;
   tmp14 = A(1,0)*tmp5 - A(1,1)*tmp8 + A(1,3)*tmp10;

   B(2,2) = A(0,0)*tmp11 - A(0,1)*tmp12 + A(0,3)*tmp13 - A(0,4)*tmp14;

   tmp1  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp4  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp6  = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp7  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp8  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp9  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp10 = A(0,3)*A(1,4) - A(0,4)*A(1,3);

   tmp11 = A(2,1)*tmp10 - A(2,3)*tmp8 + A(2,4)*tmp2;
   tmp12 = A(2,1)*tmp7  - A(2,2)*tmp8 + A(2,4)*tmp3;
   tmp13 = A(2,1)*tmp1  - A(2,2)*tmp2 + A(2,3)*tmp3;
   tmp14 = A(2,0)*tmp10 - A(2,3)*tmp9 + A(2,4)*tmp4;
   tmp15 = A(2,0)*tmp7  - A(2,2)*tmp9 + A(2,4)*tmp5;
   tmp16 = A(2,0)*tmp1  - A(2,2)*tmp4 + A(2,3)*tmp5;
   tmp17 = A(2,0)*tmp8  - A(2,1)*tmp9 + A(2,4)*tmp6;
   tmp18 = A(2,0)*tmp2  - A(2,1)*tmp4 + A(2,3)*tmp6;
   tmp19 = A(2,0)*tmp3  - A(2,1)*tmp5 + A(2,2)*tmp6;

   B(2,3) =   A(4,0)*tmp11 - A(4,1)*tmp14 + A(4,3)*tmp17 - A(4,4)*tmp18;
   B(2,4) = - A(3,0)*tmp11 + A(3,1)*tmp14 - A(3,3)*tmp17 + A(3,4)*tmp18;
   B(3,3) = - A(4,0)*tmp12 + A(4,1)*tmp15 - A(4,2)*tmp17 + A(4,4)*tmp19;
   B(3,4) =   A(3,0)*tmp12 - A(3,1)*tmp15 + A(3,2)*tmp17 - A(3,4)*tmp19;
   B(4,4) = - A(3,0)*tmp13 + A(3,1)*tmp16 - A(3,2)*tmp18 + A(3,3)*tmp19;

   B(0,2) = B(2,0);
   B(0,3) = B(3,0);
   B(0,4) = B(4,0);
   B(1,0) = B(0,1);
   B(1,2) = B(2,1);
   B(1,3) = B(3,1);
   B(1,4) = B(4,1);
   B(3,2) = B(2,3);
   B(4,2) = B(2,4);
   B(4,3) = B(3,4);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) + A(0,4)*B(4,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given Hermitian dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The Hermitian dense matrix to be inverted.
// \return void
//
// This function inverts the given Hermitian dense \f$ 5 \times 5 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertHermitian5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(3,3)*A(4,4) - A(3,4)*A(4,3) );
   ET tmp2 ( A(3,2)*A(4,4) - A(3,4)*A(4,2) );
   ET tmp3 ( A(3,2)*A(4,3) - A(3,3)*A(4,2) );
   ET tmp4 ( A(3,1)*A(4,4) - A(3,4)*A(4,1) );
   ET tmp5 ( A(3,1)*A(4,3) - A(3,3)*A(4,1) );
   ET tmp6 ( A(3,1)*A(4,2) - A(3,2)*A(4,1) );
   ET tmp7 ( A(3,0)*A(4,4) - A(3,4)*A(4,0) );
   ET tmp8 ( A(3,0)*A(4,3) - A(3,3)*A(4,0) );
   ET tmp9 ( A(3,0)*A(4,2) - A(3,2)*A(4,0) );
   ET tmp10( A(3,0)*A(4,1) - A(3,1)*A(4,0) );

   ET tmp11( A(2,2)*tmp1 - A(2,3)*tmp2 + A(2,4)*tmp3  );
   ET tmp12( A(2,1)*tmp1 - A(2,3)*tmp4 + A(2,4)*tmp5  );
   ET tmp13( A(2,1)*tmp2 - A(2,2)*tmp4 + A(2,4)*tmp6  );
   ET tmp14( A(2,1)*tmp3 - A(2,2)*tmp5 + A(2,3)*tmp6  );
   ET tmp15( A(2,0)*tmp1 - A(2,3)*tmp7 + A(2,4)*tmp8  );
   ET tmp16( A(2,0)*tmp2 - A(2,2)*tmp7 + A(2,4)*tmp9  );
   ET tmp17( A(2,0)*tmp3 - A(2,2)*tmp8 + A(2,3)*tmp9  );

   B(0,0) = ET( real( A(1,1)*tmp11 - A(1,2)*tmp12 + A(1,3)*tmp13 - A(1,4)*tmp14 ) );
   B(0,1) = - A(0,1)*tmp11 + A(0,2)*tmp12 - A(0,3)*tmp13 + A(0,4)*tmp14;
   B(1,1) = ET( real( A(0,0)*tmp11 - A(0,2)*tmp15 + A(0,3)*tmp16 - A(0,4)*tmp17 ) );

   ET tmp18( A(2,0)*tmp4 - A(2,1)*tmp7 + A(2,4)*tmp10 );
   ET tmp19( A(2,0)*tmp5 - A(2,1)*tmp8 + A(2,3)*tmp10 );
   ET tmp20( A(2,0)*tmp6 - A(2,1)*tmp9 + A(2,2)*tmp10 );

   B(2,0) =   A(1,0)*tmp12 - A(1,1)*tmp15 + A(1,3)*tmp18 - A(1,4)*tmp19;
   B(2,1) = - A(0,0)*tmp12 + A(0,1)*tmp15 - A(0,3)*tmp18 + A(0,4)*tmp19;
   B(3,0) = - A(1,0)*tmp13 + A(1,1)*tmp16 - A(1,2)*tmp18 + A(1,4)*tmp20;
   B(3,1) =   A(0,0)*tmp13 - A(0,1)*tmp16 + A(0,2)*tmp18 - A(0,4)*tmp20;
   B(4,0) =   A(1,0)*tmp14 - A(1,1)*tmp17 + A(1,2)*tmp19 - A(1,3)*tmp20;
   B(4,1) = - A(0,0)*tmp14 + A(0,1)*tmp17 - A(0,2)*tmp19 + A(0,3)*tmp20;

   tmp11 = A(1,1)*tmp1 - A(1,3)*tmp4 + A(1,4)*tmp5;
   tmp12 = A(1,0)*tmp1 - A(1,3)*tmp7 + A(1,4)*tmp8;
   tmp13 = A(1,0)*tmp4 - A(1,1)*tmp7 + A(1,4)*tmp10;
   tmp14 = A(1,0)*tmp5 - A(1,1)*tmp8 + A(1,3)*tmp10;

   B(2,2) = ET( real( A(0,0)*tmp11 - A(0,1)*tmp12 + A(0,3)*tmp13 - A(0,4)*tmp14 ) );

   tmp1  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp2  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp3  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp4  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp5  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp6  = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp7  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp8  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp9  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp10 = A(0,3)*A(1,4) - A(0,4)*A(1,3);

   tmp11 = A(2,1)*tmp10 - A(2,3)*tmp8 + A(2,4)*tmp2;
   tmp12 = A(2,1)*tmp7  - A(2,2)*tmp8 + A(2,4)*tmp3;
   tmp13 = A(2,1)*tmp1  - A(2,2)*tmp2 + A(2,3)*tmp3;
   tmp14 = A(2,0)*tmp10 - A(2,3)*tmp9 + A(2,4)*tmp4;
   tmp15 = A(2,0)*tmp7  - A(2,2)*tmp9 + A(2,4)*tmp5;
   tmp16 = A(2,0)*tmp1  - A(2,2)*tmp4 + A(2,3)*tmp5;
   tmp17 = A(2,0)*tmp8  - A(2,1)*tmp9 + A(2,4)*tmp6;
   tmp18 = A(2,0)*tmp2  - A(2,1)*tmp4 + A(2,3)*tmp6;
   tmp19 = A(2,0)*tmp3  - A(2,1)*tmp5 + A(2,2)*tmp6;

   B(2,3) =   A(4,0)*tmp11 - A(4,1)*tmp14 + A(4,3)*tmp17 - A(4,4)*tmp18;
   B(2,4) = - A(3,0)*tmp11 + A(3,1)*tmp14 - A(3,3)*tmp17 + A(3,4)*tmp18;
   B(3,3) = - ET( real( A(4,0)*tmp12 - A(4,1)*tmp15 + A(4,2)*tmp17 - A(4,4)*tmp19 ) );
   B(3,4) =   A(3,0)*tmp12 - A(3,1)*tmp15 + A(3,2)*tmp17 - A(3,4)*tmp19;
   B(4,4) = - ET( real( A(3,0)*tmp13 - A(3,1)*tmp16 + A(3,2)*tmp18 - A(3,3)*tmp19 ) );

   B(0,2) = conj( B(2,0) );
   B(0,3) = conj( B(3,0) );
   B(0,4) = conj( B(4,0) );
   B(1,0) = conj( B(0,1) );
   B(1,2) = conj( B(2,1) );
   B(1,3) = conj( B(3,1) );
   B(1,4) = conj( B(4,1) );
   B(3,2) = conj( B(2,3) );
   B(4,2) = conj( B(2,4) );
   B(4,3) = conj( B(3,4) );

   const ET det( real( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) + A(0,3)*B(3,0) + A(0,4)*B(4,0) ) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The lower dense matrix to be inverted.
// \return void
//
// This function inverts the given lower dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLower5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(3,3)*A(4,4) );
   const ET tmp2( A(3,2)*A(4,4) );
   const ET tmp3( A(3,2)*A(4,3) - A(3,3)*A(4,2) );
   const ET tmp4( A(0,0)*A(1,1) );

   const ET tmp5 ( A(2,2)*tmp1 );
   const ET tmp6 ( A(2,1)*tmp1 );
   const ET tmp7 ( A(2,1)*tmp2 - A(2,2)*A(3,1)*A(4,4) );
   const ET tmp8 ( A(2,1)*tmp3 - A(2,2)*( A(3,1)*A(4,3) - A(3,3)*A(4,1) ) );
   const ET tmp9 ( A(3,2)*tmp4 );
   const ET tmp10( A(2,2)*tmp4 );

   B(0,0) =   A(1,1)*tmp5;
   B(1,0) = - A(1,0)*tmp5;
   B(2,0) =   A(1,0)*tmp6 - A(1,1)*A(2,0)*tmp1;
   B(3,0) =   A(1,1)*( A(2,0)*tmp2 - A(2,2)*A(3,0)*A(4,4) ) - A(1,0)*tmp7;
   B(4,0) =   A(1,0)*tmp8 - A(1,1)*( A(2,0)*tmp3 - A(2,2)*( A(3,0)*A(4,3) - A(3,3)*A(4,0) ) );
   B(1,1) =   A(0,0)*tmp5;
   B(2,1) = - A(0,0)*tmp6;
   B(3,1) =   A(0,0)*tmp7;
   B(4,1) = - A(0,0)*tmp8;
   B(2,2) =   A(0,0)*A(1,1)*tmp1;
   B(3,2) = - A(4,4)*tmp9;
   B(4,2) =   A(4,3)*tmp9 - A(4,2)*A(3,3)*tmp4;
   B(3,3) =   A(4,4)*tmp10;
   B(4,3) = - A(4,3)*tmp10;
   B(4,4) =   A(3,3)*tmp10;

   const ET det( B(4,4) * A(4,4) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given unilower dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The unilower dense matrix to be inverted.
// \return void
//
// This function inverts the given unilower dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLower5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(3,2)*A(4,3) - A(4,2) );
   const ET tmp2( A(2,1)*A(3,2) - A(3,1) );
   const ET tmp3( A(2,1)*tmp1 - A(3,1)*A(4,3) + A(4,1) );

   B(1,0) = - A(1,0);
   B(2,0) =   A(1,0)*A(2,1) - A(2,0);
   B(3,0) = - A(1,0)*tmp2 + A(2,0)*A(3,2) - A(3,0);
   B(4,0) =   A(1,0)*tmp3 - A(2,0)*tmp1 + A(3,0)*A(4,3) - A(4,0);
   B(2,1) = - A(2,1);
   B(3,1) =   tmp2;
   B(4,1) = - tmp3;
   B(3,2) = - A(3,2);
   B(4,2) =   A(4,3)*A(3,2) - A(4,2);
   B(4,3) = - A(4,3);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The upper dense matrix to be inverted.
// \return void
//
// This function inverts the given upper dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpper5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(3,3)*A(4,4) );
   const ET tmp2( A(0,1)*A(1,2) - A(0,2)*A(1,1) );
   const ET tmp3( A(0,0)*A(1,2) );
   const ET tmp4( A(0,0)*A(1,1) );

   const ET tmp5 ( A(2,2)*tmp1 );
   const ET tmp6 ( A(1,2)*tmp1 );
   const ET tmp7 ( A(1,1)*tmp1 );
   const ET tmp8 ( A(2,3)*tmp2 - A(2,2)*( A(0,1)*A(1,3) - A(0,3)*A(1,1) ) );
   const ET tmp9 ( A(2,3)*tmp3 - A(2,2)*A(0,0)*A(1,3) );
   const ET tmp10( A(2,3)*tmp4 );
   const ET tmp11( A(2,2)*tmp4 );

   B(0,0) =   A(1,1)*tmp5;
   B(0,1) = - A(0,1)*tmp5;
   B(1,1) =   A(0,0)*tmp5;
   B(0,2) =   A(0,1)*tmp6 - A(0,2)*tmp7;
   B(1,2) = - A(0,0)*tmp6;
   B(2,2) =   A(0,0)*tmp7;
   B(0,3) = - A(4,4)*tmp8;
   B(1,3) =   A(4,4)*tmp9;
   B(2,3) = - A(4,4)*tmp10;
   B(3,3) =   A(4,4)*tmp11;
   B(0,4) =   A(3,4)*tmp8 - A(3,3)*( A(2,4)*tmp2 - A(2,2)*( A(0,1)*A(1,4) - A(0,4)*A(1,1) ) );
   B(1,4) =   A(3,3)*( A(2,4)*tmp3 - A(2,2)*A(0,0)*A(1,4) ) - A(3,4)*tmp9;
   B(2,4) =   A(3,4)*tmp10 - A(3,3)*A(2,4)*tmp4;
   B(3,4) = - A(3,4)*tmp11;
   B(4,4) =   A(3,3)*tmp11;

   const ET det( A(0,0) * B(0,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The uniupper dense matrix to be inverted.
// \return void
//
// This function inverts the given uniupper dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpper5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,5UL,5UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(0,1)*A(1,2) - A(0,2) );
   const ET tmp2( A(2,3)*A(1,2) - A(1,3) );
   const ET tmp3( A(2,3)*tmp1 - A(0,1)*A(1,3) + A(0,3) );

   B(0,1) = - A(0,1);
   B(0,2) =   A(0,1)*A(1,2) - A(0,2);
   B(1,2) = - A(1,2);
   B(0,3) = - tmp3;
   B(1,3) =   tmp2;
   B(2,3) = - A(2,3);
   B(0,4) =   A(3,4)*tmp3 - A(2,4)*tmp1 + A(0,1)*A(1,4) - A(0,4);
   B(1,4) =   A(2,4)*A(1,2) - A(1,4) - A(3,4)*tmp2;
   B(2,4) =   A(3,4)*A(2,3) - A(2,4);
   B(3,4) = - A(3,4);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense \f$ 5 \times 5 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The diagonal dense matrix to be inverted.
// \return void
//
// This function inverts the given diagonal dense \f$ 5 \times 5 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonal5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 5UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 5UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET tmp1( A(0,0)*A(1,1) );
   const ET tmp2( A(3,3)*A(4,4) );
   const ET tmp3( A(0,0)*tmp2 );
   const ET tmp4( tmp1*A(2,2) );
   const ET tmp5( tmp4*A(3,3) );

   const ET det( tmp2*tmp4 );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );

   A(0,0) = A(1,1)*A(2,2)*tmp2*idet;
   A(1,1) = A(2,2)*tmp3*idet;
   A(2,2) = tmp1*tmp2*idet;
   A(3,3) = tmp4*A(4,4)*idet;
   A(4,4) = tmp5*idet;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense 5x5 matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense 5x5 matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert5x5( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertGeneral5x5  ( *dm ); break;
      case byLDLT     : invertSymmetric5x5( *dm ); break;
      case byLDLH     : invertHermitian5x5( *dm ); break;
      case byLLH      : invertHermitian5x5( *dm ); break;
      case asGeneral  : invertGeneral5x5  ( *dm ); break;
      case asSymmetric: invertSymmetric5x5( *dm ); break;
      case asHermitian: invertHermitian5x5( *dm ); break;
      case asLower    : invertLower5x5    ( *dm ); break;
      case asUniLower : invertUniLower5x5 ( *dm ); break;
      case asUpper    : invertUpper5x5    ( *dm ); break;
      case asUniUpper : invertUniUpper5x5 ( *dm ); break;
      case asDiagonal : invertDiagonal5x5 ( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR 6x6 MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given general dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The general dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertGeneral6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(4,4)*A(5,5) - A(4,5)*A(5,4) );
   ET tmp2 ( A(4,3)*A(5,5) - A(4,5)*A(5,3) );
   ET tmp3 ( A(4,3)*A(5,4) - A(4,4)*A(5,3) );
   ET tmp4 ( A(4,2)*A(5,5) - A(4,5)*A(5,2) );
   ET tmp5 ( A(4,2)*A(5,4) - A(4,4)*A(5,2) );
   ET tmp6 ( A(4,2)*A(5,3) - A(4,3)*A(5,2) );
   ET tmp7 ( A(4,1)*A(5,5) - A(4,5)*A(5,1) );
   ET tmp8 ( A(4,1)*A(5,4) - A(4,4)*A(5,1) );
   ET tmp9 ( A(4,1)*A(5,3) - A(4,3)*A(5,1) );
   ET tmp10( A(4,1)*A(5,2) - A(4,2)*A(5,1) );
   ET tmp11( A(4,0)*A(5,5) - A(4,5)*A(5,0) );
   ET tmp12( A(4,0)*A(5,4) - A(4,4)*A(5,0) );
   ET tmp13( A(4,0)*A(5,3) - A(4,3)*A(5,0) );
   ET tmp14( A(4,0)*A(5,2) - A(4,2)*A(5,0) );
   ET tmp15( A(4,0)*A(5,1) - A(4,1)*A(5,0) );

   ET tmp16( A(3,3)*tmp1  - A(3,4)*tmp2  + A(3,5)*tmp3  );
   ET tmp17( A(3,2)*tmp1  - A(3,4)*tmp4  + A(3,5)*tmp5  );
   ET tmp18( A(3,2)*tmp2  - A(3,3)*tmp4  + A(3,5)*tmp6  );
   ET tmp19( A(3,2)*tmp3  - A(3,3)*tmp5  + A(3,4)*tmp6  );
   ET tmp20( A(3,1)*tmp1  - A(3,4)*tmp7  + A(3,5)*tmp8  );
   ET tmp21( A(3,1)*tmp2  - A(3,3)*tmp7  + A(3,5)*tmp9  );
   ET tmp22( A(3,1)*tmp3  - A(3,3)*tmp8  + A(3,4)*tmp9  );
   ET tmp23( A(3,1)*tmp4  - A(3,2)*tmp7  + A(3,5)*tmp10 );
   ET tmp24( A(3,1)*tmp5  - A(3,2)*tmp8  + A(3,4)*tmp10 );
   ET tmp25( A(3,1)*tmp6  - A(3,2)*tmp9  + A(3,3)*tmp10 );
   ET tmp26( A(3,0)*tmp1  - A(3,4)*tmp11 + A(3,5)*tmp12 );
   ET tmp27( A(3,0)*tmp2  - A(3,3)*tmp11 + A(3,5)*tmp13 );
   ET tmp28( A(3,0)*tmp3  - A(3,3)*tmp12 + A(3,4)*tmp13 );
   ET tmp29( A(3,0)*tmp4  - A(3,2)*tmp11 + A(3,5)*tmp14 );
   ET tmp30( A(3,0)*tmp5  - A(3,2)*tmp12 + A(3,4)*tmp14 );
   ET tmp31( A(3,0)*tmp6  - A(3,2)*tmp13 + A(3,3)*tmp14 );
   ET tmp32( A(3,0)*tmp7  - A(3,1)*tmp11 + A(3,5)*tmp15 );
   ET tmp33( A(3,0)*tmp8  - A(3,1)*tmp12 + A(3,4)*tmp15 );
   ET tmp34( A(3,0)*tmp9  - A(3,1)*tmp13 + A(3,3)*tmp15 );
   ET tmp35( A(3,0)*tmp10 - A(3,1)*tmp14 + A(3,2)*tmp15 );

   ET tmp36( A(2,2)*tmp16 - A(2,3)*tmp17 + A(2,4)*tmp18 - A(2,5)*tmp19 );
   ET tmp37( A(2,1)*tmp16 - A(2,3)*tmp20 + A(2,4)*tmp21 - A(2,5)*tmp22 );
   ET tmp38( A(2,1)*tmp17 - A(2,2)*tmp20 + A(2,4)*tmp23 - A(2,5)*tmp24 );
   ET tmp39( A(2,1)*tmp18 - A(2,2)*tmp21 + A(2,3)*tmp23 - A(2,5)*tmp25 );
   ET tmp40( A(2,1)*tmp19 - A(2,2)*tmp22 + A(2,3)*tmp24 - A(2,4)*tmp25 );
   ET tmp41( A(2,0)*tmp16 - A(2,3)*tmp26 + A(2,4)*tmp27 - A(2,5)*tmp28 );
   ET tmp42( A(2,0)*tmp17 - A(2,2)*tmp26 + A(2,4)*tmp29 - A(2,5)*tmp30 );
   ET tmp43( A(2,0)*tmp18 - A(2,2)*tmp27 + A(2,3)*tmp29 - A(2,5)*tmp31 );
   ET tmp44( A(2,0)*tmp19 - A(2,2)*tmp28 + A(2,3)*tmp30 - A(2,4)*tmp31 );

   B(0,0) =   A(1,1)*tmp36 - A(1,2)*tmp37 + A(1,3)*tmp38 - A(1,4)*tmp39 + A(1,5)*tmp40;
   B(0,1) = - A(0,1)*tmp36 + A(0,2)*tmp37 - A(0,3)*tmp38 + A(0,4)*tmp39 - A(0,5)*tmp40;
   B(1,0) = - A(1,0)*tmp36 + A(1,2)*tmp41 - A(1,3)*tmp42 + A(1,4)*tmp43 - A(1,5)*tmp44;
   B(1,1) =   A(0,0)*tmp36 - A(0,2)*tmp41 + A(0,3)*tmp42 - A(0,4)*tmp43 + A(0,5)*tmp44;

   ET tmp45( A(2,0)*tmp20 - A(2,1)*tmp26 + A(2,4)*tmp32 - A(2,5)*tmp33 );
   ET tmp46( A(2,0)*tmp21 - A(2,1)*tmp27 + A(2,3)*tmp32 - A(2,5)*tmp34 );
   ET tmp47( A(2,0)*tmp22 - A(2,1)*tmp28 + A(2,3)*tmp33 - A(2,4)*tmp34 );
   ET tmp48( A(2,0)*tmp23 - A(2,1)*tmp29 + A(2,2)*tmp32 - A(2,5)*tmp35 );
   ET tmp49( A(2,0)*tmp24 - A(2,1)*tmp30 + A(2,2)*tmp33 - A(2,4)*tmp35 );

   B(2,0) =   A(1,0)*tmp37 - A(1,1)*tmp41 + A(1,3)*tmp45 - A(1,4)*tmp46 + A(1,5)*tmp47;
   B(2,1) = - A(0,0)*tmp37 + A(0,1)*tmp41 - A(0,3)*tmp45 + A(0,4)*tmp46 - A(0,5)*tmp47;
   B(3,0) = - A(1,0)*tmp38 + A(1,1)*tmp42 - A(1,2)*tmp45 + A(1,4)*tmp48 - A(1,5)*tmp49;
   B(3,1) =   A(0,0)*tmp38 - A(0,1)*tmp42 + A(0,2)*tmp45 - A(0,4)*tmp48 + A(0,5)*tmp49;

   ET tmp50( A(2,0)*tmp25 - A(2,1)*tmp31 + A(2,2)*tmp34 - A(2,3)*tmp35 );

   B(4,0) =   A(1,0)*tmp39 - A(1,1)*tmp43 + A(1,2)*tmp46 - A(1,3)*tmp48 + A(1,5)*tmp50;
   B(4,1) = - A(0,0)*tmp39 + A(0,1)*tmp43 - A(0,2)*tmp46 + A(0,3)*tmp48 - A(0,5)*tmp50;
   B(5,0) = - A(1,0)*tmp40 + A(1,1)*tmp44 - A(1,2)*tmp47 + A(1,3)*tmp49 - A(1,4)*tmp50;
   B(5,1) =   A(0,0)*tmp40 - A(0,1)*tmp44 + A(0,2)*tmp47 - A(0,3)*tmp49 + A(0,4)*tmp50;

   tmp36 = A(1,2)*tmp16 - A(1,3)*tmp17 + A(1,4)*tmp18 - A(1,5)*tmp19;
   tmp37 = A(1,1)*tmp16 - A(1,3)*tmp20 + A(1,4)*tmp21 - A(1,5)*tmp22;
   tmp38 = A(1,1)*tmp17 - A(1,2)*tmp20 + A(1,4)*tmp23 - A(1,5)*tmp24;
   tmp39 = A(1,1)*tmp18 - A(1,2)*tmp21 + A(1,3)*tmp23 - A(1,5)*tmp25;
   tmp40 = A(1,1)*tmp19 - A(1,2)*tmp22 + A(1,3)*tmp24 - A(1,4)*tmp25;
   tmp41 = A(1,0)*tmp16 - A(1,3)*tmp26 + A(1,4)*tmp27 - A(1,5)*tmp28;
   tmp42 = A(1,0)*tmp17 - A(1,2)*tmp26 + A(1,4)*tmp29 - A(1,5)*tmp30;
   tmp43 = A(1,0)*tmp18 - A(1,2)*tmp27 + A(1,3)*tmp29 - A(1,5)*tmp31;
   tmp44 = A(1,0)*tmp19 - A(1,2)*tmp28 + A(1,3)*tmp30 - A(1,4)*tmp31;
   tmp45 = A(1,0)*tmp20 - A(1,1)*tmp26 + A(1,4)*tmp32 - A(1,5)*tmp33;
   tmp46 = A(1,0)*tmp21 - A(1,1)*tmp27 + A(1,3)*tmp32 - A(1,5)*tmp34;
   tmp47 = A(1,0)*tmp22 - A(1,1)*tmp28 + A(1,3)*tmp33 - A(1,4)*tmp34;
   tmp48 = A(1,0)*tmp23 - A(1,1)*tmp29 + A(1,2)*tmp32 - A(1,5)*tmp35;
   tmp49 = A(1,0)*tmp24 - A(1,1)*tmp30 + A(1,2)*tmp33 - A(1,4)*tmp35;
   tmp50 = A(1,0)*tmp25 - A(1,1)*tmp31 + A(1,2)*tmp34 - A(1,3)*tmp35;

   B(0,2) =   A(0,1)*tmp36 - A(0,2)*tmp37 + A(0,3)*tmp38 - A(0,4)*tmp39 + A(0,5)*tmp40;
   B(1,2) = - A(0,0)*tmp36 + A(0,2)*tmp41 - A(0,3)*tmp42 + A(0,4)*tmp43 - A(0,5)*tmp44;
   B(2,2) =   A(0,0)*tmp37 - A(0,1)*tmp41 + A(0,3)*tmp45 - A(0,4)*tmp46 + A(0,5)*tmp47;
   B(3,2) = - A(0,0)*tmp38 + A(0,1)*tmp42 - A(0,2)*tmp45 + A(0,4)*tmp48 - A(0,5)*tmp49;
   B(4,2) =   A(0,0)*tmp39 - A(0,1)*tmp43 + A(0,2)*tmp46 - A(0,3)*tmp48 + A(0,5)*tmp50;
   B(5,2) = - A(0,0)*tmp40 + A(0,1)*tmp44 - A(0,2)*tmp47 + A(0,3)*tmp49 - A(0,4)*tmp50;

   tmp1  = A(0,3)*A(1,4) - A(0,4)*A(1,3);
   tmp2  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp3  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp4  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp5  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp6  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp7  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp8  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp9  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp10 = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp11 = A(0,3)*A(1,5) - A(0,5)*A(1,3);
   tmp12 = A(0,2)*A(1,5) - A(0,5)*A(1,2);
   tmp13 = A(0,1)*A(1,5) - A(0,5)*A(1,1);
   tmp14 = A(0,0)*A(1,5) - A(0,5)*A(1,0);
   tmp15 = A(0,4)*A(1,5) - A(0,5)*A(1,4);

   tmp16 = A(2,3)*tmp15 - A(2,4)*tmp11 + A(2,5)*tmp1;
   tmp17 = A(2,2)*tmp15 - A(2,4)*tmp12 + A(2,5)*tmp2;
   tmp18 = A(2,2)*tmp11 - A(2,3)*tmp12 + A(2,5)*tmp3;
   tmp19 = A(2,2)*tmp1  - A(2,3)*tmp2  + A(2,4)*tmp3;
   tmp20 = A(2,1)*tmp15 - A(2,4)*tmp13 + A(2,5)*tmp4;
   tmp21 = A(2,1)*tmp11 - A(2,3)*tmp13 + A(2,5)*tmp5;
   tmp22 = A(2,1)*tmp1  - A(2,3)*tmp4  + A(2,4)*tmp5;
   tmp23 = A(2,1)*tmp12 - A(2,2)*tmp13 + A(2,5)*tmp6;
   tmp24 = A(2,1)*tmp2  - A(2,2)*tmp4  + A(2,4)*tmp6;
   tmp25 = A(2,1)*tmp3  - A(2,2)*tmp5  + A(2,3)*tmp6;
   tmp26 = A(2,0)*tmp15 - A(2,4)*tmp14 + A(2,5)*tmp7;
   tmp27 = A(2,0)*tmp11 - A(2,3)*tmp14 + A(2,5)*tmp8;
   tmp28 = A(2,0)*tmp1  - A(2,3)*tmp7  + A(2,4)*tmp8;
   tmp29 = A(2,0)*tmp12 - A(2,2)*tmp14 + A(2,5)*tmp9;
   tmp30 = A(2,0)*tmp2  - A(2,2)*tmp7  + A(2,4)*tmp9;
   tmp31 = A(2,0)*tmp3  - A(2,2)*tmp8  + A(2,3)*tmp9;
   tmp32 = A(2,0)*tmp13 - A(2,1)*tmp14 + A(2,5)*tmp10;
   tmp33 = A(2,0)*tmp4  - A(2,1)*tmp7  + A(2,4)*tmp10;
   tmp34 = A(2,0)*tmp5  - A(2,1)*tmp8  + A(2,3)*tmp10;
   tmp35 = A(2,0)*tmp6  - A(2,1)*tmp9  + A(2,2)*tmp10;

   tmp36 = A(3,2)*tmp16 - A(3,3)*tmp17 + A(3,4)*tmp18 - A(3,5)*tmp19;
   tmp37 = A(3,1)*tmp16 - A(3,3)*tmp20 + A(3,4)*tmp21 - A(3,5)*tmp22;
   tmp38 = A(3,1)*tmp17 - A(3,2)*tmp20 + A(3,4)*tmp23 - A(3,5)*tmp24;
   tmp39 = A(3,1)*tmp18 - A(3,2)*tmp21 + A(3,3)*tmp23 - A(3,5)*tmp25;
   tmp40 = A(3,1)*tmp19 - A(3,2)*tmp22 + A(3,3)*tmp24 - A(3,4)*tmp25;
   tmp41 = A(3,0)*tmp16 - A(3,3)*tmp26 + A(3,4)*tmp27 - A(3,5)*tmp28;
   tmp42 = A(3,0)*tmp17 - A(3,2)*tmp26 + A(3,4)*tmp29 - A(3,5)*tmp30;
   tmp43 = A(3,0)*tmp18 - A(3,2)*tmp27 + A(3,3)*tmp29 - A(3,5)*tmp31;
   tmp44 = A(3,0)*tmp19 - A(3,2)*tmp28 + A(3,3)*tmp30 - A(3,4)*tmp31;

   B(0,4) = - A(5,1)*tmp36 + A(5,2)*tmp37 - A(5,3)*tmp38 + A(5,4)*tmp39 - A(5,5)*tmp40;
   B(0,5) =   A(4,1)*tmp36 - A(4,2)*tmp37 + A(4,3)*tmp38 - A(4,4)*tmp39 + A(4,5)*tmp40;
   B(1,4) =   A(5,0)*tmp36 - A(5,2)*tmp41 + A(5,3)*tmp42 - A(5,4)*tmp43 + A(5,5)*tmp44;
   B(1,5) = - A(4,0)*tmp36 + A(4,2)*tmp41 - A(4,3)*tmp42 + A(4,4)*tmp43 - A(4,5)*tmp44;

   tmp45 = A(3,0)*tmp20 - A(3,1)*tmp26 + A(3,4)*tmp32 - A(3,5)*tmp33;
   tmp46 = A(3,0)*tmp21 - A(3,1)*tmp27 + A(3,3)*tmp32 - A(3,5)*tmp34;
   tmp47 = A(3,0)*tmp22 - A(3,1)*tmp28 + A(3,3)*tmp33 - A(3,4)*tmp34;
   tmp48 = A(3,0)*tmp23 - A(3,1)*tmp29 + A(3,2)*tmp32 - A(3,5)*tmp35;
   tmp49 = A(3,0)*tmp24 - A(3,1)*tmp30 + A(3,2)*tmp33 - A(3,4)*tmp35;

   B(2,4) = - A(5,0)*tmp37 + A(5,1)*tmp41 - A(5,3)*tmp45 + A(5,4)*tmp46 - A(5,5)*tmp47;
   B(2,5) =   A(4,0)*tmp37 - A(4,1)*tmp41 + A(4,3)*tmp45 - A(4,4)*tmp46 + A(4,5)*tmp47;
   B(3,4) =   A(5,0)*tmp38 - A(5,1)*tmp42 + A(5,2)*tmp45 - A(5,4)*tmp48 + A(5,5)*tmp49;
   B(3,5) = - A(4,0)*tmp38 + A(4,1)*tmp42 - A(4,2)*tmp45 + A(4,4)*tmp48 - A(4,5)*tmp49;

   tmp50 = A(3,0)*tmp25 - A(3,1)*tmp31 + A(3,2)*tmp34 - A(3,3)*tmp35;

   B(4,4) = - A(5,0)*tmp39 + A(5,1)*tmp43 - A(5,2)*tmp46 + A(5,3)*tmp48 - A(5,5)*tmp50;
   B(4,5) =   A(4,0)*tmp39 - A(4,1)*tmp43 + A(4,2)*tmp46 - A(4,3)*tmp48 + A(4,5)*tmp50;
   B(5,4) =   A(5,0)*tmp40 - A(5,1)*tmp44 + A(5,2)*tmp47 - A(5,3)*tmp49 + A(5,4)*tmp50;
   B(5,5) = - A(4,0)*tmp40 + A(4,1)*tmp44 - A(4,2)*tmp47 + A(4,3)*tmp49 - A(4,4)*tmp50;

   tmp36 = A(4,2)*tmp16 - A(4,3)*tmp17 + A(4,4)*tmp18 - A(4,5)*tmp19;
   tmp37 = A(4,1)*tmp16 - A(4,3)*tmp20 + A(4,4)*tmp21 - A(4,5)*tmp22;
   tmp38 = A(4,1)*tmp17 - A(4,2)*tmp20 + A(4,4)*tmp23 - A(4,5)*tmp24;
   tmp39 = A(4,1)*tmp18 - A(4,2)*tmp21 + A(4,3)*tmp23 - A(4,5)*tmp25;
   tmp40 = A(4,1)*tmp19 - A(4,2)*tmp22 + A(4,3)*tmp24 - A(4,4)*tmp25;
   tmp41 = A(4,0)*tmp16 - A(4,3)*tmp26 + A(4,4)*tmp27 - A(4,5)*tmp28;
   tmp42 = A(4,0)*tmp17 - A(4,2)*tmp26 + A(4,4)*tmp29 - A(4,5)*tmp30;
   tmp43 = A(4,0)*tmp18 - A(4,2)*tmp27 + A(4,3)*tmp29 - A(4,5)*tmp31;
   tmp44 = A(4,0)*tmp19 - A(4,2)*tmp28 + A(4,3)*tmp30 - A(4,4)*tmp31;
   tmp45 = A(4,0)*tmp20 - A(4,1)*tmp26 + A(4,4)*tmp32 - A(4,5)*tmp33;
   tmp46 = A(4,0)*tmp21 - A(4,1)*tmp27 + A(4,3)*tmp32 - A(4,5)*tmp34;
   tmp47 = A(4,0)*tmp22 - A(4,1)*tmp28 + A(4,3)*tmp33 - A(4,4)*tmp34;
   tmp48 = A(4,0)*tmp23 - A(4,1)*tmp29 + A(4,2)*tmp32 - A(4,5)*tmp35;
   tmp49 = A(4,0)*tmp24 - A(4,1)*tmp30 + A(4,2)*tmp33 - A(4,4)*tmp35;
   tmp50 = A(4,0)*tmp25 - A(4,1)*tmp31 + A(4,2)*tmp34 - A(4,3)*tmp35;

   B(0,3) =   A(5,1)*tmp36 - A(5,2)*tmp37 + A(5,3)*tmp38 - A(5,4)*tmp39 + A(5,5)*tmp40;
   B(1,3) = - A(5,0)*tmp36 + A(5,2)*tmp41 - A(5,3)*tmp42 + A(5,4)*tmp43 - A(5,5)*tmp44;
   B(2,3) =   A(5,0)*tmp37 - A(5,1)*tmp41 + A(5,3)*tmp45 - A(5,4)*tmp46 + A(5,5)*tmp47;
   B(3,3) = - A(5,0)*tmp38 + A(5,1)*tmp42 - A(5,2)*tmp45 + A(5,4)*tmp48 - A(5,5)*tmp49;
   B(4,3) =   A(5,0)*tmp39 - A(5,1)*tmp43 + A(5,2)*tmp46 - A(5,3)*tmp48 + A(5,5)*tmp50;
   B(5,3) = - A(5,0)*tmp40 + A(5,1)*tmp44 - A(5,2)*tmp47 + A(5,3)*tmp49 - A(5,4)*tmp50;

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) +
                 A(0,3)*B(3,0) + A(0,4)*B(4,0) + A(0,5)*B(5,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given symmetric dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The symmetric dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense \f$ 6 \times 6 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertSymmetric6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(4,4)*A(5,5) - A(4,5)*A(5,4) );
   ET tmp2 ( A(4,3)*A(5,5) - A(4,5)*A(5,3) );
   ET tmp3 ( A(4,3)*A(5,4) - A(4,4)*A(5,3) );
   ET tmp4 ( A(4,2)*A(5,5) - A(4,5)*A(5,2) );
   ET tmp5 ( A(4,2)*A(5,4) - A(4,4)*A(5,2) );
   ET tmp6 ( A(4,2)*A(5,3) - A(4,3)*A(5,2) );
   ET tmp7 ( A(4,1)*A(5,5) - A(4,5)*A(5,1) );
   ET tmp8 ( A(4,1)*A(5,4) - A(4,4)*A(5,1) );
   ET tmp9 ( A(4,1)*A(5,3) - A(4,3)*A(5,1) );
   ET tmp10( A(4,1)*A(5,2) - A(4,2)*A(5,1) );
   ET tmp11( A(4,0)*A(5,5) - A(4,5)*A(5,0) );
   ET tmp12( A(4,0)*A(5,4) - A(4,4)*A(5,0) );
   ET tmp13( A(4,0)*A(5,3) - A(4,3)*A(5,0) );
   ET tmp14( A(4,0)*A(5,2) - A(4,2)*A(5,0) );
   ET tmp15( A(4,0)*A(5,1) - A(4,1)*A(5,0) );

   ET tmp16( A(3,3)*tmp1  - A(3,4)*tmp2  + A(3,5)*tmp3  );
   ET tmp17( A(3,2)*tmp1  - A(3,4)*tmp4  + A(3,5)*tmp5  );
   ET tmp18( A(3,2)*tmp2  - A(3,3)*tmp4  + A(3,5)*tmp6  );
   ET tmp19( A(3,2)*tmp3  - A(3,3)*tmp5  + A(3,4)*tmp6  );
   ET tmp20( A(3,1)*tmp1  - A(3,4)*tmp7  + A(3,5)*tmp8  );
   ET tmp21( A(3,1)*tmp2  - A(3,3)*tmp7  + A(3,5)*tmp9  );
   ET tmp22( A(3,1)*tmp3  - A(3,3)*tmp8  + A(3,4)*tmp9  );
   ET tmp23( A(3,1)*tmp4  - A(3,2)*tmp7  + A(3,5)*tmp10 );
   ET tmp24( A(3,1)*tmp5  - A(3,2)*tmp8  + A(3,4)*tmp10 );
   ET tmp25( A(3,1)*tmp6  - A(3,2)*tmp9  + A(3,3)*tmp10 );
   ET tmp26( A(3,0)*tmp1  - A(3,4)*tmp11 + A(3,5)*tmp12 );
   ET tmp27( A(3,0)*tmp2  - A(3,3)*tmp11 + A(3,5)*tmp13 );
   ET tmp28( A(3,0)*tmp3  - A(3,3)*tmp12 + A(3,4)*tmp13 );
   ET tmp29( A(3,0)*tmp4  - A(3,2)*tmp11 + A(3,5)*tmp14 );
   ET tmp30( A(3,0)*tmp5  - A(3,2)*tmp12 + A(3,4)*tmp14 );
   ET tmp31( A(3,0)*tmp6  - A(3,2)*tmp13 + A(3,3)*tmp14 );
   ET tmp32( A(3,0)*tmp7  - A(3,1)*tmp11 + A(3,5)*tmp15 );
   ET tmp33( A(3,0)*tmp8  - A(3,1)*tmp12 + A(3,4)*tmp15 );
   ET tmp34( A(3,0)*tmp9  - A(3,1)*tmp13 + A(3,3)*tmp15 );
   ET tmp35( A(3,0)*tmp10 - A(3,1)*tmp14 + A(3,2)*tmp15 );

   ET tmp36( A(2,2)*tmp16 - A(2,3)*tmp17 + A(2,4)*tmp18 - A(2,5)*tmp19 );
   ET tmp37( A(2,1)*tmp16 - A(2,3)*tmp20 + A(2,4)*tmp21 - A(2,5)*tmp22 );
   ET tmp38( A(2,1)*tmp17 - A(2,2)*tmp20 + A(2,4)*tmp23 - A(2,5)*tmp24 );
   ET tmp39( A(2,1)*tmp18 - A(2,2)*tmp21 + A(2,3)*tmp23 - A(2,5)*tmp25 );
   ET tmp40( A(2,1)*tmp19 - A(2,2)*tmp22 + A(2,3)*tmp24 - A(2,4)*tmp25 );
   ET tmp41( A(2,0)*tmp16 - A(2,3)*tmp26 + A(2,4)*tmp27 - A(2,5)*tmp28 );
   ET tmp42( A(2,0)*tmp17 - A(2,2)*tmp26 + A(2,4)*tmp29 - A(2,5)*tmp30 );
   ET tmp43( A(2,0)*tmp18 - A(2,2)*tmp27 + A(2,3)*tmp29 - A(2,5)*tmp31 );
   ET tmp44( A(2,0)*tmp19 - A(2,2)*tmp28 + A(2,3)*tmp30 - A(2,4)*tmp31 );

   B(0,0) =   A(1,1)*tmp36 - A(1,2)*tmp37 + A(1,3)*tmp38 - A(1,4)*tmp39 + A(1,5)*tmp40;
   B(0,1) = - A(0,1)*tmp36 + A(0,2)*tmp37 - A(0,3)*tmp38 + A(0,4)*tmp39 - A(0,5)*tmp40;
   B(1,1) =   A(0,0)*tmp36 - A(0,2)*tmp41 + A(0,3)*tmp42 - A(0,4)*tmp43 + A(0,5)*tmp44;

   ET tmp45( A(2,0)*tmp20 - A(2,1)*tmp26 + A(2,4)*tmp32 - A(2,5)*tmp33 );
   ET tmp46( A(2,0)*tmp21 - A(2,1)*tmp27 + A(2,3)*tmp32 - A(2,5)*tmp34 );
   ET tmp47( A(2,0)*tmp22 - A(2,1)*tmp28 + A(2,3)*tmp33 - A(2,4)*tmp34 );
   ET tmp48( A(2,0)*tmp23 - A(2,1)*tmp29 + A(2,2)*tmp32 - A(2,5)*tmp35 );
   ET tmp49( A(2,0)*tmp24 - A(2,1)*tmp30 + A(2,2)*tmp33 - A(2,4)*tmp35 );

   B(2,0) =   A(1,0)*tmp37 - A(1,1)*tmp41 + A(1,3)*tmp45 - A(1,4)*tmp46 + A(1,5)*tmp47;
   B(2,1) = - A(0,0)*tmp37 + A(0,1)*tmp41 - A(0,3)*tmp45 + A(0,4)*tmp46 - A(0,5)*tmp47;
   B(3,0) = - A(1,0)*tmp38 + A(1,1)*tmp42 - A(1,2)*tmp45 + A(1,4)*tmp48 - A(1,5)*tmp49;
   B(3,1) =   A(0,0)*tmp38 - A(0,1)*tmp42 + A(0,2)*tmp45 - A(0,4)*tmp48 + A(0,5)*tmp49;

   ET tmp50( A(2,0)*tmp25 - A(2,1)*tmp31 + A(2,2)*tmp34 - A(2,3)*tmp35 );

   B(4,0) =   A(1,0)*tmp39 - A(1,1)*tmp43 + A(1,2)*tmp46 - A(1,3)*tmp48 + A(1,5)*tmp50;
   B(4,1) = - A(0,0)*tmp39 + A(0,1)*tmp43 - A(0,2)*tmp46 + A(0,3)*tmp48 - A(0,5)*tmp50;
   B(5,0) = - A(1,0)*tmp40 + A(1,1)*tmp44 - A(1,2)*tmp47 + A(1,3)*tmp49 - A(1,4)*tmp50;
   B(5,1) =   A(0,0)*tmp40 - A(0,1)*tmp44 + A(0,2)*tmp47 - A(0,3)*tmp49 + A(0,4)*tmp50;

   tmp36 = A(1,1)*tmp16 - A(1,3)*tmp20 + A(1,4)*tmp21 - A(1,5)*tmp22;
   tmp37 = A(1,1)*tmp17 - A(1,2)*tmp20 + A(1,4)*tmp23 - A(1,5)*tmp24;
   tmp38 = A(1,0)*tmp16 - A(1,3)*tmp26 + A(1,4)*tmp27 - A(1,5)*tmp28;
   tmp39 = A(1,0)*tmp17 - A(1,2)*tmp26 + A(1,4)*tmp29 - A(1,5)*tmp30;
   tmp40 = A(1,0)*tmp20 - A(1,1)*tmp26 + A(1,4)*tmp32 - A(1,5)*tmp33;
   tmp41 = A(1,0)*tmp21 - A(1,1)*tmp27 + A(1,3)*tmp32 - A(1,5)*tmp34;
   tmp42 = A(1,0)*tmp22 - A(1,1)*tmp28 + A(1,3)*tmp33 - A(1,4)*tmp34;
   tmp43 = A(1,0)*tmp23 - A(1,1)*tmp29 + A(1,2)*tmp32 - A(1,5)*tmp35;
   tmp44 = A(1,0)*tmp24 - A(1,1)*tmp30 + A(1,2)*tmp33 - A(1,4)*tmp35;

   B(2,2) =   A(0,0)*tmp36 - A(0,1)*tmp38 + A(0,3)*tmp40 - A(0,4)*tmp41 + A(0,5)*tmp42;
   B(3,2) = - A(0,0)*tmp37 + A(0,1)*tmp39 - A(0,2)*tmp40 + A(0,4)*tmp43 - A(0,5)*tmp44;

   tmp1  = A(0,3)*A(1,4) - A(0,4)*A(1,3);
   tmp2  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp3  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp4  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp5  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp6  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp7  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp8  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp9  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp10 = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp11 = A(0,3)*A(1,5) - A(0,5)*A(1,3);
   tmp12 = A(0,2)*A(1,5) - A(0,5)*A(1,2);
   tmp13 = A(0,1)*A(1,5) - A(0,5)*A(1,1);
   tmp14 = A(0,0)*A(1,5) - A(0,5)*A(1,0);
   tmp15 = A(0,4)*A(1,5) - A(0,5)*A(1,4);

   tmp16 = A(2,3)*tmp15 - A(2,4)*tmp11 + A(2,5)*tmp1;
   tmp17 = A(2,2)*tmp15 - A(2,4)*tmp12 + A(2,5)*tmp2;
   tmp18 = A(2,2)*tmp11 - A(2,3)*tmp12 + A(2,5)*tmp3;
   tmp19 = A(2,2)*tmp1  - A(2,3)*tmp2  + A(2,4)*tmp3;
   tmp20 = A(2,1)*tmp15 - A(2,4)*tmp13 + A(2,5)*tmp4;
   tmp21 = A(2,1)*tmp11 - A(2,3)*tmp13 + A(2,5)*tmp5;
   tmp22 = A(2,1)*tmp1  - A(2,3)*tmp4  + A(2,4)*tmp5;
   tmp23 = A(2,1)*tmp12 - A(2,2)*tmp13 + A(2,5)*tmp6;
   tmp24 = A(2,1)*tmp2  - A(2,2)*tmp4  + A(2,4)*tmp6;
   tmp25 = A(2,1)*tmp3  - A(2,2)*tmp5  + A(2,3)*tmp6;
   tmp26 = A(2,0)*tmp15 - A(2,4)*tmp14 + A(2,5)*tmp7;
   tmp27 = A(2,0)*tmp11 - A(2,3)*tmp14 + A(2,5)*tmp8;
   tmp28 = A(2,0)*tmp1  - A(2,3)*tmp7  + A(2,4)*tmp8;
   tmp29 = A(2,0)*tmp12 - A(2,2)*tmp14 + A(2,5)*tmp9;
   tmp30 = A(2,0)*tmp2  - A(2,2)*tmp7  + A(2,4)*tmp9;
   tmp31 = A(2,0)*tmp3  - A(2,2)*tmp8  + A(2,3)*tmp9;
   tmp32 = A(2,0)*tmp13 - A(2,1)*tmp14 + A(2,5)*tmp10;
   tmp33 = A(2,0)*tmp4  - A(2,1)*tmp7  + A(2,4)*tmp10;
   tmp34 = A(2,0)*tmp5  - A(2,1)*tmp8  + A(2,3)*tmp10;
   tmp35 = A(2,0)*tmp6  - A(2,1)*tmp9  + A(2,2)*tmp10;

   tmp36 = A(3,1)*tmp16 - A(3,3)*tmp20 + A(3,4)*tmp21 - A(3,5)*tmp22;
   tmp37 = A(3,1)*tmp17 - A(3,2)*tmp20 + A(3,4)*tmp23 - A(3,5)*tmp24;
   tmp38 = A(3,0)*tmp16 - A(3,3)*tmp26 + A(3,4)*tmp27 - A(3,5)*tmp28;
   tmp39 = A(3,0)*tmp17 - A(3,2)*tmp26 + A(3,4)*tmp29 - A(3,5)*tmp30;
   tmp40 = A(3,0)*tmp20 - A(3,1)*tmp26 + A(3,4)*tmp32 - A(3,5)*tmp33;
   tmp41 = A(3,0)*tmp21 - A(3,1)*tmp27 + A(3,3)*tmp32 - A(3,5)*tmp34;
   tmp42 = A(3,0)*tmp22 - A(3,1)*tmp28 + A(3,3)*tmp33 - A(3,4)*tmp34;
   tmp43 = A(3,0)*tmp23 - A(3,1)*tmp29 + A(3,2)*tmp32 - A(3,5)*tmp35;
   tmp44 = A(3,0)*tmp24 - A(3,1)*tmp30 + A(3,2)*tmp33 - A(3,4)*tmp35;

   B(2,4) = - A(5,0)*tmp36 + A(5,1)*tmp38 - A(5,3)*tmp40 + A(5,4)*tmp41 - A(5,5)*tmp42;
   B(2,5) =   A(4,0)*tmp36 - A(4,1)*tmp38 + A(4,3)*tmp40 - A(4,4)*tmp41 + A(4,5)*tmp42;
   B(3,4) =   A(5,0)*tmp37 - A(5,1)*tmp39 + A(5,2)*tmp40 - A(5,4)*tmp43 + A(5,5)*tmp44;
   B(3,5) = - A(4,0)*tmp37 + A(4,1)*tmp39 - A(4,2)*tmp40 + A(4,4)*tmp43 - A(4,5)*tmp44;

   tmp36 = A(3,1)*tmp18 - A(3,2)*tmp21 + A(3,3)*tmp23 - A(3,5)*tmp25;
   tmp37 = A(3,1)*tmp19 - A(3,2)*tmp22 + A(3,3)*tmp24 - A(3,4)*tmp25;
   tmp38 = A(3,0)*tmp18 - A(3,2)*tmp27 + A(3,3)*tmp29 - A(3,5)*tmp31;
   tmp39 = A(3,0)*tmp19 - A(3,2)*tmp28 + A(3,3)*tmp30 - A(3,4)*tmp31;
   tmp40 = A(3,0)*tmp25 - A(3,1)*tmp31 + A(3,2)*tmp34 - A(3,3)*tmp35;

   B(4,4) = - A(5,0)*tmp36 + A(5,1)*tmp38 - A(5,2)*tmp41 + A(5,3)*tmp43 - A(5,5)*tmp40;
   B(4,5) =   A(4,0)*tmp36 - A(4,1)*tmp38 + A(4,2)*tmp41 - A(4,3)*tmp43 + A(4,5)*tmp40;
   B(5,5) = - A(4,0)*tmp37 + A(4,1)*tmp39 - A(4,2)*tmp42 + A(4,3)*tmp44 - A(4,4)*tmp40;

   tmp36 = A(4,1)*tmp17 - A(4,2)*tmp20 + A(4,4)*tmp23 - A(4,5)*tmp24;
   tmp37 = A(4,0)*tmp17 - A(4,2)*tmp26 + A(4,4)*tmp29 - A(4,5)*tmp30;
   tmp38 = A(4,0)*tmp20 - A(4,1)*tmp26 + A(4,4)*tmp32 - A(4,5)*tmp33;
   tmp39 = A(4,0)*tmp23 - A(4,1)*tmp29 + A(4,2)*tmp32 - A(4,5)*tmp35;
   tmp40 = A(4,0)*tmp24 - A(4,1)*tmp30 + A(4,2)*tmp33 - A(4,4)*tmp35;

   B(3,3) = - A(5,0)*tmp36 + A(5,1)*tmp37 - A(5,2)*tmp38 + A(5,4)*tmp39 - A(5,5)*tmp40;

   B(0,2) = B(2,0);
   B(0,3) = B(3,0);
   B(0,4) = B(4,0);
   B(0,5) = B(5,0);
   B(1,0) = B(0,1);
   B(1,2) = B(2,1);
   B(1,3) = B(3,1);
   B(1,4) = B(4,1);
   B(1,5) = B(5,1);
   B(2,3) = B(3,2);
   B(4,2) = B(2,4);
   B(4,3) = B(3,4);
   B(5,2) = B(2,5);
   B(5,3) = B(3,5);
   B(5,4) = B(4,5);

   const ET det( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) +
                 A(0,3)*B(3,0) + A(0,4)*B(4,0) + A(0,5)*B(5,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given Hermitian dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The Hermitian dense matrix to be inverted.
// \return void
//
// This function inverts the given Hermitian dense \f$ 6 \times 6 \f$ matrix via the rule of
// Sarrus. The matrix inversion fails if the given matrix is singular and not invertible. In
// this case a \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertHermitian6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   ET tmp1 ( A(4,4)*A(5,5) - A(4,5)*A(5,4) );
   ET tmp2 ( A(4,3)*A(5,5) - A(4,5)*A(5,3) );
   ET tmp3 ( A(4,3)*A(5,4) - A(4,4)*A(5,3) );
   ET tmp4 ( A(4,2)*A(5,5) - A(4,5)*A(5,2) );
   ET tmp5 ( A(4,2)*A(5,4) - A(4,4)*A(5,2) );
   ET tmp6 ( A(4,2)*A(5,3) - A(4,3)*A(5,2) );
   ET tmp7 ( A(4,1)*A(5,5) - A(4,5)*A(5,1) );
   ET tmp8 ( A(4,1)*A(5,4) - A(4,4)*A(5,1) );
   ET tmp9 ( A(4,1)*A(5,3) - A(4,3)*A(5,1) );
   ET tmp10( A(4,1)*A(5,2) - A(4,2)*A(5,1) );
   ET tmp11( A(4,0)*A(5,5) - A(4,5)*A(5,0) );
   ET tmp12( A(4,0)*A(5,4) - A(4,4)*A(5,0) );
   ET tmp13( A(4,0)*A(5,3) - A(4,3)*A(5,0) );
   ET tmp14( A(4,0)*A(5,2) - A(4,2)*A(5,0) );
   ET tmp15( A(4,0)*A(5,1) - A(4,1)*A(5,0) );

   ET tmp16( A(3,3)*tmp1  - A(3,4)*tmp2  + A(3,5)*tmp3  );
   ET tmp17( A(3,2)*tmp1  - A(3,4)*tmp4  + A(3,5)*tmp5  );
   ET tmp18( A(3,2)*tmp2  - A(3,3)*tmp4  + A(3,5)*tmp6  );
   ET tmp19( A(3,2)*tmp3  - A(3,3)*tmp5  + A(3,4)*tmp6  );
   ET tmp20( A(3,1)*tmp1  - A(3,4)*tmp7  + A(3,5)*tmp8  );
   ET tmp21( A(3,1)*tmp2  - A(3,3)*tmp7  + A(3,5)*tmp9  );
   ET tmp22( A(3,1)*tmp3  - A(3,3)*tmp8  + A(3,4)*tmp9  );
   ET tmp23( A(3,1)*tmp4  - A(3,2)*tmp7  + A(3,5)*tmp10 );
   ET tmp24( A(3,1)*tmp5  - A(3,2)*tmp8  + A(3,4)*tmp10 );
   ET tmp25( A(3,1)*tmp6  - A(3,2)*tmp9  + A(3,3)*tmp10 );
   ET tmp26( A(3,0)*tmp1  - A(3,4)*tmp11 + A(3,5)*tmp12 );
   ET tmp27( A(3,0)*tmp2  - A(3,3)*tmp11 + A(3,5)*tmp13 );
   ET tmp28( A(3,0)*tmp3  - A(3,3)*tmp12 + A(3,4)*tmp13 );
   ET tmp29( A(3,0)*tmp4  - A(3,2)*tmp11 + A(3,5)*tmp14 );
   ET tmp30( A(3,0)*tmp5  - A(3,2)*tmp12 + A(3,4)*tmp14 );
   ET tmp31( A(3,0)*tmp6  - A(3,2)*tmp13 + A(3,3)*tmp14 );
   ET tmp32( A(3,0)*tmp7  - A(3,1)*tmp11 + A(3,5)*tmp15 );
   ET tmp33( A(3,0)*tmp8  - A(3,1)*tmp12 + A(3,4)*tmp15 );
   ET tmp34( A(3,0)*tmp9  - A(3,1)*tmp13 + A(3,3)*tmp15 );
   ET tmp35( A(3,0)*tmp10 - A(3,1)*tmp14 + A(3,2)*tmp15 );

   ET tmp36( A(2,2)*tmp16 - A(2,3)*tmp17 + A(2,4)*tmp18 - A(2,5)*tmp19 );
   ET tmp37( A(2,1)*tmp16 - A(2,3)*tmp20 + A(2,4)*tmp21 - A(2,5)*tmp22 );
   ET tmp38( A(2,1)*tmp17 - A(2,2)*tmp20 + A(2,4)*tmp23 - A(2,5)*tmp24 );
   ET tmp39( A(2,1)*tmp18 - A(2,2)*tmp21 + A(2,3)*tmp23 - A(2,5)*tmp25 );
   ET tmp40( A(2,1)*tmp19 - A(2,2)*tmp22 + A(2,3)*tmp24 - A(2,4)*tmp25 );
   ET tmp41( A(2,0)*tmp16 - A(2,3)*tmp26 + A(2,4)*tmp27 - A(2,5)*tmp28 );
   ET tmp42( A(2,0)*tmp17 - A(2,2)*tmp26 + A(2,4)*tmp29 - A(2,5)*tmp30 );
   ET tmp43( A(2,0)*tmp18 - A(2,2)*tmp27 + A(2,3)*tmp29 - A(2,5)*tmp31 );
   ET tmp44( A(2,0)*tmp19 - A(2,2)*tmp28 + A(2,3)*tmp30 - A(2,4)*tmp31 );

   B(0,0) = ET( real( A(1,1)*tmp36 - A(1,2)*tmp37 + A(1,3)*tmp38 - A(1,4)*tmp39 + A(1,5)*tmp40 ) );
   B(0,1) = - A(0,1)*tmp36 + A(0,2)*tmp37 - A(0,3)*tmp38 + A(0,4)*tmp39 - A(0,5)*tmp40;
   B(1,1) = ET( real( A(0,0)*tmp36 - A(0,2)*tmp41 + A(0,3)*tmp42 - A(0,4)*tmp43 + A(0,5)*tmp44 ) );

   ET tmp45( A(2,0)*tmp20 - A(2,1)*tmp26 + A(2,4)*tmp32 - A(2,5)*tmp33 );
   ET tmp46( A(2,0)*tmp21 - A(2,1)*tmp27 + A(2,3)*tmp32 - A(2,5)*tmp34 );
   ET tmp47( A(2,0)*tmp22 - A(2,1)*tmp28 + A(2,3)*tmp33 - A(2,4)*tmp34 );
   ET tmp48( A(2,0)*tmp23 - A(2,1)*tmp29 + A(2,2)*tmp32 - A(2,5)*tmp35 );
   ET tmp49( A(2,0)*tmp24 - A(2,1)*tmp30 + A(2,2)*tmp33 - A(2,4)*tmp35 );

   B(2,0) =   A(1,0)*tmp37 - A(1,1)*tmp41 + A(1,3)*tmp45 - A(1,4)*tmp46 + A(1,5)*tmp47;
   B(2,1) = - A(0,0)*tmp37 + A(0,1)*tmp41 - A(0,3)*tmp45 + A(0,4)*tmp46 - A(0,5)*tmp47;
   B(3,0) = - A(1,0)*tmp38 + A(1,1)*tmp42 - A(1,2)*tmp45 + A(1,4)*tmp48 - A(1,5)*tmp49;
   B(3,1) =   A(0,0)*tmp38 - A(0,1)*tmp42 + A(0,2)*tmp45 - A(0,4)*tmp48 + A(0,5)*tmp49;

   ET tmp50( A(2,0)*tmp25 - A(2,1)*tmp31 + A(2,2)*tmp34 - A(2,3)*tmp35 );

   B(4,0) =   A(1,0)*tmp39 - A(1,1)*tmp43 + A(1,2)*tmp46 - A(1,3)*tmp48 + A(1,5)*tmp50;
   B(4,1) = - A(0,0)*tmp39 + A(0,1)*tmp43 - A(0,2)*tmp46 + A(0,3)*tmp48 - A(0,5)*tmp50;
   B(5,0) = - A(1,0)*tmp40 + A(1,1)*tmp44 - A(1,2)*tmp47 + A(1,3)*tmp49 - A(1,4)*tmp50;
   B(5,1) =   A(0,0)*tmp40 - A(0,1)*tmp44 + A(0,2)*tmp47 - A(0,3)*tmp49 + A(0,4)*tmp50;

   tmp36 = A(1,1)*tmp16 - A(1,3)*tmp20 + A(1,4)*tmp21 - A(1,5)*tmp22;
   tmp37 = A(1,1)*tmp17 - A(1,2)*tmp20 + A(1,4)*tmp23 - A(1,5)*tmp24;
   tmp38 = A(1,0)*tmp16 - A(1,3)*tmp26 + A(1,4)*tmp27 - A(1,5)*tmp28;
   tmp39 = A(1,0)*tmp17 - A(1,2)*tmp26 + A(1,4)*tmp29 - A(1,5)*tmp30;
   tmp40 = A(1,0)*tmp20 - A(1,1)*tmp26 + A(1,4)*tmp32 - A(1,5)*tmp33;
   tmp41 = A(1,0)*tmp21 - A(1,1)*tmp27 + A(1,3)*tmp32 - A(1,5)*tmp34;
   tmp42 = A(1,0)*tmp22 - A(1,1)*tmp28 + A(1,3)*tmp33 - A(1,4)*tmp34;
   tmp43 = A(1,0)*tmp23 - A(1,1)*tmp29 + A(1,2)*tmp32 - A(1,5)*tmp35;
   tmp44 = A(1,0)*tmp24 - A(1,1)*tmp30 + A(1,2)*tmp33 - A(1,4)*tmp35;

   B(2,2) = ET( real( A(0,0)*tmp36 - A(0,1)*tmp38 + A(0,3)*tmp40 - A(0,4)*tmp41 + A(0,5)*tmp42 ) );
   B(3,2) = - A(0,0)*tmp37 + A(0,1)*tmp39 - A(0,2)*tmp40 + A(0,4)*tmp43 - A(0,5)*tmp44;

   tmp1  = A(0,3)*A(1,4) - A(0,4)*A(1,3);
   tmp2  = A(0,2)*A(1,4) - A(0,4)*A(1,2);
   tmp3  = A(0,2)*A(1,3) - A(0,3)*A(1,2);
   tmp4  = A(0,1)*A(1,4) - A(0,4)*A(1,1);
   tmp5  = A(0,1)*A(1,3) - A(0,3)*A(1,1);
   tmp6  = A(0,1)*A(1,2) - A(0,2)*A(1,1);
   tmp7  = A(0,0)*A(1,4) - A(0,4)*A(1,0);
   tmp8  = A(0,0)*A(1,3) - A(0,3)*A(1,0);
   tmp9  = A(0,0)*A(1,2) - A(0,2)*A(1,0);
   tmp10 = A(0,0)*A(1,1) - A(0,1)*A(1,0);
   tmp11 = A(0,3)*A(1,5) - A(0,5)*A(1,3);
   tmp12 = A(0,2)*A(1,5) - A(0,5)*A(1,2);
   tmp13 = A(0,1)*A(1,5) - A(0,5)*A(1,1);
   tmp14 = A(0,0)*A(1,5) - A(0,5)*A(1,0);
   tmp15 = A(0,4)*A(1,5) - A(0,5)*A(1,4);

   tmp16 = A(2,3)*tmp15 - A(2,4)*tmp11 + A(2,5)*tmp1;
   tmp17 = A(2,2)*tmp15 - A(2,4)*tmp12 + A(2,5)*tmp2;
   tmp18 = A(2,2)*tmp11 - A(2,3)*tmp12 + A(2,5)*tmp3;
   tmp19 = A(2,2)*tmp1  - A(2,3)*tmp2  + A(2,4)*tmp3;
   tmp20 = A(2,1)*tmp15 - A(2,4)*tmp13 + A(2,5)*tmp4;
   tmp21 = A(2,1)*tmp11 - A(2,3)*tmp13 + A(2,5)*tmp5;
   tmp22 = A(2,1)*tmp1  - A(2,3)*tmp4  + A(2,4)*tmp5;
   tmp23 = A(2,1)*tmp12 - A(2,2)*tmp13 + A(2,5)*tmp6;
   tmp24 = A(2,1)*tmp2  - A(2,2)*tmp4  + A(2,4)*tmp6;
   tmp25 = A(2,1)*tmp3  - A(2,2)*tmp5  + A(2,3)*tmp6;
   tmp26 = A(2,0)*tmp15 - A(2,4)*tmp14 + A(2,5)*tmp7;
   tmp27 = A(2,0)*tmp11 - A(2,3)*tmp14 + A(2,5)*tmp8;
   tmp28 = A(2,0)*tmp1  - A(2,3)*tmp7  + A(2,4)*tmp8;
   tmp29 = A(2,0)*tmp12 - A(2,2)*tmp14 + A(2,5)*tmp9;
   tmp30 = A(2,0)*tmp2  - A(2,2)*tmp7  + A(2,4)*tmp9;
   tmp31 = A(2,0)*tmp3  - A(2,2)*tmp8  + A(2,3)*tmp9;
   tmp32 = A(2,0)*tmp13 - A(2,1)*tmp14 + A(2,5)*tmp10;
   tmp33 = A(2,0)*tmp4  - A(2,1)*tmp7  + A(2,4)*tmp10;
   tmp34 = A(2,0)*tmp5  - A(2,1)*tmp8  + A(2,3)*tmp10;
   tmp35 = A(2,0)*tmp6  - A(2,1)*tmp9  + A(2,2)*tmp10;

   tmp36 = A(3,1)*tmp16 - A(3,3)*tmp20 + A(3,4)*tmp21 - A(3,5)*tmp22;
   tmp37 = A(3,1)*tmp17 - A(3,2)*tmp20 + A(3,4)*tmp23 - A(3,5)*tmp24;
   tmp38 = A(3,0)*tmp16 - A(3,3)*tmp26 + A(3,4)*tmp27 - A(3,5)*tmp28;
   tmp39 = A(3,0)*tmp17 - A(3,2)*tmp26 + A(3,4)*tmp29 - A(3,5)*tmp30;
   tmp40 = A(3,0)*tmp20 - A(3,1)*tmp26 + A(3,4)*tmp32 - A(3,5)*tmp33;
   tmp41 = A(3,0)*tmp21 - A(3,1)*tmp27 + A(3,3)*tmp32 - A(3,5)*tmp34;
   tmp42 = A(3,0)*tmp22 - A(3,1)*tmp28 + A(3,3)*tmp33 - A(3,4)*tmp34;
   tmp43 = A(3,0)*tmp23 - A(3,1)*tmp29 + A(3,2)*tmp32 - A(3,5)*tmp35;
   tmp44 = A(3,0)*tmp24 - A(3,1)*tmp30 + A(3,2)*tmp33 - A(3,4)*tmp35;

   B(2,4) = - A(5,0)*tmp36 + A(5,1)*tmp38 - A(5,3)*tmp40 + A(5,4)*tmp41 - A(5,5)*tmp42;
   B(2,5) =   A(4,0)*tmp36 - A(4,1)*tmp38 + A(4,3)*tmp40 - A(4,4)*tmp41 + A(4,5)*tmp42;
   B(3,4) =   A(5,0)*tmp37 - A(5,1)*tmp39 + A(5,2)*tmp40 - A(5,4)*tmp43 + A(5,5)*tmp44;
   B(3,5) = - A(4,0)*tmp37 + A(4,1)*tmp39 - A(4,2)*tmp40 + A(4,4)*tmp43 - A(4,5)*tmp44;

   tmp36 = A(3,1)*tmp18 - A(3,2)*tmp21 + A(3,3)*tmp23 - A(3,5)*tmp25;
   tmp37 = A(3,1)*tmp19 - A(3,2)*tmp22 + A(3,3)*tmp24 - A(3,4)*tmp25;
   tmp38 = A(3,0)*tmp18 - A(3,2)*tmp27 + A(3,3)*tmp29 - A(3,5)*tmp31;
   tmp39 = A(3,0)*tmp19 - A(3,2)*tmp28 + A(3,3)*tmp30 - A(3,4)*tmp31;
   tmp40 = A(3,0)*tmp25 - A(3,1)*tmp31 + A(3,2)*tmp34 - A(3,3)*tmp35;

   B(4,4) = - ET( real( A(5,0)*tmp36 - A(5,1)*tmp38 + A(5,2)*tmp41 - A(5,3)*tmp43 + A(5,5)*tmp40 ) );
   B(4,5) =   A(4,0)*tmp36 - A(4,1)*tmp38 + A(4,2)*tmp41 - A(4,3)*tmp43 + A(4,5)*tmp40;
   B(5,5) = - ET( real( A(4,0)*tmp37 - A(4,1)*tmp39 + A(4,2)*tmp42 - A(4,3)*tmp44 + A(4,4)*tmp40 ) );

   tmp36 = A(4,1)*tmp17 - A(4,2)*tmp20 + A(4,4)*tmp23 - A(4,5)*tmp24;
   tmp37 = A(4,0)*tmp17 - A(4,2)*tmp26 + A(4,4)*tmp29 - A(4,5)*tmp30;
   tmp38 = A(4,0)*tmp20 - A(4,1)*tmp26 + A(4,4)*tmp32 - A(4,5)*tmp33;
   tmp39 = A(4,0)*tmp23 - A(4,1)*tmp29 + A(4,2)*tmp32 - A(4,5)*tmp35;
   tmp40 = A(4,0)*tmp24 - A(4,1)*tmp30 + A(4,2)*tmp33 - A(4,4)*tmp35;

   B(3,3) = - ET( real( A(5,0)*tmp36 - A(5,1)*tmp37 + A(5,2)*tmp38 - A(5,4)*tmp39 + A(5,5)*tmp40 ) );

   B(0,2) = conj( B(2,0) );
   B(0,3) = conj( B(3,0) );
   B(0,4) = conj( B(4,0) );
   B(0,5) = conj( B(5,0) );
   B(1,0) = conj( B(0,1) );
   B(1,2) = conj( B(2,1) );
   B(1,3) = conj( B(3,1) );
   B(1,4) = conj( B(4,1) );
   B(1,5) = conj( B(5,1) );
   B(2,3) = conj( B(3,2) );
   B(4,2) = conj( B(2,4) );
   B(4,3) = conj( B(3,4) );
   B(5,2) = conj( B(2,5) );
   B(5,3) = conj( B(3,5) );
   B(5,4) = conj( B(4,5) );

   const ET det( real( A(0,0)*B(0,0) + A(0,1)*B(1,0) + A(0,2)*B(2,0) +
                       A(0,3)*B(3,0) + A(0,4)*B(4,0) + A(0,5)*B(5,0) ) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The lower dense matrix to be inverted.
// \return void
//
// This function inverts the given lower dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLower6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(4,4)*A(5,5) );
   const ET tmp2( A(4,3)*A(5,5) );
   const ET tmp3( A(4,3)*A(5,4) - A(4,4)*A(5,3) );

   const ET tmp4( A(3,3)*tmp1 );
   const ET tmp5( A(3,2)*tmp1 );
   const ET tmp6( A(3,2)*tmp2 - A(3,3)*A(4,2)*A(5,5) );
   const ET tmp7( A(3,2)*tmp3 - A(3,3)*( A(4,2)*A(5,4) - A(4,4)*A(5,2) ) );
   const ET tmp8( A(0,0)*A(1,1)*A(2,2) );

   const ET tmp9 ( A(2,2)*tmp4 );
   const ET tmp10( A(2,1)*tmp4 );
   const ET tmp11( A(2,1)*tmp5 - A(2,2)*A(3,1)*tmp1 );
   const ET tmp12( A(2,1)*tmp6 - A(2,2)*( A(3,1)*tmp2 - A(3,3)*A(4,1)*A(5,5) ) );
   const ET tmp13( A(2,1)*tmp7 - A(2,2)*( A(3,1)*tmp3 - A(3,3)*( A(4,1)*A(5,4) - A(4,4)*A(5,1) ) ) );
   const ET tmp14( A(4,4)*tmp8 );
   const ET tmp15( A(4,3)*tmp8 );
   const ET tmp16( A(3,3)*tmp8 );

   B(0,0) =   A(1,1)*tmp9;
   B(1,0) = - A(1,0)*tmp9;
   B(2,0) =   A(1,0)*tmp10 - A(1,1)*A(2,0)*tmp4;
   B(3,0) = - A(1,0)*tmp11 + A(1,1)*( A(2,0)*tmp5 - A(2,2)*A(3,0)*tmp1 );
   B(4,0) =   A(1,0)*tmp12 - A(1,1)*( A(2,0)*tmp6 - A(2,2)*( A(3,0)*tmp2 - A(3,3)*A(4,0)*A(5,5) ) );
   B(5,0) = - A(1,0)*tmp13 + A(1,1)*( A(2,0)*tmp7 - A(2,2)*( A(3,0)*tmp3 - A(3,3)*( A(4,0)*A(5,4) - A(4,4)*A(5,0) ) ) );
   B(1,1) =   A(0,0)*tmp9;
   B(2,1) = - A(0,0)*tmp10;
   B(3,1) =   A(0,0)*tmp11;
   B(4,1) = - A(0,0)*tmp12;
   B(5,1) =   A(0,0)*tmp13;
   B(2,2) =   A(0,0)*A(1,1)*tmp4;
   B(3,2) = - A(0,0)*A(1,1)*tmp5;
   B(4,2) =   A(0,0)*A(1,1)*tmp6;
   B(5,2) = - A(0,0)*A(1,1)*tmp7;
   B(3,3) =   A(5,5)*tmp14;
   B(4,3) = - A(5,5)*tmp15;
   B(5,3) =   A(5,4)*tmp15 - A(5,3)*tmp14;
   B(4,4) =   A(5,5)*tmp16 - A(5,3)*A(3,5)*tmp8;
   B(5,4) = - A(5,4)*tmp16;
   B(5,5) =   A(4,4)*tmp16;

   const ET det( B(5,5)*A(5,5) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given unilower dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The unilower dense matrix to be inverted.
// \return void
//
// This function inverts the given unilower dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLower6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(4,3)*A(5,4) - A(5,3) );
   const ET tmp2( A(3,2)*A(4,3) - A(4,2) );
   const ET tmp3( A(3,2)*tmp1 - A(4,2)*A(5,4) + A(5,2) );
   const ET tmp4( A(2,1)*A(3,2) - A(3,1) );
   const ET tmp5( A(2,1)*tmp2 - A(3,1)*A(4,3) + A(4,1) );
   const ET tmp6( A(2,1)*tmp3 - A(3,1)*tmp1 + A(4,1)*A(5,4) - A(5,1) );

   B(1,0) = - A(1,0);
   B(2,0) =   A(1,0)*A(2,1) - A(2,0);
   B(3,0) = - A(1,0)*tmp4 + A(2,0)*A(3,2) - A(3,0);
   B(4,0) =   A(1,0)*tmp5 - A(2,0)*tmp2 + A(3,0)*A(4,3) - A(4,0);
   B(5,0) = - A(1,0)*tmp6 + A(2,0)*tmp3 - A(3,0)*tmp1 + A(4,0)*A(5,4) - A(5,0);
   B(2,1) = - A(2,1);
   B(3,1) =   tmp4;
   B(4,1) = - tmp5;
   B(5,1) =   tmp6;
   B(3,2) = - A(3,2);
   B(4,2) =   tmp2;
   B(5,2) = - tmp3;
   B(4,3) =   A(5,3)*A(4,5) - A(4,3);
   B(5,3) =   A(5,4)*A(4,3) - A(5,3);
   B(5,4) =   A(5,3)*A(3,4) - A(5,4);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The upper dense matrix to be inverted.
// \return void
//
// This function inverts the given upper dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpper6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(0,1)*A(1,2) - A(0,2)*A(1,1) );
   const ET tmp2( A(0,0)*A(1,2) );
   const ET tmp3( A(0,0)*A(1,1) );

   const ET tmp4( A(3,3)*A(4,4)*A(5,5)  );
   const ET tmp5( A(2,3)*tmp1 - A(2,2)*( A(0,1)*A(1,3) - A(0,3)*A(1,1) ) );
   const ET tmp6( A(2,3)*tmp2 - A(0,0)*A(1,3)*A(2,2) );
   const ET tmp7( A(2,3)*tmp3 );
   const ET tmp8( A(2,2)*tmp3 );

   const ET tmp9 ( A(2,2)*tmp4 );
   const ET tmp10( A(1,2)*tmp4 );
   const ET tmp11( A(1,1)*tmp4 );
   const ET tmp12( A(3,3)*( A(2,4)*tmp1 - A(2,2)*( A(0,1)*A(1,4) - A(0,4)*A(1,1) ) ) - A(3,4)*tmp5 );
   const ET tmp13( A(3,3)*( A(2,4)*tmp2 - A(0,0)*A(1,4)*A(2,2) ) - A(3,4)*tmp6 );
   const ET tmp14( A(3,3)*A(2,4)*tmp3 - A(3,4)*tmp7 );
   const ET tmp15( - A(3,4)*tmp8 );
   const ET tmp16( - A(3,3)*tmp8 );

   B(0,0) =   A(1,1)*tmp9;
   B(0,1) = - A(0,1)*tmp9;
   B(1,1) =   A(0,0)*tmp9;
   B(0,2) =   A(0,1)*tmp10 - A(0,2)*tmp11;
   B(1,2) = - A(0,0)*tmp10;
   B(2,2) =   A(0,0)*tmp11;
   B(0,3) = - A(5,5)*A(4,4)*tmp5;
   B(1,3) =   A(5,5)*A(4,4)*tmp6;
   B(2,3) = - A(5,5)*A(4,4)*tmp7;
   B(3,3) =   A(5,5)*A(4,4)*tmp8;
   B(0,4) = - A(5,5)*tmp12;
   B(1,4) =   A(5,5)*tmp13;
   B(2,4) = - A(5,5)*tmp14;
   B(3,4) =   A(5,5)*tmp15;
   B(4,4) = - A(5,5)*tmp16;
   B(0,5) = - A(4,4)*( A(3,3)*( A(2,5)*tmp1 - A(2,2)*( A(0,1)*A(1,5) - A(0,5)*A(1,1) ) ) - A(3,5)*tmp5 ) + A(4,5)*tmp12;
   B(1,5) =   A(4,4)*( A(3,3)*( A(2,5)*tmp2 - A(0,0)*A(1,5)*A(2,2) ) - A(3,5)*tmp6 ) - A(4,5)*tmp13;
   B(2,5) = - A(4,4)*( A(3,3)*A(2,5)*tmp3 - A(3,5)*tmp7 ) + A(4,5)*tmp14;
   B(3,5) = - A(4,4)*A(3,5)*tmp8 - A(4,5)*tmp15;
   B(4,5) =   A(4,5)*tmp16;
   B(5,5) = - A(4,4)*tmp16;

   const ET det( A(0,0)*B(0,0) );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   B /= det;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given uniupper dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The uniupper dense matrix to be inverted.
// \return void
//
// This function inverts the given uniupper dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpper6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   const StaticMatrix<ET,6UL,6UL,SO> A( *dm );
   MT& B( *dm );

   const ET tmp1( A(0,1)*A(1,2) - A(0,2) );
   const ET tmp2( A(2,3)*tmp1 - A(0,1)*A(1,3) + A(0,3) );
   const ET tmp3( A(2,3)*A(1,2) - A(1,3) );
   const ET tmp4( A(2,4)*tmp1 - A(0,1)*A(1,4) + A(0,4) - A(3,4)*tmp2 );
   const ET tmp5( A(2,4)*A(1,2) - A(1,4) - A(3,4)*tmp3 );
   const ET tmp6( A(2,4) - A(3,4)*A(2,3) );

   B(0,1) = - A(0,1);
   B(0,2) =   A(0,1)*A(1,2) - A(0,2);
   B(1,2) = - A(1,2);
   B(0,3) = - tmp2;
   B(1,3) =   tmp3;
   B(2,3) = - A(2,3);
   B(0,4) = - tmp4;
   B(1,4) =   tmp5;
   B(2,4) = - tmp6;
   B(3,4) = - A(3,4);
   B(0,5) = - A(2,5)*tmp1 + A(0,1)*A(1,5) - A(0,5) + A(3,5)*tmp2 + A(4,5)*tmp4;
   B(1,5) =   A(2,5)*A(1,2) - A(1,5) - A(3,5)*tmp3 - A(4,5)*tmp5;
   B(2,5) = - A(2,5) + A(3,5)*A(2,3) + A(4,5)*tmp6;
   B(3,5) = - A(3,5) + A(4,5)*A(3,4);
   B(4,5) = - A(4,5);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense \f$ 6 \times 6 \f$ matrix.
// \ingroup dense_matrix
//
// \param dm The diagonal dense matrix to be inverted.
// \return void
//
// This function inverts the given diagonal dense \f$ 6 \times 6 \f$ matrix via the rule of Sarrus.
// The matrix inversion fails if the given matrix is singular and not invertible. In this case a
// \a std::invalid_argument exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonal6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( (*dm).rows()    == 6UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*dm).columns() == 6UL, "Invalid number of columns detected" );

   using ET = ElementType_t<MT>;

   MT& A( *dm );

   const ET tmp1( A(0,0)*A(1,1) );
   const ET tmp2( A(3,3)*A(4,4) );
   const ET tmp3( tmp1*A(2,2) );
   const ET tmp4( tmp2*A(5,5) );
   const ET tmp5( A(0,0)*tmp4 );
   const ET tmp6( tmp3*A(3,3) );

   const ET det( tmp3*tmp4 );

   if( !isDivisor( det ) ) {
      BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
   }

   const ET idet( ET(1) / det );

   A(0,0) = A(1,1)*A(2,2)*tmp4*idet;
   A(1,1) = tmp5*A(2,2)*idet;
   A(2,2) = tmp1*tmp4*idet;
   A(3,3) = tmp3*A(4,4)*A(5,5)*idet;
   A(4,4) = tmp6*A(5,5)*idet;
   A(5,5) = tmp2*tmp3*idet;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense 6x6 matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense 6x6 matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert6x6( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertGeneral6x6  ( *dm ); break;
      case byLDLT     : invertSymmetric6x6( *dm ); break;
      case byLDLH     : invertHermitian6x6( *dm ); break;
      case byLLH      : invertHermitian6x6( *dm ); break;
      case asGeneral  : invertGeneral6x6  ( *dm ); break;
      case asSymmetric: invertSymmetric6x6( *dm ); break;
      case asHermitian: invertHermitian6x6( *dm ); break;
      case asLower    : invertLower6x6    ( *dm ); break;
      case asUniLower : invertUniLower6x6 ( *dm ); break;
      case asUpper    : invertUpper6x6    ( *dm ); break;
      case asUniUpper : invertUniUpper6x6 ( *dm ); break;
      case asDiagonal : invertDiagonal6x6 ( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS FOR NxN MATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place LU-based inversion of the given general dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given general dense matrix by means of an LU decomposition. The
// matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertByLU( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   const size_t n( min( (*dm).rows(), (*dm).columns() ) );
   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[n] );

   getrf( *dm, ipiv.get() );
   getri( *dm, ipiv.get() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place Bunch-Kaufman-based inversion of the given symmetric dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense matrix by means of a Bunch-Kaufman decomposition.
// The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertByLDLT( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_USER_ASSERT( isSymmetric( *dm ), "Invalid non-symmetric matrix detected" );

   const char uplo( ( SO )?( 'L' ):( 'U' ) );
   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[(*dm).rows()] );

   sytrf( *dm, uplo, ipiv.get() );
   sytri( *dm, uplo, ipiv.get() );

   if( SO ) {
      for( size_t i=1UL; i<(*dm).rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            (*dm)(j,i) = (*dm)(i,j);
         }
      }
   }
   else {
      for( size_t j=1UL; j<(*dm).columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            (*dm)(j,i) = (*dm)(i,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place Bunch-Kaufman-based inversion of the given symmetric dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given symmetric dense matrix by means of a Bunch-Kaufman decomposition.
// The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline auto invertByLDLH( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsBuiltin_v< ElementType_t<MT> > >
{
   invertByLDLT( *dm );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place Bunch-Kaufman-based inversion of the given Hermitian dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given Hermitian dense matrix by means of a Bunch-Kaufman decomposition.
// The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline auto invertByLDLH( DenseMatrix<MT,SO>& dm )
   -> EnableIf_t< IsComplex_v< ElementType_t<MT> > >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_USER_ASSERT( isHermitian( *dm ), "Invalid non-Hermitian matrix detected" );

   const char uplo( ( SO )?( 'L' ):( 'U' ) );
   const std::unique_ptr<blas_int_t[]> ipiv( new blas_int_t[(*dm).rows()] );

   hetrf( *dm, uplo, ipiv.get() );
   hetri( *dm, uplo, ipiv.get() );

   if( SO ) {
      for( size_t i=1UL; i<(*dm).rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            (*dm)(j,i) = conj( (*dm)(i,j) );
         }
      }
   }
   else {
      for( size_t j=1UL; j<(*dm).columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            (*dm)(j,i) = conj( (*dm)(i,j) );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place Cholesky-based inversion of the given symmetric PD dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense matrix by means of a Cholesky decomposition. The matrix
// inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertByLLH( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_USER_ASSERT( isHermitian( *dm ), "Invalid non-symmetric matrix detected" );

   const char uplo( ( SO )?( 'L' ):( 'U' ) );

   potrf( *dm, uplo );
   potri( *dm, uplo );

   if( SO ) {
      for( size_t i=1UL; i<(*dm).rows(); ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            (*dm)(j,i) = conj( (*dm)(i,j) );
         }
      }
   }
   else {
      for( size_t j=1UL; j<(*dm).columns(); ++j ) {
         for( size_t i=0UL; i<j; ++i ) {
            (*dm)(j,i) = conj( (*dm)(i,j) );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower triangular dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given lower triangular dense matrix by means of a forward
// substitution. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertLowerNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   trtri( *dm, 'L', 'N' );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given lower unitriangular dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given lower unitriangular dense matrix by means of a forward
// substitution. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniLowerNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   trtri( *dm, 'L', 'U' );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper triangular dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given upper triangular dense matrix by means of a back substitution.
// The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUpperNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   trtri( *dm, 'U', 'N' );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given upper unitriangular dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given upper unitriangular dense matrix by means of a back
// substitution. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertUniUpperNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   trtri( *dm, 'U', 'U' );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given diagonal dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given diagonal dense matrix by means of computing the reciprocal of
// each diagonal element. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invertDiagonalNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_ADAPTOR_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   for( size_t i=0UL; i<(*dm).rows(); ++i )
   {
      if( !isDivisor( (*dm)(i,i) ) ) {
         BLAZE_THROW_DIVISION_BY_ZERO( "Inversion of singular matrix failed" );
      }

      invert( (*dm)(i,i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place inversion of the given dense NxN matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense NxN matrix via the specified matrix inversion algorithm
// \a IF. The matrix inversion fails if the given matrix is singular and not invertible. In this
// case a \a std::runtime_error exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invertNxN( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   BLAZE_INTERNAL_ASSERT( isSquare( *dm ), "Non-square matrix detected" );

   switch( IF ) {
      case byLU       : invertByLU       ( *dm ); break;
      case byLDLT     : invertByLDLT     ( *dm ); break;
      case byLDLH     : invertByLDLH     ( *dm ); break;
      case byLLH      : invertByLLH      ( *dm ); break;
      case asGeneral  : invertByLU       ( *dm ); break;
      case asSymmetric: invertByLDLT     ( *dm ); break;
      case asHermitian: invertByLDLH     ( *dm ); break;
      case asLower    : invertLowerNxN   ( *dm ); break;
      case asUniLower : invertUniLowerNxN( *dm ); break;
      case asUpper    : invertUpperNxN   ( *dm ); break;
      case asUniUpper : invertUniUpperNxN( *dm ); break;
      case asDiagonal : invertDiagonalNxN( *dm ); break;
      default: BLAZE_INTERNAL_ASSERT( false, "Unhandled case detected" );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INVERSION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Inversion functions */
//@{
template< typename MT, bool SO >
inline void invert( DenseMatrix<MT,SO>& dm );

template< InversionFlag IF, typename MT, bool SO >
inline void invert( DenseMatrix<MT,SO>& dm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense square matrix. The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< typename MT  // Type of the dense matrix
        , bool SO >    // Storage order of the dense matrix
inline void invert( DenseMatrix<MT,SO>& dm )
{
   invert< getInversionFlag<MT>() >( *dm );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the given dense matrix.
// \ingroup dense_matrix
//
// \param dm The dense matrix to be inverted.
// \return void
// \exception std::invalid_argument Invalid non-square matrix provided.
// \exception std::runtime_error Inversion of singular matrix failed.
//
// This function inverts the given dense matrix by means of the specified matrix type or matrix
// inversion algorithm \c IF (see the InversionFlag documentation):

   \code
   invert<asLower>( A );     // Inversion of a lower triangular matrix
   invert<asUniUpper>( A );  // Inversion of an upper unitriangular matrix
   invert<byLU>( A );        // Inversion by means of an LU decomposition
   invert<byLLH>( A );       // Inversion by means of a Cholesky decomposition
   ...
   \endcode

// The matrix inversion fails if ...
//
//  - ... the given matrix is not a square matrix;
//  - ... the given matrix is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or an exception is thrown.
//
// \note The matrix inversion can only be used for dense matrices with \c float, \c double,
// \c complex<float> or \c complex<double> element type. The attempt to call the function with
// matrices of any other element type results in a compile time error!
//
// \note This function can only be used if a fitting LAPACK library is available and linked to
// the executable. Otherwise a linker error will be created.
//
// \note This function does only provide the basic exception safety guarantee, i.e. in case of an
// exception \a dm may already have been modified.
*/
template< InversionFlag IF  // Inversion algorithm
        , typename MT       // Type of the dense matrix
        , bool SO >         // Storage order of the dense matrix
inline void invert( DenseMatrix<MT,SO>& dm )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_STRICTLY_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNIFORM_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_BE_BLAS_COMPATIBLE_TYPE( ElementType_t<MT> );

   if( !isSquare( *dm ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid non-square matrix provided" );
   }

   switch( (*dm).rows() ) {
      case 0UL:                       break;
      case 1UL: invert( (*dm)(0,0) ); break;
      case 2UL: invert2x2<IF>( *dm ); break;
      case 3UL: invert3x3<IF>( *dm ); break;
      case 4UL: invert4x4<IF>( *dm ); break;
      case 5UL: invert5x5<IF>( *dm ); break;
      case 6UL: invert6x6<IF>( *dm ); break;
      default : invertNxN<IF>( *dm ); break;
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *dm ), "Broken invariant detected" );
}
//*************************************************************************************************

} // namespace blaze

#endif
