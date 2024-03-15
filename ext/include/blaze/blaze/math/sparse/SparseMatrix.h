//=================================================================================================
/*!
//  \file blaze/math/sparse/SparseMatrix.h
//  \brief Header file for utility functions for sparse matrices
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the foyllowing disclaimer.
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

#ifndef _BLAZE_MATH_SPARSE_SPARSEMATRIX_H_
#define _BLAZE_MATH_SPARSE_SPARSEMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Triangular.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/expressions/SparseMatrix.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/Equal.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsFinite.h>
#include <blaze/math/shims/IsInf.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/IsReal.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsIdentity.h>
#include <blaze/math/typetraits/IsInvertible.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsSquare.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniform.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniTriangular.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
#include <blaze/math/typetraits/UnderlyingScalar.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsBuiltin.h>
#include <blaze/util/typetraits/IsComplex.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseMatrix operators */
//@{
template< typename MT, bool SO, typename ST >
auto operator*=( SparseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator*=( SparseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator/=( SparseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;

template< typename MT, bool SO, typename ST >
auto operator/=( SparseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse matrix and
//        a scalar value (\f$ A*=s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side sparse matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   if( IsRestricted_v<MT> ) {
      if( !tryMult( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted matrix" );
      }
   }

   if( !IsResizable_v< ElementType_t<MT> > && isZero( scalar ) )
   {
      reset( *mat );
   }
   else
   {
      decltype(auto) left( derestrict( *mat ) );

      const size_t iend( SO == rowMajor ? (*mat).rows() : (*mat).columns() );
      for( size_t i=0UL; i<iend; ++i ) {
         const auto last( left.end(i) );
         for( auto element=left.begin(i); element!=last; ++element ) {
            element->value() *= scalar;
         }
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a temporary sparse matrix
//        and a scalar value (\f$ A*=s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side temporary sparse matrix for the multiplication.
// \param scalar The right-hand side scalar value for the multiplication.
// \return Reference to the left-hand side sparse matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator*=( SparseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator*=( *mat, scalar );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a sparse matrix by a scalar value
//        (\f$ A/=s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side sparse matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseMatrix<MT,SO>& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   BLAZE_USER_ASSERT( !isZero( scalar ), "Division by zero detected" );

   if( IsRestricted_v<MT> ) {
      if( !tryDiv( *mat, 0UL, 0UL, (*mat).rows(), (*mat).columns(), scalar ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid scaling of restricted matrix" );
      }
   }

   using ScalarType = If_t< IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > ||
                            IsFloatingPoint_v< UnderlyingBuiltin_t<ST> >
                          , If_t< IsComplex_v< UnderlyingScalar_t<MT> > && IsBuiltin_v<ST>
                                , DivTrait_t< UnderlyingBuiltin_t<MT>, ST >
                                , DivTrait_t< UnderlyingScalar_t<MT>, ST > >
                          , ST >;

   decltype(auto) left( derestrict( *mat ) );

   if( IsInvertible_v<ScalarType> ) {
      const ScalarType tmp( ScalarType(1)/static_cast<ScalarType>( scalar ) );
      const size_t iend( SO == rowMajor ? (*mat).rows() : (*mat).columns() );
      for( size_t i=0UL; i<iend; ++i ) {
         const auto last( left.end(i) );
         for( auto element=left.begin(i); element!=last; ++element ) {
            element->value() *= tmp;
         }
      }
   }
   else {
      const size_t iend( SO == rowMajor ? (*mat).rows() : (*mat).columns() );
      for( size_t i=0UL; i<iend; ++i ) {
         const auto last( left.end(i) );
         for( auto element=left.begin(i); element!=last; ++element ) {
            element->value() /= scalar;
         }
      }
   }

   BLAZE_INTERNAL_ASSERT( isIntact( *mat ), "Invariant violation detected" );

   return *mat;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a temporary sparse matrix by a scalar
//        value (\f$ A/=s \f$).
// \ingroup sparse_matrix
//
// \param mat The left-hand side temporary sparse matrix for the division.
// \param scalar The right-hand side scalar value for the division.
// \return Reference to the left-hand side sparse matrix.
// \exception std::invalid_argument Invalid scaling of restricted matrix.
//
// In case the matrix \a MT is restricted and the assignment would violate an invariant of the
// matrix, a \a std::invalid_argument exception is thrown.
//
// \note A division by zero is only checked by an user assert.
*/
template< typename MT    // Type of the left-hand side sparse matrix
        , bool SO        // Storage order
        , typename ST >  // Data type of the right-hand side scalar
inline auto operator/=( SparseMatrix<MT,SO>&& mat, ST scalar )
   -> EnableIf_t< IsScalar_v<ST>, MT& >
{
   return operator/=( *mat, scalar );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name SparseMatrix functions */
//@{
template< typename MT, bool SO >
bool isnan( const SparseMatrix<MT,SO>& sm );

template< typename MT, bool SO >
bool isinf( const SparseMatrix<MT,SO>& sm );

template< typename MT, bool SO >
bool isfinite( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isSymmetric( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isHermitian( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniform( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isZero( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isLower( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniLower( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isStrictlyLower( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUpper( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isUniUpper( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isStrictlyUpper( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isDiagonal( const SparseMatrix<MT,SO>& sm );

template< RelaxationFlag RF, typename MT, bool SO >
bool isIdentity( const SparseMatrix<MT,SO>& sm );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse matrix for not-a-number elements.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked for not-a-number elements.
// \return \a true if at least one element of the sparse matrix is not-a-number, \a false otherwise.
//
// This function checks the sparse matrix for not-a-number (NaN) elements. If at least one
// element of the matrix is not-a-number, the function returns \a true, otherwise it returns
// \a false.

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isnan( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isnan( const SparseMatrix<MT,SO>& sm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return false;

   CompositeType_t<MT> A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            if( isnan( element->value() ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
            if( isnan( element->value() ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse matrix for infinite elements.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked for infinite elements.
// \return \a true if at least one element of the sparse matrix is infinite, \a false otherwise.
//
// This function checks the sparse matrix for infinite (inf) elements. If at least one
// element of the matrix is infinite, the function returns \a true, otherwise it returns
// \a false.

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isinf( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isinf( const SparseMatrix<MT,SO>& sm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return false;

   CompositeType_t<MT> A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            if( isinf( element->value() ) ) return true;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
            if( isinf( element->value() ) ) return true;
      }
   }

   return false;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks the given sparse matrix for finite elements.
// \ingroup sparse_matrix
//
// \param sm The sparsre matrix to be checked for finite elements.
// \return \a true if all elements of the matrix are finite, \a false otherwise.
//
// This function checks if all elements of the sparse matrix are finite elements (i.e. normal,
// subnormal or zero elements, but not infinite or NaN). If all elements of the matrix are
// finite, the function returns \a true, otherwise it returns \a false.

   \code
   blaze::CompressedMatrix<double> A( 3UL, 4UL );
   // ... Initialization
   if( isfinite( A ) ) { ... }
   \endcode
*/
template< typename MT  // Type of the sparse matrix
        , bool SO >    // Storage order
bool isfinite( const SparseMatrix<MT,SO>& sm )
{
   if( !IsFloatingPoint_v< UnderlyingBuiltin_t<MT> > )
      return true;

   CompositeType_t<MT> A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            if( !isfinite( element->value() ) ) return false;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
            if( !isfinite( element->value() ) ) return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is symmetric.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is symmetric, \a false if not.
//
// This function checks if the given sparse matrix is symmetric. The matrix is considered to be
// symmetric if it is a square matrix whose transpose is equal to itself (\f$ A = A^T \f$). The
// following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isSymmetric( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isSymmetric<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a symmetric matrix:

   \code
   if( isSymmetric( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isSymmetric( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsSymmetric_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsUniform_v<MT> || (*sm).rows() < 2UL )
      return true;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            const size_t j( element->index() );

            if( i == j || isDefault<RF>( element->value() ) )
               continue;

            const auto pos( A.find( j, i ) );
            if( pos == A.end(j) || !equal<RF>( pos->value(), element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            const size_t i( element->index() );

            if( j == i || isDefault<RF>( element->value() ) )
               continue;

            const auto pos( A.find( j, i ) );
            if( pos == A.end(i) || !equal<RF>( pos->value(), element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is Hermitian.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is Hermitian, \a false if not.
//
// This function checks if the given sparse matrix is an Hermitian matrix. The matrix is considered
// to be an Hermitian matrix if it is a square matrix whose conjugate transpose is equal to itself
// (\f$ A = \overline{A^T} \f$), i.e. each matrix element \f$ a_{ij} \f$ is equal to the complex
// conjugate of the element \f$ a_{ji} \f$. The following code example demonstrates the use of the
// function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isHermitian( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isHermitian<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an Hermitian matrix:

   \code
   if( isHermitian( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isHermitian( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using ET  = ElementType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsHermitian_v<MT> )
      return true;

   if( !IsScalar_v<ET> || !isSquare( *sm ) )
      return false;

   if( IsBuiltin_v<ET> && IsUniform_v<MT> )
      return true;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            const size_t j( element->index() );

            if( isDefault<RF>( element->value() ) )
               continue;

            if( i == j && !isReal<RF>( element->value() ) )
               return false;

            const auto pos( A.find( j, i ) );
            if( pos == A.end(j) || !equal<RF>( pos->value(), conj( element->value() ) ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            const size_t i( element->index() );

            if( isDefault<RF>( element->value() ) )
               continue;

            if( j == i && !isReal<RF>( element->value() ) )
               return false;

            const auto pos( A.find( j, i ) );
            if( pos == A.end(i) || !equal<RF>( pos->value(), conj( element->value() ) ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given row-major triangular sparse matrix is a uniform matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the sparse matrix
bool isUniform_backend( const SparseMatrix<MT,false>& sm, TrueType )
{
   BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*sm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*sm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t ibegin( ( IsStrictlyLower_v<MT> )?( 1UL ):( 0UL ) );
   const size_t iend  ( ( IsStrictlyUpper_v<MT> )?( (*sm).rows()-1UL ):( (*sm).rows() ) );

   for( size_t i=ibegin; i<iend; ++i ) {
      for( auto element=(*sm).begin(i); element!=(*sm).end(i); ++element ) {
         if( !isDefault<RF>( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given column-major triangular sparse matrix is a uniform matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the sparse matrix
bool isUniform_backend( const SparseMatrix<MT,true>& sm, TrueType )
{
   BLAZE_CONSTRAINT_MUST_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*sm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*sm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t jbegin( ( IsStrictlyUpper_v<MT> )?( 1UL ):( 0UL ) );
   const size_t jend  ( ( IsStrictlyLower_v<MT> )?( (*sm).columns()-1UL ):( (*sm).columns() ) );

   for( size_t j=jbegin; j<jend; ++j ) {
      for( auto element=(*sm).begin(j); element!=(*sm).end(j); ++element ) {
         if( !isDefault<RF>( element->value() ) )
            return false;
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given row-major general sparse matrix is a uniform matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the sparse matrix
bool isUniform_backend( const SparseMatrix<MT,false>& sm, FalseType )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*sm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*sm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t maxElements( (*sm).rows() * (*sm).columns() );

   if( (*sm).nonZeros() != maxElements )
   {
      for( size_t i=0UL; i<(*sm).rows(); ++i ) {
         for( auto element=(*sm).begin(i); element!=(*sm).end(i); ++element ) {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else
   {
      BLAZE_INTERNAL_ASSERT( (*sm).find(0UL,0UL) != (*sm).end(0UL), "Missing element detected" );

      const auto& cmp( (*sm)(0UL,0UL) );

      for( size_t i=0UL; i<(*sm).rows(); ++i ) {
         for( auto element=(*sm).begin(i); element!=(*sm).end(i); ++element ) {
            if( !equal<RF>( element->value(), cmp ) )
               return false;
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checks if the given column-major general sparse matrix is a uniform matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT >      // Type of the sparse matrix
bool isUniform_backend( const SparseMatrix<MT,true>& sm, FalseType )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRIANGULAR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MT );

   BLAZE_INTERNAL_ASSERT( (*sm).rows()    != 0UL, "Invalid number of rows detected"    );
   BLAZE_INTERNAL_ASSERT( (*sm).columns() != 0UL, "Invalid number of columns detected" );

   const size_t maxElements( (*sm).rows() * (*sm).columns() );

   if( (*sm).nonZeros() != maxElements )
   {
      for( size_t j=0UL; j<(*sm).columns(); ++j ) {
         for( auto element=(*sm).begin(j); element!=(*sm).end(j); ++element ) {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else
   {
      BLAZE_INTERNAL_ASSERT( (*sm).find(0UL,0UL) != (*sm).end(0UL), "Missing element detected" );

      const auto& cmp( (*sm)(0UL,0UL) );

      for( size_t j=0UL; j<(*sm).columns(); ++j ) {
         for( auto element=(*sm).begin(j); element!=(*sm).end(j); ++element ) {
            if( !equal<RF>( element->value(), cmp ) )
               return false;
         }
      }
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a uniform matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a uniform matrix, \a false if not.
//
// This function checks if the given sparse matrix is a uniform matrix. The matrix is considered
// to be uniform if all its elements are identical. The following code example demonstrates the
// use of the function:

   \code
   blaze::DynamicMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniform( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniform<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a uniform matrix:

   \code
   if( isUniform( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isUniform( const SparseMatrix<MT,SO>& sm )
{
   if( IsUniform_v<MT> ||
       (*sm).rows() == 0UL || (*sm).columns() == 0UL ||
       ( (*sm).rows() == 1UL && (*sm).columns() == 1UL ) )
      return true;

   if( IsUniTriangular_v<MT> )
      return false;

   CompositeType_t<MT> A( *sm );  // Evaluation of the sparse matrix operand

   return isUniform_backend<RF>( A, typename IsTriangular<MT>::Type() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a zero matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a zero matrix, \a false if not.
//
// This function checks if the given sparse matrix is a zero matrix. The matrix is considered to
// be zero if all its elements are zero. The following code example demonstrates the use of the
// function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isZero( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isZero<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results is a zero matrix:

   \code
   if( isZero( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isZero( const SparseMatrix<MT,SO>& sm )
{
   const size_t M( (*sm).rows()    );
   const size_t N( (*sm).columns() );

   if( IsZero_v<MT> || M == 0UL || N == 0UL )
      return true;

   if( IsUniTriangular_v<MT> )
      return false;

   CompositeType_t<MT> A( *sm );  // Evaluation of the sparse matrix operand

   const size_t iend( SO == rowMajor ? A.rows() : A.columns() );

   for( size_t i=0UL; i<iend; ++i ) {
      for( auto element=A.begin(i); element!=A.end(i); ++element ) {
         if( !isZero<RF>( element->value() ) ) {
            return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a lower triangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a lower triangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is a lower triangular matrix. The matrix is
// considered to be lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        l_{0,0} & 0       & 0       & \cdots & 0       \\
                        l_{1,0} & l_{1,1} & 0       & \cdots & 0       \\
                        l_{2,0} & l_{2,1} & l_{2,2} & \cdots & 0       \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & l_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially lower triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a lower triangular matrix:

   \code
   if( isLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isLower( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsLower_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsZero_v<MT> || (*sm).rows() < 2UL )
      return true;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows()-1UL; ++i ) {
         for( auto element=A.lowerBound(i,i+1UL); element!=A.end(i); ++element )
         {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=1UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            if( element->index() >= j )
               break;

            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a lower unitriangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a lower unitriangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is a lower unitriangular matrix. The matrix is
// considered to be lower unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 1       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 1       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 1      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a lower unitriangular matrix:

   \code
   if( isUniLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isUniLower( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUniLower_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i )
      {
         auto element( A.lowerBound(i,i) );

         if( element == A.end(i) || element->index() != i || !isOne<RF>( element->value() ) )
            return false;

         ++element;

         for( ; element!=A.end(i); ++element ) {
            if( !isZero<RF>( element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j )
      {
         bool hasDiagonalElement( false );

         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            if( element->index() >= j ) {
               if( element->index() != j || !isOne<RF>( element->value() ) )
                  return false;
               hasDiagonalElement = true;
               break;
            }

            if( !isZero<RF>( element->value() ) )
               return false;
         }

         if( !hasDiagonalElement ) {
            return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a strictly lower triangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a strictly lower triangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is a strictly lower triangular matrix. The
// matrix is considered to be strictly lower triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0       & 0       & 0       & \cdots & 0      \\
                        l_{1,0} & 0       & 0       & \cdots & 0      \\
                        l_{2,0} & l_{2,1} & 0       & \cdots & 0      \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots \\
                        l_{N,0} & l_{N,1} & l_{N,2} & \cdots & 0      \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyLower( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyLower<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a strictly lower triangular
// matrix:

   \code
   if( isStrictlyLower( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isStrictlyLower( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsStrictlyLower_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsZero_v<MT> || (*sm).rows() < 2UL )
      return true;

   if( IsUniLower_v<MT> || IsUniUpper_v<MT> )
      return false;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.lowerBound(i,i); element!=A.end(i); ++element )
         {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            if( element->index() > j )
               break;

            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is an upper triangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is an upper triangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is an upper triangular matrix. The matrix is
// considered to be upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        u_{0,0} & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0       & u_{1,1} & u_{1,2} & \cdots & u_{1,N} \\
                        0       & 0       & u_{2,2} & \cdots & u_{2,N} \\
                        \vdots  & \vdots  & \vdots  & \ddots & \vdots  \\
                        0       & 0       & 0       & \cdots & u_{N,N} \\
                        \end{array}\right).\f]

// \f$ 0 \times 0 \f$ or \f$ 1 \times 1 \f$ matrices are considered as trivially upper triangular.
// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper triangular matrix:

   \code
   if( isUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isUpper( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUpper_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsZero_v<MT> || (*sm).rows() < 2UL )
      return true;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=1UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            if( element->index() >= i )
               break;

            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns()-1UL; ++j ) {
         for( auto element=A.lowerBound(j+1UL,j); element!=A.end(j); ++element )
         {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is an upper unitriangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is an upper unitriangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is an upper unitriangular matrix. The matrix is
// considered to be upper unitriangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 1       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 1       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 1       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isUniUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an upper unitriangular matrix:

   \code
   if( isUniUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isUniUpper( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsUniUpper_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i )
      {
         bool hasDiagonalElement( false );

         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            if( element->index() >= i ) {
               if( element->index() != i || !isOne<RF>( element->value() ) )
                  return false;
               hasDiagonalElement = true;
               break;
            }
            else if( !isZero<RF>( element->value() ) ) {
               return false;
            }
         }

         if( !hasDiagonalElement ) {
            return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j )
      {
         auto element( A.lowerBound(j,j) );

         if( element == A.end(j) || element->index() != j || !isOne<RF>( element->value() ) )
            return false;

         ++element;

         for( ; element!=A.end(j); ++element ) {
            if( !isZero<RF>( element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given sparse matrix is a strictly upper triangular matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is a strictly upper triangular matrix, \a false if not.
//
// This function checks if the given sparse matrix is a strictly upper triangular matrix. The
// matrix is considered to be strictly upper triangular if it is a square matrix of the form

                        \f[\left(\begin{array}{*{5}{c}}
                        0      & u_{0,1} & u_{0,2} & \cdots & u_{0,N} \\
                        0      & 0       & u_{1,2} & \cdots & u_{1,N} \\
                        0      & 0       & 0       & \cdots & u_{2,N} \\
                        \vdots & \vdots  & \vdots  & \ddots & \vdots  \\
                        0      & 0       & 0       & \cdots & 0       \\
                        \end{array}\right).\f]

// The following code example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isStrictlyUpper( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isStrictlyUpper<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a strictly upper triangular
// matrix:

   \code
   if( isStrictlyUpper( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isStrictlyUpper( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsStrictlyUpper_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsZero_v<MT> || (*sm).rows() < 2UL )
      return true;

   if( IsUniLower_v<MT> || IsUniUpper_v<MT> )
      return false;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            if( element->index() > i )
               break;

            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.lowerBound(j,j); element!=A.end(j); ++element )
         {
            if( !isDefault<RF>( element->value() ) )
               return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give sparse matrix is diagonal.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is diagonal, \a false if not.
//
// This function tests whether the matrix is diagonal, i.e. if the non-diagonal elements are
// default elements. In case of integral or floating point data types, a diagonal matrix has
// the form

                        \f[\left(\begin{array}{*{5}{c}}
                        aa     & 0      & 0      & \cdots & 0  \\
                        0      & bb     & 0      & \cdots & 0  \\
                        0      & 0      & cc     & \cdots & 0  \\
                        \vdots & \vdots & \vdots & \ddots & 0  \\
                        0      & 0      & 0      & 0      & xx \\
                        \end{array}\right)\f]

// The following example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isDiagonal( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDiagonal<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in a diagonal matrix:

   \code
   if( isDiagonal( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isDiagonal( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsDiagonal_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   if( IsZero_v<MT> || (*sm).rows() < 2UL )
      return true;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i ) {
         for( auto element=A.begin(i); element!=A.end(i); ++element )
            if( element->index() != i && !isDefault<RF>( element->value() ) )
               return false;
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j ) {
         for( auto element=A.begin(j); element!=A.end(j); ++element )
            if( element->index() != j && !isDefault<RF>( element->value() ) )
               return false;
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the give sparse matrix is an identity matrix.
// \ingroup sparse_matrix
//
// \param sm The sparse matrix to be checked.
// \return \a true if the matrix is an identity matrix, \a false if not.
//
// This function tests whether the matrix is an identity matrix, i.e. if the diagonal elements
// are 1 and the non-diagonal elements are 0. In case of integral or floating point data types,
// an identity matrix has the form

                        \f[\left(\begin{array}{*{5}{c}}
                        1      & 0      & 0      & \cdots & 0 \\
                        0      & 1      & 0      & \cdots & 0 \\
                        0      & 0      & 1      & \cdots & 0 \\
                        \vdots & \vdots & \vdots & \ddots & 0 \\
                        0      & 0      & 0      & 0      & 1 \\
                        \end{array}\right)\f]

// The following example demonstrates the use of the function:

   \code
   blaze::CompressedMatrix<int,blaze::rowMajor> A, B;
   // ... Initialization
   if( isIdentity( A ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isIdentity<relaxed>( A ) ) { ... }
   \endcode

// It is also possible to check if a matrix expression results in an identity matrix:

   \code
   if( isIdentity( A * B ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary matrix.
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename MT        // Type of the sparse matrix
        , bool SO >          // Storage order
bool isIdentity( const SparseMatrix<MT,SO>& sm )
{
   using RT  = ResultType_t<MT>;
   using RN  = ReturnType_t<MT>;
   using CT  = CompositeType_t<MT>;
   using Tmp = If_t< IsExpression_v<RN>, const RT, CT >;

   if( IsIdentity_v<MT> )
      return true;

   if( !isSquare( *sm ) )
      return false;

   Tmp A( *sm );  // Evaluation of the sparse matrix operand

   if( SO == rowMajor ) {
      for( size_t i=0UL; i<A.rows(); ++i )
      {
         bool hasDiagonalElement( false );

         for( auto element=A.begin(i); element!=A.end(i); ++element )
         {
            if( element->index() == i ) {
               if( !isOne<RF>( element->value() ) )
                  return false;
               hasDiagonalElement = true;
            }
            else if( !isZero<RF>( element->value() ) ) {
               return false;
            }
         }

         if( !hasDiagonalElement ) {
            return false;
         }
      }
   }
   else {
      for( size_t j=0UL; j<A.columns(); ++j )
      {
         bool hasDiagonalElement( false );

         for( auto element=A.begin(j); element!=A.end(j); ++element )
         {
            if( element->index() == j ) {
               if( !isOne<RF>( element->value() ) )
                  return false;
               hasDiagonalElement = true;
            }
            else if( !isZero<RF>( element->value() ) ) {
               return false;
            }
         }

         if( !hasDiagonalElement ) {
            return false;
         }
      }
   }

   return true;
}
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a single element, a range or selection of elements from the given sparse matrix.
// \ingroup sparse_matrix
//
// \param sm The given sparse matrix.
// \param args The runtime arguments for the erase call.
// \return The result of the according erase member function.
//
// This function represents an abstract interface for erasing a single element, a range of
// elements or a selection of elements from the given sparse matrix. It forwards the given
// arguments to the according \a erase() member function of the sparse matrix and returns
// the result of the function call.
*/
template< typename MT         // Type of the sparse matrix
        , bool SO             // Storage order
        , typename... Args >  // Type of the erase arguments
auto erase( SparseMatrix<MT,SO>& sm, Args&&... args )
   -> decltype( (*sm).erase( std::forward<Args>( args )... ) )
{
   return (*sm).erase( std::forward<Args>( args )... );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
