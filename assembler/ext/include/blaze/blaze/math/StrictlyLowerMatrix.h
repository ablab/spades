//=================================================================================================
/*!
//  \file blaze/math/StrictlyLowerMatrix.h
//  \brief Header file for the complete StrictlyLowerMatrix implementation
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

#ifndef _BLAZE_MATH_STRICTLYLOWERMATRIX_H_
#define _BLAZE_MATH_STRICTLYLOWERMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cmath>
#include <vector>
#include <blaze/math/Aliases.h>
#include <blaze/math/adaptors/LowerMatrix.h>
#include <blaze/math/adaptors/StrictlyLowerMatrix.h>
#include <blaze/math/adaptors/UniLowerMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Resizable.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/DenseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/SparseMatrix.h>
#include <blaze/math/StrictlyUpperMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/ZeroMatrix.h>
#include <blaze/util/Indices.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Random.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for StrictlyLowerMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of StrictlyLowerMatrix.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
class Rand< StrictlyLowerMatrix<MT,SO,DF> >
{
 public:
   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \return The generated random matrix.
   */
   inline const StrictlyLowerMatrix<MT,SO,DF> generate() const
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_RESIZABLE_TYPE( MT );

      StrictlyLowerMatrix<MT,SO,DF> matrix;
      randomize( matrix );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \return The generated random matrix.
   */
   inline const StrictlyLowerMatrix<MT,SO,DF> generate( size_t n ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE( MT );

      StrictlyLowerMatrix<MT,SO,DF> matrix( n );
      randomize( matrix );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \return The generated random matrix.
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   inline const StrictlyLowerMatrix<MT,SO,DF> generate( size_t n, size_t nonzeros ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE    ( MT );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      StrictlyLowerMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, nonzeros );

      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   */
   template< typename Arg >  // Min/max argument type
   inline const StrictlyLowerMatrix<MT,SO,DF> generate( const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_RESIZABLE_TYPE( MT );

      StrictlyLowerMatrix<MT,SO,DF> matrix;
      randomize( matrix, min, max );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   */
   template< typename Arg >  // Min/max argument type
   inline const StrictlyLowerMatrix<MT,SO,DF>
      generate( size_t n, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE( MT );

      StrictlyLowerMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, min, max );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random StrictlyLowerMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename Arg >  // Min/max argument type
   inline const StrictlyLowerMatrix<MT,SO,DF>
      generate( size_t n, size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE    ( MT );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      StrictlyLowerMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, nonzeros, min, max );

      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix ) const
   {
      randomize( matrix, typename IsDenseMatrix<MT>::Type() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a row-major sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   inline void randomize( StrictlyLowerMatrix<MT,false,DF>& matrix, size_t nonzeros ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL || n == 1UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );
      matrix.finalize( 0UL );

      std::vector<size_t> dist( n );

      for( size_t nz=0UL; nz<nonzeros; ) {
         const size_t index = rand<size_t>( 1UL, n-1UL );
         if( dist[index] == index ) continue;
         ++dist[index];
         ++nz;
      }

      for( size_t i=1UL; i<n; ++i ) {
         const Indices<size_t> indices( 0UL, i-1UL, dist[i] );
         for( size_t j : indices ) {
            matrix.append( i, j, rand<ET>() );
         }
         matrix.finalize( i );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a column-major sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   inline void randomize( StrictlyLowerMatrix<MT,true,DF>& matrix, size_t nonzeros ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL || n == 1UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );

      std::vector<size_t> dist( n-1UL );

      for( size_t nz=0UL; nz<nonzeros; ) {
         const size_t index = rand<size_t>( 0UL, n-2UL );
         if( dist[index] == n - index - 1UL ) continue;
         ++dist[index];
         ++nz;
      }

      for( size_t j=0UL; j<n-1UL; ++j ) {
         const Indices<size_t> indices( j+1UL, n-1UL, dist[j] );
         for( size_t i : indices ) {
            matrix.append( i, j, rand<ET>() );
         }
         matrix.finalize( j );
      }

      matrix.finalize( n-1UL );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix, const Arg& min, const Arg& max ) const
   {
      randomize( matrix, min, max, typename IsDenseMatrix<MT>::Type() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a row-major sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( StrictlyLowerMatrix<MT,false,DF>& matrix,
                          size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL || n == 1UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );
      matrix.finalize( 0UL );

      std::vector<size_t> dist( n );

      for( size_t nz=0UL; nz<nonzeros; ) {
         const size_t index = rand<size_t>( 1UL, n-1UL );
         if( dist[index] == index ) continue;
         ++dist[index];
         ++nz;
      }

      for( size_t i=1UL; i<n; ++i ) {
         const Indices<size_t> indices( 0UL, i-1UL, dist[i] );
         for( size_t j : indices ) {
            matrix.append( i, j, rand<ET>( min, max ) );
         }
         matrix.finalize( i );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a column-major sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( StrictlyLowerMatrix<MT,true,DF>& matrix,
                          size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > StrictlyLowerMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL || n == 1UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );

      std::vector<size_t> dist( n-1UL );

      for( size_t nz=0UL; nz<nonzeros; ) {
         const size_t index = rand<size_t>( 0UL, n-2UL );
         if( dist[index] == n - index - 1UL ) continue;
         ++dist[index];
         ++nz;
      }

      for( size_t j=0UL; j<n-1UL; ++j ) {
         const Indices<size_t> indices( j+1UL, n-1UL, dist[j] );
         for( size_t i : indices ) {
            matrix.append( i, j, rand<ET>( min, max ) );
         }
         matrix.finalize( j );
      }

      matrix.finalize( n-1UL );
   }
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*!\brief Randomization of a dense StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix, TrueType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      for( size_t i=1UL; i<n; ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            matrix(i,j) = rand<ET>();
         }
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix, FalseType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      const size_t n( matrix.rows() );

      if( n == 0UL || n == 1UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.2*n*n ) ) );

      randomize( matrix, nonzeros );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a dense StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix,
                          const Arg& min, const Arg& max, TrueType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      for( size_t i=1UL; i<n; ++i ) {
         for( size_t j=0UL; j<i; ++j ) {
            matrix(i,j) = rand<ET>( min, max );
         }
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse StrictlyLowerMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( StrictlyLowerMatrix<MT,SO,DF>& matrix,
                          const Arg& min, const Arg& max, FalseType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      const size_t n( matrix.rows() );

      if( n == 0UL || n == 1UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.2*n*n ) ) );

      randomize( matrix, nonzeros, min, max );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
