//=================================================================================================
/*!
//  \file blaze/math/DiagonalMatrix.h
//  \brief Header file for the complete DiagonalMatrix implementation
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

#ifndef _BLAZE_MATH_DIAGONALMATRIX_H_
#define _BLAZE_MATH_DIAGONALMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/adaptors/DiagonalMatrix.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Resizable.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/DenseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/SparseMatrix.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/UnderlyingBuiltin.h>
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
/*!\brief Specialization of the Rand class template for DiagonalMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of DiagonalMatrix.
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
class Rand< DiagonalMatrix<MT,SO,DF> >
{
 public:
   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \return The generated random matrix.
   */
   inline const DiagonalMatrix<MT,SO,DF> generate() const
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_RESIZABLE_TYPE( MT );

      DiagonalMatrix<MT,SO,DF> matrix;
      randomize( matrix );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \return The generated random matrix.
   */
   inline const DiagonalMatrix<MT,SO,DF> generate( size_t n ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE( MT );

      DiagonalMatrix<MT,SO,DF> matrix( n );
      randomize( matrix );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \return The generated random matrix.
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   inline const DiagonalMatrix<MT,SO,DF> generate( size_t n, size_t nonzeros ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE    ( MT );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      if( nonzeros > n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      DiagonalMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, nonzeros );

      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   */
   template< typename Arg >  // Min/max argument type
   inline const DiagonalMatrix<MT,SO,DF> generate( const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_NOT_BE_RESIZABLE_TYPE( MT );

      DiagonalMatrix<MT,SO,DF> matrix;
      randomize( matrix, min, max );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   */
   template< typename Arg >  // Min/max argument type
   inline const DiagonalMatrix<MT,SO,DF>
      generate( size_t n, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE( MT );

      DiagonalMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, min, max );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random DiagonalMatrix.
   //
   // \param n The number of rows and columns of the random matrix.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename Arg >  // Min/max argument type
   inline const DiagonalMatrix<MT,SO,DF>
      generate( size_t n, size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_RESIZABLE_TYPE    ( MT );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      if( nonzeros > DiagonalMatrix<MT,SO,DF>::maxNonZeros( n ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      DiagonalMatrix<MT,SO,DF> matrix( n );
      randomize( matrix, nonzeros, min, max );

      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix ) const
   {
      randomize( matrix, typename IsDenseMatrix<MT>::Type() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix, size_t nonzeros ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );

      Indices<size_t> indices( 0UL, n-1UL, nonzeros );
      size_t i( 0UL );

      for( size_t index : indices ) {
         for( ; i<index; ++i )
            matrix.finalize( i );
         matrix.append( i, i, rand<ET>() );
      }

      for( ; i<n; ++i ) {
         matrix.finalize( i );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix,
                          const Arg& min, const Arg& max ) const
   {
      randomize( matrix, min, max, typename IsDenseMatrix<MT>::Type() );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param nonzeros The number of non-zero elements of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix,
                          size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      if( nonzeros > n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( n == 0UL ) return;

      matrix.reset();
      matrix.reserve( nonzeros );

      Indices<size_t> indices( 0UL, n-1UL, nonzeros );
      size_t i( 0UL );

      for( size_t index : indices ) {
         for( ; i<index; ++i )
            matrix.finalize( i );
         matrix.append( i, i, rand<ET>( min, max ) );
      }

      for( ; i<n; ++i ) {
         matrix.finalize( i );
      }
   }
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   /*!\brief Randomization of a dense DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix, TrueType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      for( size_t i=0UL; i<n; ++i ) {
         matrix(i,i) = rand<ET>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix, FalseType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      const size_t n( matrix.rows() );

      if( n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, n ) );

      randomize( matrix, nonzeros );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a dense DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix,
                          const Arg& min, const Arg& max, TrueType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( MT );

      using ET = ElementType_t<MT>;

      const size_t n( matrix.rows() );

      for( size_t i=0UL; i<n; ++i ) {
         matrix(i,i) = rand<ET>( min, max );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse DiagonalMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( DiagonalMatrix<MT,SO,DF>& matrix,
                          const Arg& min, const Arg& max, FalseType ) const
   {
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( MT );

      const size_t n( matrix.rows() );

      if( n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, n ) );

      randomize( matrix, nonzeros, min, max );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAKE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setup of a random symmetric DiagonalMatrix.
//
// \param matrix The matrix to be randomized.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
void makeSymmetric( DiagonalMatrix<MT,SO,DF>& matrix )
{
   const size_t n( matrix.rows() );

   reset( matrix );

   for( size_t i=0UL; i<n; ++i ) {
      matrix(i,i) = rand< ElementType_t<MT> >();
   }

   BLAZE_INTERNAL_ASSERT( isSymmetric( matrix ), "Non-symmetric matrix detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setup of a random symmetric DiagonalMatrix.
//
// \param matrix The matrix to be randomized.
// \param min The smallest possible value for a matrix element.
// \param max The largest possible value for a matrix element.
// \return void
*/
template< typename MT     // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename Arg >  // Min/max argument type
void makeSymmetric( DiagonalMatrix<MT,SO,DF>& matrix, const Arg& min, const Arg& max )
{
   using Type = ElementType_t<MT>;

   const size_t n( matrix.rows() );

   reset( matrix );

   for( size_t i=0UL; i<n; ++i ) {
      matrix(i,i) = rand<Type>( min, max );
   }

   BLAZE_INTERNAL_ASSERT( isSymmetric( matrix ), "Non-symmetric matrix detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setup of a random Hermitian DiagonalMatrix.
//
// \param matrix The matrix to be randomized.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
void makeHermitian( DiagonalMatrix<MT,SO,DF>& matrix )
{
   using Type = UnderlyingBuiltin_t< ElementType_t<MT> >;

   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( Type );

   const size_t n( matrix.rows() );

   reset( matrix );

   for( size_t i=0UL; i<n; ++i ) {
      matrix(i,i) = rand<Type>();
   }

   BLAZE_INTERNAL_ASSERT( isHermitian( matrix ), "Non-Hermitian matrix detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setup of a random Hermitian DiagonalMatrix.
//
// \param matrix The matrix to be randomized.
// \param min The smallest possible value for a matrix element.
// \param max The largest possible value for a matrix element.
// \return void
*/
template< typename MT     // Type of the adapted matrix
        , bool SO         // Storage order of the adapted matrix
        , bool DF         // Density flag
        , typename Arg >  // Min/max argument type
void makeHermitian( DiagonalMatrix<MT,SO,DF>& matrix, const Arg& min, const Arg& max )
{
   using Type = UnderlyingBuiltin_t< ElementType_t<MT> >;

   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE( Type );

   const size_t n( matrix.rows() );

   reset( matrix );

   for( size_t i=0UL; i<n; ++i ) {
      matrix(i,i) = rand<Type>( min, max );
   }

   BLAZE_INTERNAL_ASSERT( isHermitian( matrix ), "Non-Hermitian matrix detected" );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setup of a random (Hermitian) positive definite DiagonalMatrix.
//
// \param matrix The matrix to be randomized.
// \return void
*/
template< typename MT  // Type of the adapted matrix
        , bool SO      // Storage order of the adapted matrix
        , bool DF >    // Density flag
void makePositiveDefinite( DiagonalMatrix<MT,SO,DF>& matrix )
{
   makeHermitian( matrix );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
