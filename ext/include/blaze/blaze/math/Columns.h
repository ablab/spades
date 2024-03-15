//=================================================================================================
/*!
//  \file blaze/math/Columns.h
//  \brief Header file for the complete Columns implementation
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

#ifndef _BLAZE_MATH_COLUMNS_H_
#define _BLAZE_MATH_COLUMNS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Columns.h>
#include <blaze/math/constraints/SparseMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/smp/DenseMatrix.h>
#include <blaze/math/smp/SparseMatrix.h>
#include <blaze/math/views/Columns.h>
#include <blaze/math/views/Rows.h>
#include <blaze/util/Random.h>
#include <blaze/util/typetraits/RemoveReference.h>


namespace blaze {

//=================================================================================================
//
//  RAND SPECIALIZATION FOR DENSE COLUMN SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for dense column selections.
// \ingroup random
//
// This specialization of the Rand class randomizes dense column selections.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
class Rand< Columns<MT,SO,true,SF,CCAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a dense column selection.
   //
   // \param columns The column selection to be randomized.
   // \return void
   */
   template< typename CT >  // Type of the column selection
   inline void randomize( CT&& columns ) const
   {
      using blaze::randomize;

      using ColumnsType = RemoveReference_t<CT>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ColumnsType );

      if( SO == false ) {
         for( size_t i=0UL; i<columns.rows(); ++i ) {
            for( size_t j=0UL; j<columns.columns(); ++j ) {
               randomize( columns(i,j) );
            }
         }
      }
      else {
         for( size_t j=0UL; j<columns.columns(); ++j ) {
            for( size_t i=0UL; i<columns.rows(); ++i ) {
               randomize( columns(i,j) );
            }
         }
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a dense column selection.
   //
   // \param columns The column selection to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename CT     // Type of the column selection
           , typename Arg >  // Min/max argument type
   inline void randomize( CT&& columns, const Arg& min, const Arg& max ) const
   {
      using blaze::randomize;

      using ColumnsType = RemoveReference_t<CT>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( ColumnsType );

      if( SO == false ) {
         for( size_t i=0UL; i<columns.rows(); ++i ) {
            for( size_t j=0UL; j<columns.columns(); ++j ) {
               randomize( columns(i,j), min, max );
            }
         }
      }
      else {
         for( size_t j=0UL; j<columns.columns(); ++j ) {
            for( size_t i=0UL; i<columns.rows(); ++i ) {
               randomize( columns(i,j), min, max );
            }
         }
      }
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  RAND SPECIALIZATION FOR SPARSE COLUMN SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for sparse column selections.
// \ingroup random
//
// This specialization of the Rand class randomizes sparse column selections.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CCAs >  // Compile time column arguments
class Rand< Columns<MT,SO,false,SF,CCAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a sparse column selection.
   //
   // \param columns The column selection to be randomized.
   // \return void
   */
   template< typename CT >  // Type of the column selection
   inline void randomize( CT&& columns ) const
   {
      using ColumnsType = RemoveReference_t<CT>;
      using ElementType = ElementType_t<ColumnsType>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ColumnsType );

      const size_t m( columns.rows()    );
      const size_t n( columns.columns() );

      if( m == 0UL || n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

      columns.reset();
      columns.reserve( nonzeros );

      while( columns.nonZeros() < nonzeros ) {
         columns( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse column selection.
   //
   // \param columns The column selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random column selection.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename CT >  // Type of the column selection
   inline void randomize( CT&& columns, size_t nonzeros ) const
   {
      using ColumnsType = RemoveReference_t<CT>;
      using ElementType = ElementType_t<ColumnsType>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ColumnsType );

      const size_t m( columns.rows()    );
      const size_t n( columns.columns() );

      if( nonzeros > m*n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( m == 0UL || n == 0UL ) return;

      columns.reset();
      columns.reserve( nonzeros );

      while( columns.nonZeros() < nonzeros ) {
         columns( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse column selection.
   //
   // \param columns The column selection to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename CT     // Type of the column selection
           , typename Arg >  // Min/max argument type
   inline void randomize( CT&& columns, const Arg& min, const Arg& max ) const
   {
      using ColumnsType = RemoveReference_t<CT>;
      using ElementType = ElementType_t<ColumnsType>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ColumnsType );

      const size_t m( columns.rows()    );
      const size_t n( columns.columns() );

      if( m == 0UL || n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

      columns.reset();
      columns.reserve( nonzeros );

      while( columns.nonZeros() < nonzeros ) {
         columns( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse column selection.
   //
   // \param columns The column selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random column selection.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename CT     // Type of the column selection
           , typename Arg >  // Min/max argument type
   inline void randomize( CT&& columns, size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      using ColumnsType = RemoveReference_t<CT>;
      using ElementType = ElementType_t<ColumnsType>;

      BLAZE_CONSTRAINT_MUST_BE_COLUMNS_TYPE( ColumnsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( ColumnsType );

      const size_t m( columns.rows()    );
      const size_t n( columns.columns() );

      if( nonzeros > m*n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( m == 0UL || n == 0UL ) return;

      columns.reset();
      columns.reserve( nonzeros );

      while( columns.nonZeros() < nonzeros ) {
         columns( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
