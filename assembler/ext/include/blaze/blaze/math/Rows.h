//=================================================================================================
/*!
//  \file blaze/math/Rows.h
//  \brief Header file for the complete Rows implementation
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

#ifndef _BLAZE_MATH_ROWS_H_
#define _BLAZE_MATH_ROWS_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/Rows.h>
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
//  RAND SPECIALIZATION FOR DENSE ROW SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for dense row selections.
// \ingroup random
//
// This specialization of the Rand class randomizes dense row selections.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
class Rand< Rows<MT,SO,true,SF,CRAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a dense row selection.
   //
   // \param rows The row selection to be randomized.
   // \return void
   */
   template< typename RT >  // Type of the row selection
   inline void randomize( RT&& rows ) const
   {
      using blaze::randomize;

      using RowsType = RemoveReference_t<RT>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RowsType );

      if( SO == true ) {
         for( size_t i=0UL; i<rows.rows(); ++i ) {
            for( size_t j=0UL; j<rows.columns(); ++j ) {
               randomize( rows(i,j) );
            }
         }
      }
      else {
         for( size_t j=0UL; j<rows.columns(); ++j ) {
            for( size_t i=0UL; i<rows.rows(); ++i ) {
               randomize( rows(i,j) );
            }
         }
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a dense row selection.
   //
   // \param rows The row selection to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename RT     // Type of the row selection
           , typename Arg >  // Min/max argument type
   inline void randomize( RT&& rows, const Arg& min, const Arg& max ) const
   {
      using blaze::randomize;

      using RowsType = RemoveReference_t<RT>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE( RowsType );

      if( SO == true ) {
         for( size_t i=0UL; i<rows.rows(); ++i ) {
            for( size_t j=0UL; j<rows.columns(); ++j ) {
               randomize( rows(i,j), min, max );
            }
         }
      }
      else {
         for( size_t j=0UL; j<rows.columns(); ++j ) {
            for( size_t i=0UL; i<rows.rows(); ++i ) {
               randomize( rows(i,j), min, max );
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
//  RAND SPECIALIZATION FOR SPARSE ROW SELECTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the Rand class template for sparse row selections.
// \ingroup random
//
// This specialization of the Rand class randomizes sparse row selections.
*/
template< typename MT         // Type of the matrix
        , bool SO             // Storage order
        , bool SF             // Symmetry flag
        , typename... CRAs >  // Compile time row arguments
class Rand< Rows<MT,SO,false,SF,CRAs...> >
{
 public:
   //**********************************************************************************************
   /*!\brief Randomization of a sparse row selection.
   //
   // \param rows The row selection to be randomized.
   // \return void
   */
   template< typename RT >       // Type of the row selection
   inline void randomize( RT&& rows ) const
   {
      using RowsType    = RemoveReference_t<RT>;
      using ElementType = ElementType_t<RowsType>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RowsType );

      const size_t m( rows.rows()    );
      const size_t n( rows.columns() );

      if( m == 0UL || n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

      rows.reset();
      rows.reserve( nonzeros );

      while( rows.nonZeros() < nonzeros ) {
         rows( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse row selection.
   //
   // \param rows The row selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random row selection.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename RT >  // Type of the row selection
   inline void randomize( RT&& rows, size_t nonzeros ) const
   {
      using RowsType    = RemoveReference_t<RT>;
      using ElementType = ElementType_t<RowsType>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RowsType );

      const size_t m( rows.rows()    );
      const size_t n( rows.columns() );

      if( nonzeros > m*n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( m == 0UL || n == 0UL ) return;

      rows.reset();
      rows.reserve( nonzeros );

      while( rows.nonZeros() < nonzeros ) {
         rows( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>();
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse row selection.
   //
   // \param rows The row selection to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename RT     // Type of the row selection
           , typename Arg >  // Min/max argument type
   inline void randomize( RT&& rows, const Arg& min, const Arg& max ) const
   {
      using RowsType    = RemoveReference_t<RT>;
      using ElementType = ElementType_t<RowsType>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RowsType );

      const size_t m( rows.rows()    );
      const size_t n( rows.columns() );

      if( m == 0UL || n == 0UL ) return;

      const size_t nonzeros( rand<size_t>( 1UL, std::ceil( 0.5*m*n ) ) );

      rows.reset();
      rows.reserve( nonzeros );

      while( rows.nonZeros() < nonzeros ) {
         rows( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a sparse row selection.
   //
   // \param rows The row selection to be randomized.
   // \param nonzeros The number of non-zero elements of the random row selection.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   // \exception std::invalid_argument Invalid number of non-zero elements.
   */
   template< typename RT    // Type of the row selection
           , typename Arg > // Min/max argument type
   inline void randomize( RT&& rows, size_t nonzeros, const Arg& min, const Arg& max ) const
   {
      using RowsType    = RemoveReference_t<RT>;
      using ElementType = ElementType_t<RowsType>;

      BLAZE_CONSTRAINT_MUST_BE_ROWS_TYPE( RowsType );
      BLAZE_CONSTRAINT_MUST_BE_SPARSE_MATRIX_TYPE( RowsType );

      const size_t m( rows.rows()    );
      const size_t n( rows.columns() );

      if( nonzeros > m*n ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of non-zero elements" );
      }

      if( m == 0UL || n == 0UL ) return;

      rows.reset();
      rows.reserve( nonzeros );

      while( rows.nonZeros() < nonzeros ) {
         rows( rand<size_t>( 0UL, m-1UL ), rand<size_t>( 0UL, n-1UL ) ) = rand<ElementType>( min, max );
      }
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
