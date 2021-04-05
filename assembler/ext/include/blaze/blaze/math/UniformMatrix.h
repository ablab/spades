//=================================================================================================
/*!
//  \file blaze/math/UniformMatrix.h
//  \brief Header file for the complete UniformMatrix implementation
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

#ifndef _BLAZE_MATH_UNIFORMMATRIX_H_
#define _BLAZE_MATH_UNIFORMMATRIX_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/dense/UniformMatrix.h>
#include <blaze/math/DenseMatrix.h>
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
/*!\brief Specialization of the Rand class template for UniformMatrix.
// \ingroup random
//
// This specialization of the Rand class creates random instances of UniformMatrix.
*/
template< typename Type   // Data type of the matrix
        , bool SO         // Storage order
        , typename Tag >  // Type tag
class Rand< UniformMatrix<Type,SO,Tag> >
{
 public:
   //**********************************************************************************************
   /*!\brief Generation of a random UniformMatrix.
   //
   // \param m The number of rows of the random matrix.
   // \param n The number of columns of the random matrix.
   // \return The generated random matrix.
   */
   inline const UniformMatrix<Type,SO,Tag> generate( size_t m, size_t n ) const
   {
      UniformMatrix<Type,SO,Tag> matrix( m, n );
      randomize( matrix );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Generation of a random UniformMatrix.
   //
   // \param m The number of rows of the random matrix.
   // \param n The number of columns of the random matrix.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return The generated random matrix.
   */
   template< typename Arg >  // Min/max argument type
   inline const UniformMatrix<Type,SO,Tag>
      generate( size_t m, size_t n, const Arg& min, const Arg& max ) const
   {
      UniformMatrix<Type,SO,Tag> matrix( m, n );
      randomize( matrix, min, max );
      return matrix;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a UniformMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \return void
   */
   inline void randomize( UniformMatrix<Type,SO,Tag>& matrix ) const
   {
      matrix = rand<Type>();
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Randomization of a UniformMatrix.
   //
   // \param matrix The matrix to be randomized.
   // \param min The smallest possible value for a matrix element.
   // \param max The largest possible value for a matrix element.
   // \return void
   */
   template< typename Arg >  // Min/max argument type
   inline void randomize( UniformMatrix<Type,SO,Tag>& matrix,
                          const Arg& min, const Arg& max ) const
   {
      matrix = rand<Type>( min, max );
   }
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
