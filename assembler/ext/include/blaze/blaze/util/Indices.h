//=================================================================================================
/*!
//  \file blaze/util/Indices.h
//  \brief Header file for the Indices class
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

#ifndef _BLAZE_UTIL_INDICES_H_
#define _BLAZE_UTIL_INDICES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <blaze/util/Exception.h>
#include <blaze/util/Random.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Auxiliary class for the generation of random indices.
//
// This auxiliary class can be used to generate a set of random indices.
*/
template< typename T >  // Type of the indices
class Indices
{
 public:
   //**Type definitions****************************************************************************
   //! Iterator over the generated indices.
   using ConstIterator = typename std::vector<T>::const_iterator;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline Indices( T min, T max, T number );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t        size () const noexcept;
   inline ConstIterator begin() const noexcept;
   inline ConstIterator end  () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   std::vector<T> indices_;  //!< The generated indices.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The constructor for the Indices class.
//
// \param min The lower limit of the random indices.
// \param max The upper limit of the random indices.
// \param number The number of random indices to generate.
// \exception std::invalid_argument Invalid index range.
// \exception std::invalid_argument Invalid number of indices.
//
// This constructor initializes an Indices object by generating \a number random, unique indices
// in the range \a min to \a max. In case \a number is larger than the possible number of incides
// in the specified range, a \a std::invalid_argument exception is thrown.
*/
template< typename T >  // Type of the indices
inline Indices<T>::Indices( T min, T max, T number )
   : indices_()  // The generated indices
{
   if( max < min ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid index range" );
   }

   const T maxNumber( max + T(1) - min );

   if( number > maxNumber ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid number of indices" );
   }

   if( number == 0UL ) {
      return;
   }

   if( number <= T( maxNumber * 0.5 ) )
   {
      indices_.reserve( number );

      while( indices_.size() < number )
      {
         const T value = rand<T>(min,max);
         BLAZE_INTERNAL_ASSERT( min <= value && value <= max, "Invalid index detected" );
         const auto pos = std::lower_bound( indices_.begin(), indices_.end(), value );

         if( pos == indices_.end() || *pos != value ) {
            indices_.insert( pos, value );
         }
      }
   }
   else
   {
      indices_.resize( maxNumber );
      std::iota( indices_.begin(), indices_.end(), min );

      while( indices_.size() > number )
      {
         const T value = rand<T>(min,max);
         BLAZE_INTERNAL_ASSERT( min <= value && value <= max, "Invalid index detected" );
         const auto pos = std::lower_bound( indices_.begin(), indices_.end(), value );

         if( pos != indices_.end() && *pos == value ) {
            indices_.erase( pos );
         }
      }
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the total number of random indices.
//
// \return The total number of random indices.
*/
template< typename T >  // Type of the indices
inline size_t Indices<T>::size() const noexcept
{
   return indices_.size();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the beginning of the vector.
//
// \return Iterator to the beginning of the vector.
*/
template< typename T >  // Type of the indices
inline typename Indices<T>::ConstIterator Indices<T>::begin() const noexcept
{
   return indices_.begin();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last element of the vector.
//
// \return Iterator just past the last element of the vector.
*/
template< typename T >  // Type of the indices
inline typename Indices<T>::ConstIterator Indices<T>::end() const noexcept
{
   return indices_.end();
}
//*************************************************************************************************

} // namespace blaze

#endif
