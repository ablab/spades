//=================================================================================================
/*!
//  \file blaze/math/views/rows/RowsData.h
//  \brief Header file for the implementation of the RowsData class template
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

#ifndef _BLAZE_MATH_VIEWS_ROWS_ROWSDATA_H_
#define _BLAZE_MATH_VIEWS_ROWS_ROWSDATA_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <array>
#include <blaze/system/Standard.h>
#include <blaze/util/Assert.h>
#include <blaze/util/IntegerSequence.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/SmallArray.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Auxiliary class template for the data members of the Rows class.
// \ingroup rows
//
// The auxiliary RowsData class template represents an abstraction of the data members of the
// Rows class template. The necessary set of data members is selected depending on the number
// of compile time row arguments. The basic implementation of RowsData adapts the class template
// to the requirements of multiple compile time row arguments.
*/
template< typename... CRAs >  // Compile time row arguments
class RowsData
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR INDEX SEQUENCES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the RowsData class template for index sequences.
// \ingroup rows
//
// This specialization of RowsData adapts the class template to the requirements of a non-zero
// number of compile time indices.
*/
template< size_t I        // First row index
        , size_t... Is >  // Remaining row indices
class RowsData< index_sequence<I,Is...> >
{
 protected:
   //**Compile time flags**************************************************************************
   static constexpr size_t N = sizeof...( Is ) + 1UL;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Compile time flags**************************************************************************
   //! Compilation flag for compile time optimization.
   /*! The \a compileTimeArgs compilation flag indicates whether the view has been created by
       means of compile time arguments and whether these arguments can be queried at compile
       time. In that case, the \a compileTimeArgs compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool compileTimeArgs = true;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RRAs >
   explicit inline RowsData( RRAs... args ) noexcept;

   RowsData( const RowsData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~RowsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   RowsData& operator=( const RowsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr decltype(auto) idces() noexcept;
   static constexpr size_t         idx  ( size_t i ) noexcept;
   static constexpr size_t         rows () noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = std::array<size_t,N>;  //!< Type of the container for row indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static constexpr Indices indices_{ { I, Is... } };  //!< The indices of the rows in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if BLAZE_CPP14_MODE
// Definition and initialization of the static member variables
template< size_t I        // First row index
        , size_t... Is >  // Remaining row indices
constexpr typename RowsData< index_sequence<I,Is...> >::Indices
   RowsData< index_sequence<I,Is...> >::indices_;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for RowsData.
//
// \param args The optional row arguments.
*/
template< size_t I            // First row index
        , size_t... Is >      // Remaining row indices
template< typename... RRAs >  // Optional row arguments
inline RowsData< index_sequence<I,Is...> >::RowsData( RRAs... args ) noexcept
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified rows in the underlying matrix.
//
// \return A representation of the indices of the specified rows.
*/
template< size_t I        // First row index
        , size_t... Is >  // Remaining row indices
constexpr decltype(auto) RowsData< index_sequence<I,Is...> >::idces() noexcept
{
   return index_sequence<I,Is...>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified row in the underlying matrix.
//
// \param i Access index for the row.
// \return The index of the specified row.
*/
template< size_t I        // First row index
        , size_t... Is >  // Remaining row indices
constexpr size_t RowsData< index_sequence<I,Is...> >::idx( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows.
//
// \return The number of rows.
*/
template< size_t I        // First row index
        , size_t... Is >  // Remaining row indices
constexpr size_t RowsData< index_sequence<I,Is...> >::rows() noexcept
{
   return N;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR INDEX PRODUCING CALLABLES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the RowsData class template for index producing callables.
// \ingroup rows
//
// This specialization of RowsData adapts the class template to the requirements of index
// producing callables.
*/
template< typename P >  // Type of the index producer
class RowsData<P>
{
 protected:
   //**Compile time flags**************************************************************************
   static constexpr size_t N = 0UL;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Compile time flags**************************************************************************
   //! Compilation flag for compile time optimization.
   /*! The \a compileTimeArgs compilation flag indicates whether the view has been created by
       means of compile time arguments and whether these arguments can be queried at compile
       time. In that case, the \a compileTimeArgs compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool compileTimeArgs = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RRAs >
   inline RowsData( P p, size_t n, RRAs... args ) noexcept;

   RowsData( const RowsData& ) = default;
   RowsData( RowsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~RowsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   RowsData& operator=( const RowsData& ) = delete;
   RowsData& operator=( RowsData&& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces() const noexcept;
   inline size_t         idx  ( size_t i ) const noexcept;
   inline size_t         rows () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   P      p_;  //!< The callable producing the indices.
   size_t n_;  //!< The total number of indices.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for RowsData.
//
// \param p A callable producing the indices.
// \param n The total number of indices.
// \param args The optional row arguments.
*/
template< typename P >        // Type of the index producer
template< typename... RRAs >  // Optional row arguments
inline RowsData<P>::RowsData( P p, size_t n, RRAs... args ) noexcept
   : p_( p )
   , n_( n )
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified rows in the underlying matrix.
//
// \return A representation of the indices of the specified rows.
*/
template< typename P >  // Type of the index producer
inline decltype(auto) RowsData<P>::idces() const noexcept
{
   return std::make_pair( p_, n_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified row in the underlying matrix.
//
// \param i Access index for the row.
// \return The index of the specified row.
*/
template< typename P >  // Type of the index producer
inline size_t RowsData<P>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return p_(i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows.
//
// \return The number of rows.
*/
template< typename P >  // Type of the index producer
inline size_t RowsData<P>::rows() const noexcept
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME ROW ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the RowsData class template for zero compile time row arguments.
// \ingroup rows
//
// This specialization of RowsData adapts the class template to the requirements of zero compile
// time row arguments.
*/
template<>
class RowsData<>
{
 protected:
   //**Compile time flags**************************************************************************
   static constexpr size_t N = 0UL;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Compile time flags**************************************************************************
   //! Compilation flag for compile time optimization.
   /*! The \a compileTimeArgs compilation flag indicates whether the view has been created by
       means of compile time arguments and whether these arguments can be queried at compile
       time. In that case, the \a compileTimeArgs compilation flag is set to \a true, otherwise
       it is set to \a false. */
   static constexpr bool compileTimeArgs = false;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename T, typename... RRAs >
   inline RowsData( T* indices, size_t n, RRAs... args );

   RowsData( const RowsData& ) = default;
   RowsData( RowsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~RowsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   RowsData& operator=( const RowsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces() const noexcept;
   inline size_t         idx  ( size_t i ) const noexcept;
   inline size_t         rows () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = SmallArray<size_t,8UL>;  //!< Type of the container for row indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Indices indices_;  //!< The indices of the rows in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for RowsData.
//
// \param indices Pointer to the first index of the selected rows.
// \param n The total number of indices.
// \param args The optional row arguments.
*/
template< typename T          // Type of the row indices
        , typename... RRAs >  // Optional row arguments
inline RowsData<>::RowsData( T* indices, size_t n, RRAs... args )
   : indices_( indices, indices+n )  // The indices of the rows in the matrix
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified rows in the underlying matrix.
//
// \return A representation of the indices of the specified rows.
*/
inline decltype(auto) RowsData<>::idces() const noexcept
{
   return const_cast<const Indices&>( indices_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified row in the underlying matrix.
//
// \param i Access index for the row.
// \return The index of the specified row.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active.
*/
inline size_t RowsData<>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of rows.
//
// \return The number of rows.
*/
inline size_t RowsData<>::rows() const noexcept
{
   return indices_.size();
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Compares the indices of two RowsData instances.
// \ingroup rows
//
// \param lhs The left-hand side instance for the comparison.
// \param rhs The right-hand side instance for the comparison.
// \return \a true if the indices of both instances are equal, \a false if not.
*/
template< typename... CRAs1, typename... CRAs2 >
constexpr bool
   compareIndices( const RowsData<CRAs1...>& lhs, const RowsData<CRAs2...>& rhs ) noexcept
{
   if( lhs.rows() != rhs.rows() )
      return false;

   for( size_t i=0UL; i<lhs.rows(); ++i ) {
      if( lhs.idx(i) != rhs.idx(i) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
