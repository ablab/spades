//=================================================================================================
/*!
//  \file blaze/math/views/columns/ColumnsData.h
//  \brief Header file for the implementation of the ColumnsData class template
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

#ifndef _BLAZE_MATH_VIEWS_COLUMNS_COLUMNSDATA_H_
#define _BLAZE_MATH_VIEWS_COLUMNS_COLUMNSDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Columns class.
// \ingroup columns
//
// The auxiliary ColumnsData class template represents an abstraction of the data members of the
// Columns class template. The necessary set of data members is selected depending on the number
// of compile time column arguments. The basic implementation of ColumnsData adapts the class
// template to the requirements of multiple compile time column arguments.
*/
template< typename... CCAs >  // Compile time column arguments
class ColumnsData
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
/*!\brief Specialization of the ColumnsData class template for index sequences.
// \ingroup columns
//
// This specialization of ColumnsData adapts the class template to the requirements of a non-zero
// number of compile time indices.
*/
template< size_t I        // First column index
        , size_t... Is >  // Remaining column indices
class ColumnsData< index_sequence<I,Is...> >
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
   template< typename... RCAs >
   explicit inline ColumnsData( RCAs... args ) noexcept;

   ColumnsData( const ColumnsData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnsData& operator=( const ColumnsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr decltype(auto) idces  () noexcept;
   static constexpr size_t         idx    ( size_t i ) noexcept;
   static constexpr size_t         columns() noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = std::array<size_t,N>;  //!< Type of the container for column indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static constexpr Indices indices_{ { I, Is... } };  //!< The indices of the columns in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if BLAZE_CPP14_MODE
// Definition and initialization of the static member variables
template< size_t I        // First column index
        , size_t... Is >  // Remaining column indices
constexpr typename ColumnsData< index_sequence<I,Is...> >::Indices
   ColumnsData< index_sequence<I,Is...> >::indices_;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnsData.
//
// \param args The optional column arguments.
*/
template< size_t I            // First column index
        , size_t... Is >      // Remaining column indices
template< typename... RCAs >  // Optional column arguments
inline ColumnsData< index_sequence<I,Is...> >::ColumnsData( RCAs... args ) noexcept
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified columns in the underlying matrix.
//
// \return A representation of the indices of the specified columns.
*/
template< size_t I        // First column index
        , size_t... Is >  // Remaining column indices
constexpr decltype(auto) ColumnsData< index_sequence<I,Is...> >::idces() noexcept
{
   return index_sequence<I,Is...>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified column in the underlying matrix.
//
// \param i Access index for the column.
// \return The index of the specified column.
*/
template< size_t I        // First column index
        , size_t... Is >  // Remaining column indices
constexpr size_t ColumnsData< index_sequence<I,Is...> >::idx( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns.
//
// \return The number of columns.
*/
template< size_t I        // First column index
        , size_t... Is >  // Remaining column indices
constexpr size_t ColumnsData< index_sequence<I,Is...> >::columns() noexcept
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
/*!\brief Specialization of the ColumnsData class template for index producing callables.
// \ingroup columns
//
// This specialization of ColumnsData adapts the class template to the requirements of index
// producing callables.
*/
template< typename P >  // Type of the index producer
class ColumnsData<P>
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
   template< typename... RCAs >
   inline ColumnsData( P p, size_t n, RCAs... args ) noexcept;

   ColumnsData( const ColumnsData& ) = default;
   ColumnsData( ColumnsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnsData& operator=( const ColumnsData& ) = delete;
   ColumnsData& operator=( ColumnsData&& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces  () const noexcept;
   inline size_t         idx    ( size_t i ) const noexcept;
   inline size_t         columns() const noexcept;
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
/*!\brief The constructor for ColumnsData.
//
// \param p A callable producing the indices.
// \param n The total number of indices.
// \param args The optional column arguments.
*/
template< typename P >        // Type of the index producer
template< typename... RCAs >  // Optional column arguments
inline ColumnsData<P>::ColumnsData( P p, size_t n, RCAs... args ) noexcept
   : p_( p )
   , n_( n )
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified columns in the underlying matrix.
//
// \return A representation of the indices of the specified columns.
*/
template< typename P >  // Type of the index producer
inline decltype(auto) ColumnsData<P>::idces() const noexcept
{
   return std::make_pair( p_, n_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified column in the underlying matrix.
//
// \param i Access index for the column.
// \return The index of the specified column.
*/
template< typename P >  // Type of the index producer
inline size_t ColumnsData<P>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   return p_(i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns.
//
// \return The number of columns.
*/
template< typename P >  // Type of the index producer
inline size_t ColumnsData<P>::columns() const noexcept
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME COLUMN ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ColumnsData class template for zero compile time column arguments.
// \ingroup columns
//
// This specialization of ColumnsData adapts the class template to the requirements of zero compile
// time column arguments.
*/
template<>
class ColumnsData<>
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
   template< typename T, typename... RCAs >
   inline ColumnsData( T* indices, size_t n, RCAs... args );

   ColumnsData( const ColumnsData& ) = default;
   ColumnsData( ColumnsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ColumnsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ColumnsData& operator=( const ColumnsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces  () const noexcept;
   inline size_t         idx    ( size_t i ) const noexcept;
   inline size_t         columns() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = SmallArray<size_t,8UL>;  //!< Type of the container for column indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Indices indices_;  //!< The indices of the columns in the matrix.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ColumnsData.
//
// \param indices Pointer to the first index of the selected columns.
// \param n The total number of indices.
// \param args The optional column arguments.
*/
template< typename T          // Type of the column indices
        , typename... RCAs >  // Optional column arguments
inline ColumnsData<>::ColumnsData( T* indices, size_t n, RCAs... args )
   : indices_( indices, indices+n )  // The indices of the columns in the matrix
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified columns in the underlying matrix.
//
// \return A representation of the indices of the specified columns.
*/
inline decltype(auto) ColumnsData<>::idces() const noexcept
{
   return const_cast<const Indices&>( indices_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified column in the underlying matrix.
//
// \param i Access index for the column.
// \return The index of the specified column.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active.
*/
inline size_t ColumnsData<>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < columns(), "Invalid column access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of columns.
//
// \return The number of columns.
*/
inline size_t ColumnsData<>::columns() const noexcept
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
/*!\brief Compares the indices of two ColumsData instances.
// \ingroup columns
//
// \param lhs The left-hand side instance for the comparison.
// \param rhs The right-hand side instance for the comparison.
// \return \a true if the indices of both instances are equal, \a false if not.
*/
template< typename... CRAs1, typename... CRAs2 >
constexpr bool
   compareIndices( const ColumnsData<CRAs1...>& lhs, const ColumnsData<CRAs2...>& rhs ) noexcept
{
   if( lhs.columns() != rhs.columns() )
      return false;

   for( size_t i=0UL; i<lhs.columns(); ++i ) {
      if( lhs.idx(i) != rhs.idx(i) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
