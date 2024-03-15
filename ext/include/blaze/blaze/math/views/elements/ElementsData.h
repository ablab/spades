//=================================================================================================
/*!
//  \file blaze/math/views/elements/ElementsData.h
//  \brief Header file for the implementation of the ElementsData class template
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTSDATA_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTSDATA_H_


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
/*!\brief Auxiliary class template for the data members of the Elements class.
// \ingroup elements
//
// The auxiliary ElementsData class template represents an abstraction of the data members of the
// Elements class template. The necessary set of data members is selected depending on the number
// of compile time element arguments. The basic implementation of ElementsData adapts the class
// template to the requirements of multiple compile time element arguments.
*/
template< typename... CEAs >  // Compile time element arguments
class ElementsData
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
/*!\brief Specialization of the ElementsData class template for index sequences.
// \ingroup elements
//
// This specialization of ElementsData adapts the class template to the requirements of a non-zero
// number of compile time indices.
*/
template< size_t I        // First element index
        , size_t... Is >  // Remaining element indices
class ElementsData< index_sequence<I,Is...> >
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
   template< typename... REAs >
   explicit inline ElementsData( REAs... args ) noexcept;

   ElementsData( const ElementsData& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ElementsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ElementsData& operator=( const ElementsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   static constexpr decltype(auto) idces() noexcept;
   static constexpr size_t         idx  ( size_t i ) noexcept;
   static constexpr size_t         size () noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = std::array<size_t,N>;  //!< Type of the container for element indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   static constexpr Indices indices_{ { I, Is... } };  //!< The indices of the elements in the vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
#if BLAZE_CPP14_MODE
// Definition and initialization of the static member variables
template< size_t I        // First element index
        , size_t... Is >  // Remaining element indices
constexpr typename ElementsData< index_sequence<I,Is...> >::Indices
   ElementsData< index_sequence<I,Is...> >::indices_;
#endif
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ElementsData.
//
// \param args The optional element arguments.
*/
template< size_t I            // First element index
        , size_t... Is >      // Remaining element indices
template< typename... REAs >  // Optional element arguments
inline ElementsData< index_sequence<I,Is...> >::ElementsData( REAs... args ) noexcept
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified elements in the underlying vector.
//
// \return A representation of the indices of the specified elements.
*/
template< size_t I        // First element index
        , size_t... Is >  // Remaining element indices
constexpr decltype(auto) ElementsData< index_sequence<I,Is...> >::idces() noexcept
{
   return index_sequence<I,Is...>();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified element in the underlying vector.
//
// \param i Access index for the element.
// \return The index of the specified element.
*/
template< size_t I        // First element index
        , size_t... Is >  // Remaining element indices
constexpr size_t ElementsData< index_sequence<I,Is...> >::idx( size_t i ) noexcept
{
   BLAZE_USER_ASSERT( i < size(), "Invalid element access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of elements.
//
// \return The number of elements.
*/
template< size_t I        // First element index
        , size_t... Is >  // Remaining element indices
constexpr size_t ElementsData< index_sequence<I,Is...> >::size() noexcept
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
/*!\brief Specialization of the ElementsData class template for index producing callables.
// \ingroup elements
//
// This specialization of ElementsData adapts the class template to the requirements of index
// producing callables.
*/
template< typename P >  // Type of the index producer
class ElementsData<P>
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
   template< typename... REAs >
   inline ElementsData( P p, size_t n, REAs... args ) noexcept;

   ElementsData( const ElementsData& ) = default;
   ElementsData( ElementsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ElementsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ElementsData& operator=( const ElementsData& ) = delete;
   ElementsData& operator=( ElementsData&& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces() const noexcept;
   inline size_t         idx  ( size_t i ) const noexcept;
   inline size_t         size () const noexcept;
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
/*!\brief The constructor for ElementsData.
//
// \param p A callable producing the indices.
// \param n The total number of indices.
// \param args The optional element arguments.
*/
template< typename P >        // Type of the index producer
template< typename... REAs >  // Optional element arguments
inline ElementsData<P>::ElementsData( P p, size_t n, REAs... args ) noexcept
   : p_( p )
   , n_( n )
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified elements in the underlying vector.
//
// \return A representation of the indices of the specified elements.
*/
template< typename P >  // Type of the index producer
inline decltype(auto) ElementsData<P>::idces() const noexcept
{
   return std::make_pair( p_, n_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified element in the underlying vector.
//
// \param i Access index for the element.
// \return The index of the specified element.
*/
template< typename P >  // Type of the index producer
inline size_t ElementsData<P>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < size(), "Invalid element access index" );
   return p_(i);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of elements.
//
// \return The number of elements.
*/
template< typename P >  // Type of the index producer
inline size_t ElementsData<P>::size() const noexcept
{
   return n_;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ZERO COMPILE TIME ELEMENT ARGUMENTS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of the ElementsData class template for zero compile time element arguments.
// \ingroup elements
//
// This specialization of ElementsData adapts the class template to the requirements of zero
// compile time element arguments.
*/
template<>
class ElementsData<>
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
   template< typename T, typename... REAs >
   inline ElementsData( T* indices, size_t n, REAs... args );

   ElementsData( const ElementsData& ) = default;
   ElementsData( ElementsData&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ElementsData() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ElementsData& operator=( const ElementsData& ) = delete;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline decltype(auto) idces() const noexcept;
   inline size_t         idx  ( size_t i ) const noexcept;
   inline size_t         size () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Type definitions****************************************************************************
   using Indices = SmallArray<size_t,8UL>;  //!< Type of the container for element indices.
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Indices indices_;  //!< The indices of the elements in the vector.
   //@}
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief The constructor for ElementsData.
//
// \param indices Pointer to the first index of the selected elements.
// \param n The total number of indices.
// \param args The optional element arguments.
*/
template< typename T          // Type of the element indices
        , typename... REAs >  // Optional element arguments
inline ElementsData<>::ElementsData( T* indices, size_t n, REAs... args )
   : indices_( indices, indices+n )  // The indices of the elements in the vector
{
   MAYBE_UNUSED( args... );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns a representation of the indices of the specified elements in the underlying vector.
//
// \return A representation of the indices of the specified elements.
*/
inline decltype(auto) ElementsData<>::idces() const noexcept
{
   return const_cast<const Indices&>( indices_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the index of the specified element in the underlying vector.
//
// \param i Access index for the element.
// \return The index of the specified element.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active.
*/
inline size_t ElementsData<>::idx( size_t i ) const noexcept
{
   BLAZE_USER_ASSERT( i < size(), "Invalid element access index" );
   return indices_[i];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of elements.
//
// \return The number of elements.
*/
inline size_t ElementsData<>::size() const noexcept
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
/*!\brief Compares the indices of two ElementsData instances.
// \ingroup elements
//
// \param lhs The left-hand side instance for the comparison.
// \param rhs The right-hand side instance for the comparison.
// \return \a true if the indices of both instances are equal, \a false if not.
*/
template< typename... CEAs1, typename... CEAs2 >
constexpr bool
   compareIndices( const ElementsData<CEAs1...>& lhs, const ElementsData<CEAs2...>& rhs ) noexcept
{
   if( lhs.size() != rhs.size() )
      return false;

   for( size_t i=0UL; i<lhs.size(); ++i ) {
      if( lhs.idx(i) != rhs.idx(i) )
         return false;
   }

   return true;
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
