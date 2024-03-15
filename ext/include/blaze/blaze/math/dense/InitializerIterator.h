//=================================================================================================
/*!
//  \file blaze/math/dense/InitializerIterator.h
//  \brief Header file for the InitializerIterator class template
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//
//  * The names of its contributors may not be used to endorse or promote products derived
//    from this software without specific prior written permission.
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

#ifndef _BLAZE_MATH_DENSE_INITIALIZERITERATOR_H_
#define _BLAZE_MATH_DENSE_INITIALIZERITERATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/InitializerList.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of an iterator for (extended) initializer lists.
// \ingroup math
//
// The InitializerIterator represents a generic random-access iterator for (extended) initializer
// lists. It can be used for initializer lists representing dense vectors and specific rows of
// dense matrices.
*/
template< typename Type >  // Type of the elements
class InitializerIterator
{
 public:
   //**Type definitions****************************************************************************
   using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
   using ValueType        = Type;                             //!< Type of the underlying elements.
   using PointerType      = const Type*;                      //!< Pointer return type.
   using ReferenceType    = const Type&;                      //!< Reference return type.
   using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

   // STL iterator requirements
   using iterator_category = IteratorCategory;  //!< The iterator category.
   using value_type        = ValueType;         //!< Type of the underlying elements.
   using pointer           = PointerType;       //!< Pointer return type.
   using reference         = ReferenceType;     //!< Reference return type.
   using difference_type   = DifferenceType;    //!< Difference between two iterators.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit inline InitializerIterator() noexcept;
   explicit inline InitializerIterator( size_t index, initializer_list<Type> list ) noexcept;

   InitializerIterator( const InitializerIterator& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~InitializerIterator() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline InitializerIterator& operator+=( ptrdiff_t inc ) noexcept;
   inline InitializerIterator& operator-=( ptrdiff_t dec ) noexcept;

   InitializerIterator& operator=( const InitializerIterator& ) = default;
   //@}
   //**********************************************************************************************

   //**Increment/decrement operators***************************************************************
   /*!\name Increment/decrement operators */
   //@{
   inline       InitializerIterator& operator++()      noexcept;
   inline const InitializerIterator  operator++( int ) noexcept;
   inline       InitializerIterator& operator--()      noexcept;
   inline const InitializerIterator  operator--( int ) noexcept;
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline ReferenceType operator[]( size_t index ) const noexcept;
   inline ReferenceType operator* () const noexcept;
   inline PointerType   operator->() const noexcept;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t                 index() const noexcept;
   inline initializer_list<Type> list () const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   size_t                 index_;  //!< Current index of the iterator within the initializer list.
   initializer_list<Type> list_;   //!< The initializer list to be traversed.

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type >  // Type of the elements
const Type InitializerIterator<Type>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for the InitializerIterator class.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>::InitializerIterator() noexcept
   : index_()  // Current index of the iterator within the initializer list
   , list_ ()  // The initializer list to be traversed
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the InitializerIterator class.
//
// \param index The initial index of the iterator within the initializer list.
// \param list The initializer list to be traversed.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>::InitializerIterator( size_t index, initializer_list<Type> list ) noexcept
   : index_( index )  // Current index of the iterator within the initializer list
   , list_ ( list  )  // The initializer list to be traversed
{}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Addition assignment operator.
//
// \param inc The increment of the iterator.
// \return The incremented iterator.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>& InitializerIterator<Type>::operator+=( ptrdiff_t inc ) noexcept
{
   index_ += inc;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator.
//
// \param dec The decrement of the iterator.
// \return The decremented iterator.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>& InitializerIterator<Type>::operator-=( ptrdiff_t dec ) noexcept
{
   index_ -= dec;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  INCREMENT/DECREMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Pre-increment operator.
//
// \return Reference to the incremented iterator.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>& InitializerIterator<Type>::operator++() noexcept
{
   ++index_;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Post-increment operator.
//
// \return The previous position of the iterator.
*/
template< typename Type >  // Type of the elements
inline const InitializerIterator<Type> InitializerIterator<Type>::operator++( int ) noexcept
{
   return InitializerIterator( index_++, list_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Pre-decrement operator.
//
// \return Reference to the decremented iterator.
*/
template< typename Type >  // Type of the elements
inline InitializerIterator<Type>& InitializerIterator<Type>::operator--() noexcept
{
   --index_;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Post-decrement operator.
//
// \return The previous position of the iterator.
*/
template< typename Type >  // Type of the elements
inline const InitializerIterator<Type> InitializerIterator<Type>::operator--( int ) noexcept
{
   return InitializerIterator( index_--, list_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Direct access to the underlying elements.
//
// \param index Access index.
// \return Reference to the accessed value.
*/
template< typename Type >  // Type of the elements
inline typename InitializerIterator<Type>::ReferenceType
   InitializerIterator<Type>::operator[]( size_t index ) const noexcept
{
   if( index < list_.size() )
      return list_.begin()[index];
   else
      return zero_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the element at the current iterator position.
//
// \return The resulting value.
*/
template< typename Type >  // Type of the elements
inline typename InitializerIterator<Type>::ReferenceType
   InitializerIterator<Type>::operator*() const noexcept
{
   if( index_ < list_.size() )
      return list_.begin()[index_];
   else
      return zero_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the element at the current iterator position.
//
// \return Pointer to the element at the current iterator position.
*/
template< typename Type >  // Type of the elements
inline typename InitializerIterator<Type>::PointerType
   InitializerIterator<Type>::operator->() const noexcept
{
   if( index_ < list_.size() )
      return list_.begin() + index_;
   else
      return &zero_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Low-level access to the underlying index member of the iterator.
//
// \return Index of the current iterator position.
*/
template< typename Type >  // Type of the elements
inline size_t InitializerIterator<Type>::index() const noexcept
{
   return index_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level access to the underlying list member of the iterator.
//
// \return The underlying initializer list.
*/
template< typename Type >  // Type of the elements
inline initializer_list<Type> InitializerIterator<Type>::list() const noexcept
{
   return list_;
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name InitializerIterator operators */
//@{
template< typename Type >
bool operator==( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
bool operator!=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
bool operator<( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
bool operator>( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
bool operator<=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
bool operator>=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;

template< typename Type >
const InitializerIterator<Type> operator+( const InitializerIterator<Type>& it, ptrdiff_t inc ) noexcept;

template< typename Type >
const InitializerIterator<Type> operator+( ptrdiff_t inc, const InitializerIterator<Type>& it ) noexcept;

template< typename Type >
const InitializerIterator<Type> operator-( const InitializerIterator<Type>& it, ptrdiff_t dec ) noexcept;

template< typename Type >
ptrdiff_t operator-( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the iterators refer to the same element, \a false if not.
*/
template< typename Type >  // Type of the elements
inline bool operator==( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() == rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the iterators don't refer to the same element, \a false if they do.
*/
template< typename Type >  // Type of the elements
inline bool operator!=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() != rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is smaller, \a false if not.
*/
template< typename Type >  // Type of the elements
inline bool operator<( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() < rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is greater, \a false if not.
*/
template< typename Type >  // Type of the elements
inline bool operator>( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() > rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
*/
template< typename Type >  // Type of the elements
inline bool operator<=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() <= rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two InitializerIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is greater or equal, \a false if not.
*/
template< typename Type >  // Type of the elements
inline bool operator>=( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() >= rhs.index();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between a InitializerIterator and an integral value.
//
// \param it The iterator to be incremented.
// \param inc The number of elements the iterator is incremented.
// \return The incremented iterator.
*/
template< typename Type >  // Type of the elements
inline const InitializerIterator<Type> operator+( const InitializerIterator<Type>& it, ptrdiff_t inc ) noexcept
{
   return InitializerIterator<Type>( it.index() + inc, it.list() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between an integral value and a InitializerIterator.
//
// \param inc The number of elements the iterator is incremented.
// \param it The iterator to be incremented.
// \return The incremented iterator.
*/
template< typename Type >  // Type of the elements
inline const InitializerIterator<Type> operator+( ptrdiff_t inc, const InitializerIterator<Type>& it ) noexcept
{
   return InitializerIterator<Type>( it.index() + inc, it.list() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between a InitializerIterator and an integral value.
//
// \param it The iterator to be decremented.
// \param dec The number of elements the iterator is decremented.
// \return The decremented iterator.
*/
template< typename Type >  // Type of the elements
inline const InitializerIterator<Type> operator-( const InitializerIterator<Type>& it, ptrdiff_t dec ) noexcept
{
   return InitializerIterator<Type>( it.index() - dec, it.list() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating the number of elements between two iterators.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return The number of elements between the two iterators.
*/
template< typename Type >  // Type of the elements
inline ptrdiff_t operator-( const InitializerIterator<Type>& lhs, const InitializerIterator<Type>& rhs ) noexcept
{
   return lhs.index() - rhs.index();
}
//*************************************************************************************************

} // namespace blaze

#endif
