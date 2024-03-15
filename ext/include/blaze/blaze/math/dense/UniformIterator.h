//=================================================================================================
/*!
//  \file blaze/math/dense/UniformIterator.h
//  \brief Header file for the UniformIterator class template
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

#ifndef _BLAZE_MATH_DENSE_UNIFORMITERATOR_H_
#define _BLAZE_MATH_DENSE_UNIFORMITERATOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/SIMD.h>
#include <blaze/util/Assert.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of a generic iterator for uniform vectors and matrices.
// \ingroup math
//
// The UniformIterator represents a generic random-access iterator that can be used for uniform
// vectors and specific rows/columns of uniform matrices.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
class UniformIterator
{
 public:
   //**Type definitions****************************************************************************
   using IteratorCategory = std::random_access_iterator_tag;  //!< The iterator category.
   using ValueType        = Type;                             //!< Type of the underlying elements.
   using PointerType      = Type*;                            //!< Pointer return type.
   using ReferenceType    = Type&;                            //!< Reference return type.
   using DifferenceType   = ptrdiff_t;                        //!< Difference between two iterators.

   // STL iterator requirements
   using iterator_category = IteratorCategory;  //!< The iterator category.
   using value_type        = ValueType;         //!< Type of the underlying elements.
   using pointer           = PointerType;       //!< Pointer return type.
   using reference         = ReferenceType;     //!< Reference return type.
   using difference_type   = DifferenceType;    //!< Difference between two iterators.

   //! SIMD type of the elements.
   using SIMDType = SIMDTrait_t<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   explicit constexpr UniformIterator() noexcept;
   explicit constexpr UniformIterator( Type* ptr, size_t index ) noexcept;

   template< typename Other, AlignmentFlag AF2 >
   constexpr UniformIterator( const UniformIterator<Other,AF2>& it ) noexcept;

   UniformIterator( const UniformIterator& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~UniformIterator() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   constexpr UniformIterator& operator+=( ptrdiff_t inc ) noexcept;
   constexpr UniformIterator& operator-=( ptrdiff_t inc ) noexcept;

   UniformIterator& operator=( const UniformIterator& ) = default;
   //@}
   //**********************************************************************************************

   //**Increment/decrement operators***************************************************************
   /*!\name Increment/decrement operators */
   //@{
   constexpr UniformIterator&      operator++()      noexcept;
   constexpr const UniformIterator operator++( int ) noexcept;
   constexpr UniformIterator&      operator--()      noexcept;
   constexpr const UniformIterator operator--( int ) noexcept;
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   constexpr ReferenceType operator[]( size_t index ) const noexcept;
   constexpr ReferenceType operator* () const noexcept;
   constexpr PointerType   operator->() const noexcept;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   constexpr PointerType ptr() const noexcept;
   constexpr size_t      idx() const noexcept;
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   inline const SIMDType load () const noexcept;
   inline const SIMDType loada() const noexcept;
   inline const SIMDType loadu() const noexcept;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   PointerType ptr_;    //!< Pointer to the element.
   size_t      index_;  //!< Index of the current element.
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
/*!\brief Default constructor for the UniformIterator class.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>::UniformIterator() noexcept
   : ptr_  ( nullptr )  // Pointer to the element
   , index_( 0UL )      // Index of the current element
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for the UniformIterator class.
//
// \param ptr Pointer to the element.
// \param index The index of the initial element.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>::UniformIterator( Type* ptr, size_t index ) noexcept
   : ptr_  ( ptr )    // Pointer to the element
   , index_( index )  // Index of the current element
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different UniformIterator instances.
//
// \param it The foreign UniformIterator instance to be copied.
*/
template< typename Type        // Type of the elements
        , AlignmentFlag AF >   // Alignment flag
template< typename Other       // Type of the foreign elements
        , AlignmentFlag AF2 >  // Alignment flag of the foreign iterator
constexpr UniformIterator<Type,AF>::UniformIterator( const UniformIterator<Other,AF2>& it ) noexcept
   : ptr_  ( it.ptr() )  // Pointer to the element
   , index_( it.idx() )  // Index of the current element
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
// \return Reference to the incremented iterator.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>&
   UniformIterator<Type,AF>::operator+=( ptrdiff_t inc ) noexcept
{
   index_ += inc;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator.
//
// \param dec The decrement of the iterator.
// \return Reference to the decremented iterator.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>&
   UniformIterator<Type,AF>::operator-=( ptrdiff_t dec ) noexcept
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
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>&
   UniformIterator<Type,AF>::operator++() noexcept
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
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr const UniformIterator<Type,AF>
   UniformIterator<Type,AF>::operator++( int ) noexcept
{
   return UniformIterator( ptr_, index_++ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Pre-decrement operator.
//
// \return Reference to the decremented iterator.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr UniformIterator<Type,AF>&
   UniformIterator<Type,AF>::operator--() noexcept
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
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr const UniformIterator<Type,AF>
   UniformIterator<Type,AF>::operator--( int ) noexcept
{
   return UniformIterator( ptr_, index_-- );
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
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr typename UniformIterator<Type,AF>::ReferenceType
   UniformIterator<Type,AF>::operator[]( size_t index ) const noexcept
{
   MAYBE_UNUSED( index );

   return *ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the element at the current iterator position.
//
// \return Reference to the current element.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr typename UniformIterator<Type,AF>::ReferenceType
   UniformIterator<Type,AF>::operator*() const noexcept
{
   return *ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the element at the current iterator position.
//
// \return Pointer to the element at the current iterator position.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr typename UniformIterator<Type,AF>::PointerType
   UniformIterator<Type,AF>::operator->() const noexcept
{
   return ptr_;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Low-level access to the current memory location of the iterator.
//
// \return Pointer to the current memory location.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr typename UniformIterator<Type,AF>::PointerType
   UniformIterator<Type,AF>::ptr() const noexcept
{
   return ptr_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Low-level access to the underlying index of the iterator.
//
// \return Index to the current element.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
constexpr size_t
   UniformIterator<Type,AF>::idx() const noexcept
{
   return index_;
}
//*************************************************************************************************





//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Load of the SIMD element at the current iterator position.
//
// \return The loaded SIMD element.
//
// This function performs a load of the SIMD element of the current element. This function must
// \b NOT be called explicitly! It is used internally for the performance optimized evaluation
// of expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
inline const typename UniformIterator<Type,AF>::SIMDType
   UniformIterator<Type,AF>::load() const noexcept
{
   if( AF )
      return loada();
   else
      return loadu();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Aligned load of the SIMD element at the current iterator position.
//
// \return The loaded SIMD element.
//
// This function performs an aligned load of the SIMD element of the current element. This
// function must \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result in erroneous
// results and/or in compilation errors.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
inline const typename UniformIterator<Type,AF>::SIMDType
   UniformIterator<Type,AF>::loada() const noexcept
{
   return blaze::set( *ptr_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Unaligned load of the SIMD element at the current iterator position.
//
// \return The loaded SIMD element.
//
// This function performs an unaligned load of the SIMD element of the current element. This
// function must \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result in erroneous
// results and/or in compilation errors.
*/
template< typename Type       // Type of the elements
        , AlignmentFlag AF >  // Alignment flag
inline const typename UniformIterator<Type,AF>::SIMDType
   UniformIterator<Type,AF>::loadu() const noexcept
{
   return blaze::set( *ptr_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name UniformIterator operators */
//@{
template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator==( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator!=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator<( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator>( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator<=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename T1, AlignmentFlag AF1, typename T2, AlignmentFlag AF2 >
constexpr bool
   operator>=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept;

template< typename Type, AlignmentFlag AF >
constexpr const UniformIterator<Type,AF>
   operator+( const UniformIterator<Type,AF>& it, ptrdiff_t inc ) noexcept;

template< typename Type, AlignmentFlag AF >
constexpr const UniformIterator<Type,AF>
   operator+( ptrdiff_t inc, const UniformIterator<Type,AF>& it ) noexcept;

template< typename Type, AlignmentFlag AF >
constexpr const UniformIterator<Type,AF>
   operator-( const UniformIterator<Type,AF>& it, ptrdiff_t inc ) noexcept;

template< typename Type, AlignmentFlag AF >
constexpr ptrdiff_t
   operator-( const UniformIterator<Type,AF>& lhs, const UniformIterator<Type,AF>& rhs ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the iterators refer to the same element, \a false if not.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator==( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() == rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the iterators don't refer to the same element, \a false if they do.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator!=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() != rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is smaller, \a false if not.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator<( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() == rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is greater, \a false if not.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator>( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() > rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is less or equal, \a false if not.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator<=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() <= rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return \a true if the left-hand side iterator is greater or equal, \a false if not.
*/
template< typename T1          // Element type of the left-hand side iterator
        , AlignmentFlag AF1    // Alignment flag of the left-hand side iterator
        , typename T2          // Element type of the right-hand side iterator
        , AlignmentFlag AF2 >  // Alignment flag of the right-hand side iterator
constexpr bool
   operator>=( const UniformIterator<T1,AF1>& lhs, const UniformIterator<T2,AF2>& rhs ) noexcept
{
   return lhs.idx() >= rhs.idx();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between a UniformIterator and an integral value.
//
// \param it The iterator to be incremented.
// \param inc The number of elements the iterator is incremented.
// \return The incremented iterator.
*/
template< typename Type       // Element type of the iterator
        , AlignmentFlag AF >  // Alignment flag of the iterator
constexpr const UniformIterator<Type,AF>
   operator+( const UniformIterator<Type,AF>& it, ptrdiff_t inc ) noexcept
{
   return UniformIterator<Type,AF>( it.ptr(), it.idx() + inc );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between an integral value and a UniformIterator.
//
// \param inc The number of elements the iterator is incremented.
// \param it The iterator to be incremented.
// \return The incremented iterator.
*/
template< typename Type       // Element type of the iterator
        , AlignmentFlag AF >  // Alignment flag of the iterator
constexpr const UniformIterator<Type,AF>
   operator+( ptrdiff_t inc, const UniformIterator<Type,AF>& it ) noexcept
{
   return UniformIterator<Type,AF>( it.ptr(), it.idx() + inc );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between a UniformIterator and an integral value.
//
// \param it The iterator to be decremented.
// \param dec The number of elements the iterator is decremented.
// \return The decremented iterator.
*/
template< typename Type       // Element type of the iterator
        , AlignmentFlag AF >  // Alignment flag of the iterator
constexpr const UniformIterator<Type,AF>
   operator-( const UniformIterator<Type,AF>& it, ptrdiff_t dec ) noexcept
{
   return UniformIterator<Type,AF>( it.ptr(), it.idx() - dec );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating the number of elements between two UniformIterator objects.
//
// \param lhs The left-hand side iterator.
// \param rhs The right-hand side iterator.
// \return The number of elements between the two iterators.
*/
template< typename Type       // Element type of the iterator
        , AlignmentFlag AF >  // Alignment flag of the iterator
constexpr ptrdiff_t
   operator-( const UniformIterator<Type,AF>& lhs, const UniformIterator<Type,AF>& rhs ) noexcept
{
   return lhs.idx() - rhs.idx();
}
//*************************************************************************************************

} // namespace blaze

#endif
