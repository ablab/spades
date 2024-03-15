//=================================================================================================
/*!
//  \file blaze/math/views/elements/Dense.h
//  \brief Elements specialization for dense vectors
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_DENSE_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Elements.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Subvector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/elements/BaseTemplate.h>
#include <blaze/math/views/elements/ElementsData.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/Assert.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DENSE VECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for dense vectors.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the requirements of dense vectors.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements<VT,TF,true,CEAs...>
   : public View< DenseVector< Elements<VT,TF,true,CEAs...>, TF > >
   , private ElementsData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = ElementsData<CEAs...>;                //!< The type of the ElementsData base class.
   using Operand  = If_t< IsExpression_v<VT>, VT, VT& >;  //!< Composite data type of the vector expression.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //!< Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Elements instance.
   using This = Elements<VT,TF,true,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = VT;                           //!< The type viewed by this Elements instance.
   using ResultType    = ElementsTrait_t<VT,N>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;            //!< Type of the elements.
   using ReturnType    = ReturnType_t<VT>;             //!< Return type for expression template evaluations
   using CompositeType = const Elements&;              //!< Data type for composite expression templates.

   //! Reference to a constant element value.
   using ConstReference = ConstReference_t<VT>;

   //! Reference to a non-constant element value.
   using Reference = If_t< IsConst_v<VT>, ConstReference, Reference_t<VT> >;

   //! Pointer to a constant element value.
   using ConstPointer = ConstPointer_t<VT>;

   //! Pointer to a non-constant element value.
   using Pointer = If_t< IsConst_v<VT> || !HasMutableDataAccess_v<VT>, ConstPointer, Pointer_t<VT> >;
   //**********************************************************************************************

   //**ElementsIterator class definition***********************************************************
   /*!\brief Iterator over the element selection of the dense vector.
   */
   template< typename ElementsType    // Type of the element selection
           , typename IteratorType >  // Type of the dense vector iterator
   class ElementsIterator
   {
    public:
      //**Type definitions*************************************************************************
      //! The iterator category.
      using IteratorCategory = typename std::iterator_traits<IteratorType>::iterator_category;

      //! Type of the underlying elements.
      using ValueType = typename std::iterator_traits<IteratorType>::value_type;

      //! Pointer return type.
      using PointerType = typename std::iterator_traits<IteratorType>::pointer;

      //! Reference return type.
      using ReferenceType = typename std::iterator_traits<IteratorType>::reference;

      //! Difference between two iterators.
      using DifferenceType = typename std::iterator_traits<IteratorType>::difference_type;

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Default constructor of the ElementsIterator class.
      */
      inline ElementsIterator()
         : elements_( nullptr )  // Pointer to the element selection
         , index_   ( 0UL )      // Index of the current element of the element selection
         , pos_     ()           // Iterator to the current element
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ElementsIterator class.
      //
      // \param elements Pointer to the element selection.
      // \param index Index of the initial element.
      */
      inline ElementsIterator( ElementsType* elements, size_t index )
         : elements_( elements )  // Pointer to the element selection
         , index_   ( index )     // Index of the current element of the element selection
         , pos_     ()            // Iterator to the current element
      {
         if( index_ < elements_->size() )
            pos_ = elements_->operand().begin() + elements_->idx( index_ );
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ElementsIterator instances.
      //
      // \param it The elements iterator to be copied
      */
      template< typename ElementsType2, typename IteratorType2 >
      inline ElementsIterator( const ElementsIterator<ElementsType2,IteratorType2>& it )
         : elements_( it.elements_ )  // Pointer to the element selection
         , index_   ( it.index_ )     // Index of the current element of the element selection
         , pos_     ( it.pos_ )       // Iterator to the current element
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline ElementsIterator& operator+=( size_t inc ) {
         using blaze::reset;
         index_ += inc;
         if( index_ < elements_->size() )
            pos_ = elements_->operand().begin() + elements_->idx( index_ );
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline ElementsIterator& operator-=( size_t dec ) {
         using blaze::reset;
         index_ -= dec;
         if( index_ < elements_->size() )
            pos_ = elements_->operand().begin() + elements_->idx( index_ );
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ElementsIterator& operator++() {
         using blaze::reset;
         ++index_;
         if( index_ < elements_->size() )
            pos_ = elements_->operand().begin() + elements_->idx( index_ );
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ElementsIterator operator++( int ) {
         const ElementsIterator tmp( *this );
         ++(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline ElementsIterator& operator--() {
         using blaze::reset;
         --index_;
         if( index_ < elements_->size() )
            pos_ = elements_->operand().begin() + elements_->idx( index_ );
         else reset( pos_ );
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const ElementsIterator operator--( int ) {
         const ElementsIterator tmp( *this );
         --(*this);
         return tmp;
      }
      //*******************************************************************************************

      //**Subscript operator***********************************************************************
      /*!\brief Direct access to the elements.
      //
      // \param index Access index.
      // \return Reference to the accessed value.
      */
      inline ReferenceType operator[]( size_t index ) const {
         const IteratorType pos( elements_->operand().begin() + elements_->idx( index_ + index ) );
         return *pos;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReferenceType operator*() const {
         return *pos_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return Pointer to the element at the current iterator position.
      */
      inline PointerType operator->() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const ElementsIterator& rhs ) const {
         return index_ == rhs.index_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const ElementsIterator& rhs ) const {
         return index_ != rhs.index_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const ElementsIterator& rhs ) const {
         return index_ < rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const ElementsIterator& rhs ) const {
         return index_ > rhs.index_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const ElementsIterator& rhs ) const {
         return index_ <= rhs.index_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two ElementsIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const ElementsIterator& rhs ) const {
         return index_ >= rhs.index_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ElementsIterator& rhs ) const {
         return index_ - rhs.index_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a ElementsIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const ElementsIterator operator+( const ElementsIterator& it, size_t inc ) {
         return ElementsIterator( it.elements_, it.index_ + inc );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a ElementsIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const ElementsIterator operator+( size_t inc, const ElementsIterator& it ) {
         return ElementsIterator( it.elements_, it.index_ + inc );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a ElementsIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const ElementsIterator operator-( const ElementsIterator& it, size_t dec ) {
         return ElementsIterator( it.elements_, it.index_ - dec );
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      ElementsType* elements_;  //!< Pointer to the element selection.
      size_t        index_;     //!< Index of the current element of the element selection.
      IteratorType  pos_;       //!< Iterator to the current element.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename ElementsType2, typename IteratorType2 > friend class ElementsIterator;
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = ElementsIterator< const This, ConstIterator_t<VT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<VT>, ConstIterator, ElementsIterator< This, Iterator_t<VT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... REAs >
   explicit inline Elements( VT& vector, REAs... args );

   Elements( const Elements& ) = default;
   Elements( Elements&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Elements() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index );
   inline ConstReference operator[]( size_t index ) const;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Iterator       begin ();
   inline ConstIterator  begin () const;
   inline ConstIterator  cbegin() const;
   inline Iterator       end   ();
   inline ConstIterator  end   () const;
   inline ConstIterator  cend  () const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                            inline Elements& operator= ( const ElementType& rhs );
                            inline Elements& operator= ( initializer_list<ElementType> list );
                            inline Elements& operator= ( const Elements& rhs );
   template< typename VT2 > inline Elements& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Elements& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Elements& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Elements& operator*=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Elements& operator/=( const DenseVector<VT2,TF>& rhs );
   template< typename VT2 > inline Elements& operator%=( const Vector<VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::idces;
   using DataType::idx;
   using DataType::size;

   inline VT&       operand() noexcept;
   inline const VT& operand() const noexcept;

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Elements& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   template< typename VT2 > inline void assign    ( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 > inline void assign    ( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 > inline void addAssign ( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 > inline void addAssign ( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 > inline void subAssign ( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 > inline void subAssign ( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 > inline void multAssign( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 > inline void multAssign( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 > inline void divAssign ( const DenseVector <VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand vector_;  //!< The vector containing the elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBVECTOR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_ELEMENTS_TYPE   ( VT );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( VT, TF );
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Constructor for elements on dense vectors.
//
// \param vector The vector containing the elements.
// \param args The runtime element arguments.
// \exception std::invalid_argument Invalid element access index.
//
// By default, the provided element indices are checked at runtime. In case any element is not
// properly specified (i.e. if any specified index is greater than the number of elements of the
// given vector) a \a std::invalid_argument exception is thrown. The checks can be skipped by
// providing the optional \a blaze::unchecked argument.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename... REAs >  // Optional arguments
inline Elements<VT,TF,true,CEAs...>::Elements( VT& vector, REAs... args )
   : DataType( args... )  // Base class initialization
   , vector_ ( vector  )  // The vector containing the elements
{
   if( isChecked( args... ) ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( vector_.size() <= idx(i) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the elements.
//
// \param index Access index. The index must be smaller than the number of elements.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::Reference
   Elements<VT,TF,true,CEAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid element access index" );
   return vector_[idx(index)];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the elements.
//
// \param index Access index. The index must be smaller than the number of elements.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstReference
   Elements<VT,TF,true,CEAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid element access index" );
   return const_cast<const VT&>( vector_ )[idx(index)];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the elements.
//
// \param index Access index. The index must be smaller than the number of elements.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid element access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::Reference
   Elements<VT,TF,true,CEAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid element access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the elements.
//
// \param index Access index. The index must be smaller than the number of elements.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid element access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstReference
   Elements<VT,TF,true,CEAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid element access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the element selection.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the elements.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::Pointer
   Elements<VT,TF,true,CEAs...>::data() noexcept
{
   return vector_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the element selection.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the elements.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstPointer
   Elements<VT,TF,true,CEAs...>::data() const noexcept
{
   return vector_.data() + idx(0UL);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the element selection.
//
// \return Iterator to the first element of the element selection.
//
// This function returns an iterator to the first element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::Iterator
   Elements<VT,TF,true,CEAs...>::begin()
{
   return Iterator( this, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the element selection.
//
// \return Iterator to the first element of the element selection.
//
// This function returns an iterator to the first element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstIterator
   Elements<VT,TF,true,CEAs...>::begin() const
{
   return ConstIterator( this, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the element selection.
//
// \return Iterator to the first element of the element selection.
//
// This function returns an iterator to the first element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstIterator
   Elements<VT,TF,true,CEAs...>::cbegin() const
{
   return ConstIterator( this, 0UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the element selection.
//
// \return Iterator just past the last element of the element selection.
//
// This function returns an iterator just past the last element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::Iterator
   Elements<VT,TF,true,CEAs...>::end()
{
   return Iterator( this, size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the element selection.
//
// \return Iterator just past the last element of the element selection.
//
// This function returns an iterator just past the last element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstIterator
   Elements<VT,TF,true,CEAs...>::end() const
{
   return ConstIterator( this, size() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the element selection.
//
// \return Iterator just past the last element of the element selection.
//
// This function returns an iterator just past the last element of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,true,CEAs...>::ConstIterator
   Elements<VT,TF,true,CEAs...>::cend() const
{
   return ConstIterator( this, size() );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Homogenous assignment to all elements.
//
// \param rhs Scalar value to be assigned to all elements.
// \return Reference to the assigned element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator=( const ElementType& rhs )
{
   decltype(auto) left( derestrict( vector_ ) );

   for( size_t i=0UL; i<size(); ++i ) {
      const size_t index( idx(i) );
      if( !IsRestricted_v<VT> || trySet( vector_, index, rhs ) )
         left[index] = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all elements.
//
// \param list The initializer list.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Invalid assignment to element selection.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// This assignment operator offers the option to directly assign to all elements of the element
// selection by means of an initializer list. The elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the element selection, a \a std::invalid_argument
// exception is thrown. Also, if the underlying vector \a VT is restricted and the assignment
// would violate an invariant of the vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to elements" );
   }

   const InitializerVector<ElementType,TF> tmp( list, size() );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !trySet( vector_, idx(i), tmp[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( vector_ ) );
   for( size_t i=0UL; i<size(); ++i ) {
      left[idx(i)] = tmp[i];
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Elements.
//
// \param rhs Dense element selection to be copied.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two element selections don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator=( const Elements& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && compareIndices( *this, rhs ) ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !trySet( vector_, idx(i), rhs[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different vectors.
//
// \param rhs Vector to be assigned.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !trySet( vector_, idx(i), right[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      smpAssign( left, tmp );
   }
   else {
      if( IsSparseVector_v<VT2> )
         reset();
      smpAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the dense element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator+=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !tryAdd( vector_, idx(i), right[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      smpAddAssign( left, tmp );
   }
   else {
      smpAddAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the dense element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator-=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !trySub( vector_, idx(i), right[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      smpSubAssign( left, tmp );
   }
   else {
      smpSubAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Multiplication assignment operator for the multiplication of a vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be multiplied with the dense element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator*=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !tryMult( vector_, idx(i), right[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      smpMultAssign( left, tmp );
   }
   else {
      smpMultAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator/=( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !tryDiv( vector_, idx(i), right[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      smpDivAssign( left, tmp );
   }
   else {
      smpDivAssign( left, right );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::operator%=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   using CrossType = CrossTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size() != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   const CrossType tmp( *this % (*rhs) );

   if( IsRestricted_v<VT> ) {
      for( size_t i=0UL; i<size(); ++i ) {
         if( !trySet( vector_, idx(i), tmp[i] ) ) {
            BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
         }
      }
   }

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the vector containing the elements.
//
// \return The vector containing the elements.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline VT& Elements<VT,TF,true,CEAs...>::operand() noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the vector containing the elements.
//
// \return The vector containing the elements.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline const VT& Elements<VT,TF,true,CEAs...>::operand() const noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the element selection.
//
// \return The minimum capacity of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline size_t Elements<VT,TF,true,CEAs...>::spacing() const noexcept
{
   return size();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the element selection.
//
// \return The maximum capacity of the element selection.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline size_t Elements<VT,TF,true,CEAs...>::capacity() const noexcept
{
   return size();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the element selection.
//
// \return The number of non-zero elements in the element selection.
//
// Note that the number of non-zero elements is always less than or equal to the current number
// of elements.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline size_t Elements<VT,TF,true,CEAs...>::nonZeros() const
{
   size_t nonzeros( 0 );

   for( size_t i=0UL; i<size(); ++i ) {
      if( !isDefault( vector_[idx(i)] ) )
         ++nonzeros;
   }

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline void Elements<VT,TF,true,CEAs...>::reset()
{
   using blaze::clear;

   for( size_t i=0UL; i<size(); ++i )
      clear( vector_[idx(i)] );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the element selection by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the scaling.
// \return Reference to the element selection.
//
// This function scales the element selection by applying the given scalar value \a scalar to each
// element of the element selection. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the scalar value
inline Elements<VT,TF,true,CEAs...>&
   Elements<VT,TF,true,CEAs...>::scale( const Other& scalar )
{
   for( size_t i=0UL; i<size(); ++i )
      vector_[idx(i)] *= scalar;
   return *this;
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the element selection can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this element selection, \a false if not.
//
// This function returns whether the given address can alias with the element selection.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the foreign expression
inline bool Elements<VT,TF,true,CEAs...>::canAlias( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the element selection is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this element selection, \a false if not.
//
// This function returns whether the given address is aliased with the element selection.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the foreign expression
inline bool Elements<VT,TF,true,CEAs...>::isAliased( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the element selection is properly aligned in memory.
//
// \return \a true in case the element selection is aligned, \a false if not.
//
// This function returns whether the element selection is guaranteed to be properly aligned in
// memory, i.e. whether the beginning and the end of the selection are guaranteed to conform to
// the alignment restrictions of the underlying element type.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline bool Elements<VT,TF,true,CEAs...>::isAligned() const noexcept
{
   return false;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the element selection can be used in SMP assignments.
//
// \return \a true in case the element selection can be used in SMP assignments, \a false if not.
//
// This function returns whether the element selection can be used in SMP assignments. In contrast
// to the \a smpAssignable member enumeration, which is based solely on compile time information,
// this function additionally provides runtime information (as for instance the current size of
// the element selection).
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline bool Elements<VT,TF,true,CEAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,true,CEAs...>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[idx(i    )] = (*rhs)[i    ];
      vector_[idx(i+1UL)] = (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[idx(ipos)] = (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Elements<VT,TF,true,CEAs...>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[idx(element->index())] = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,true,CEAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[idx(i    )] += (*rhs)[i    ];
      vector_[idx(i+1UL)] += (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[idx(ipos)] += (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Elements<VT,TF,true,CEAs...>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[idx(element->index())] += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,true,CEAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[idx(i    )] -= (*rhs)[i    ];
      vector_[idx(i+1UL)] -= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[idx(ipos)] -= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Elements<VT,TF,true,CEAs...>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[idx(element->index())] -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,true,CEAs...>::multAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[idx(i    )] *= (*rhs)[i    ];
      vector_[idx(i+1UL)] *= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[idx(ipos)] *= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the multiplication assignment of a sparse vector.
//
// \param rhs The right-hand side sparse vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Elements<VT,TF,true,CEAs...>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( vector_[idx(i)] );
      vector_[idx(i)] *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( vector_[idx(i)] );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,true,CEAs...>::divAssign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[idx(i    )] /= (*rhs)[i    ];
      vector_[idx(i+1UL)] /= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[idx(ipos)] /= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DVECDVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for dense vector/dense vector cross products.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the special case of dense
// vector/dense vector cross products.
*/
template< typename VT1        // Type of the left-hand side dense vector
        , typename VT2        // Type of the right-hand side dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements< const DVecDVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >
   : public View< DenseVector< Elements< const DVecDVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >, TF > >
   , private ElementsData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = DVecDVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = ElementsData<CEAs...>;          //!< The type of the ElementsData base class.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Elements instance.
   using This = Elements<const CPE,TF,true,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = CPE;                          //!< The type viewed by this Elements instance.
   using ResultType    = ElementsTrait_t<RT,N>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;           //!< Type of the elements.
   using ReturnType    = ReturnType_t<CPE>;            //!< Return type for expression template evaluations
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for element selections on dense vector/dense vector cross products.
   //
   // \param vector The dense vector/dense vector cross product expression.
   // \param args The runtime element arguments.
   // \exception std::invalid_argument Invalid element access index.
   */
   template< typename... REAs >  // Optional element arguments
   explicit inline Elements( const CPE& vector, REAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The dense vector/dense vector cross product expression
   {
      if( isChecked( args... ) ) {
         for( size_t i=0UL; i<size(); ++i ) {
            if( vector_.size() <= idx(i) ) {
               BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[idx(index)];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::idces;
   using DataType::idx;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the selection of elements.
   //
   // \return The cross product expression containing the selection of elements.
   */
   inline CPE operand() const noexcept {
      return vector_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return vector_.canAlias( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return vector_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   CPE vector_;  //!< The dense vector/dense vector cross product expression.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR DVECSVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for dense vector/sparse vector cross products.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the special case of dense
// vector/sparse vector cross products.
*/
template< typename VT1        // Type of the left-hand side dense vector
        , typename VT2        // Type of the right-hand side sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements< const DVecSVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >
   : public View< DenseVector< Elements< const DVecSVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >, TF > >
   , private ElementsData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = DVecSVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = ElementsData<CEAs...>;          //!< The type of the ElementsData base class.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Elements instance.
   using This = Elements<const CPE,TF,true,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = CPE;                          //!< The type viewed by this Elements instance.
   using ResultType    = ElementsTrait_t<RT,N>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;           //!< Type of the elements.
   using ReturnType    = ReturnType_t<CPE>;            //!< Return type for expression template evaluations
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for element selections on dense vector/sparse vector cross products.
   //
   // \param vector The dense vector/sparse vector cross product expression.
   // \param args The runtime element arguments.
   // \exception std::invalid_argument Invalid element access index.
   */
   template< typename... REAs >  // Optional element arguments
   explicit inline Elements( const CPE& vector, REAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The dense vector/sparse vector cross product expression
   {
      if( isChecked( args... ) ) {
         for( size_t i=0UL; i<size(); ++i ) {
            if( vector_.size() <= idx(i) ) {
               BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[idx(index)];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::idces;
   using DataType::idx;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the selection of elements.
   //
   // \return The cross product expression containing the selection of elements.
   */
   inline CPE operand() const noexcept {
      return vector_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return vector_.canAlias( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return vector_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   CPE vector_;  //!< The dense vector/sparse vector cross product expression.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SVECDVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for sparse vector/dense vector cross products.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the special case of sparse
// vector/dense vector cross products.
*/
template< typename VT1        // Type of the left-hand side sparse vector
        , typename VT2        // Type of the right-hand side dense vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements< const SVecDVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >
   : public View< DenseVector< Elements< const SVecDVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >, TF > >
   , private ElementsData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = SVecDVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = ElementsData<CEAs...>;          //!< The type of the ElementsData base class.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Elements instance.
   using This = Elements<const CPE,TF,true,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = CPE;                          //!< The type viewed by this Elements instance.
   using ResultType    = ElementsTrait_t<RT,N>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;           //!< Type of the elements.
   using ReturnType    = ReturnType_t<CPE>;            //!< Return type for expression template evaluations
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for element selections on sparse vector/dense vector cross products.
   //
   // \param vector The sparse vector/dense vector cross product expression.
   // \param args The runtime element arguments.
   // \exception std::invalid_argument Invalid element access index.
   */
   template< typename... REAs >  // Optional element arguments
   explicit inline Elements( const CPE& vector, REAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The sparse vector/dense vector cross product expression
   {
      if( isChecked( args... ) ) {
         for( size_t i=0UL; i<size(); ++i ) {
            if( vector_.size() <= idx(i) ) {
               BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[idx(index)];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::idces;
   using DataType::idx;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the selection of elements.
   //
   // \return The cross product expression containing the selection of elements.
   */
   inline CPE operand() const noexcept {
      return vector_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return vector_.canAlias( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return vector_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   CPE vector_;  //!< The sparse vector/dense vector cross product expression.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SVECSVECCROSSEXPR
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for sparse vector/sparse vector cross products.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the special case of sparse
// vector/sparse vector cross products.
*/
template< typename VT1        // Type of the left-hand side sparse vector
        , typename VT2        // Type of the right-hand side sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements< const SVecSVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >
   : public View< DenseVector< Elements< const SVecSVecCrossExpr<VT1,VT2,TF>, TF, true, CEAs... >, TF > >
   , private ElementsData<CEAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = SVecSVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;               //!< Result type of the cross product expression.
   using DataType = ElementsData<CEAs...>;          //!< The type of the ElementsData base class.
   //**********************************************************************************************

   //**Compile time flags**************************************************************************
   using DataType::N;  //! Number of compile time indices.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Elements instance.
   using This = Elements<const CPE,TF,true,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< DenseVector<This,TF> >;

   using ViewedType    = CPE;                          //!< The type viewed by this Elements instance.
   using ResultType    = ElementsTrait_t<RT,N>;        //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;  //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;           //!< Type of the elements.
   using ReturnType    = ReturnType_t<CPE>;            //!< Return type for expression template evaluations
   using CompositeType = const ResultType;             //!< Data type for composite expression templates.
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = false;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructor*********************************************************************************
   /*!\brief Constructor for element selections on sparse vector/sparse vector cross products.
   //
   // \param vector The sparse vector/sparse vector cross product expression.
   // \param args The runtime element arguments.
   // \exception std::invalid_argument Invalid element access index.
   */
   template< typename... REAs >  // Optional element arguments
   explicit inline Elements( const CPE& vector, REAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The sparse vector/sparse vector cross product expression
   {
      if( isChecked( args... ) ) {
         for( size_t i=0UL; i<size(); ++i ) {
            if( vector_.size() <= idx(i) ) {
               BLAZE_THROW_INVALID_ARGUMENT( "Invalid element access index" );
            }
         }
      }
   }
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[idx(index)];
   }
   //**********************************************************************************************

   //**At function*********************************************************************************
   /*!\brief Checked access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   // \exception std::out_of_range Invalid vector access index.
   */
   inline ReturnType at( size_t index ) const {
      if( index >= size() ) {
         BLAZE_THROW_OUT_OF_RANGE( "Invalid vector access index" );
      }
      return (*this)[index];
   }
   //**********************************************************************************************

   //**********************************************************************************************
   using DataType::idces;
   using DataType::idx;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the selection of elements.
   //
   // \return The cross product expression containing the selection of elements.
   */
   inline CPE operand() const noexcept {
      return vector_;
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression can alias with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case the expression can alias, \a false otherwise.
   */
   template< typename T >
   inline bool canAlias( const T* alias ) const noexcept {
      return vector_.canAlias( &unview( *alias ) );
   }
   //**********************************************************************************************

   //**********************************************************************************************
   /*!\brief Returns whether the expression is aliased with the given address \a alias.
   //
   // \param alias The alias to be checked.
   // \return \a true in case an alias effect is detected, \a false otherwise.
   */
   template< typename T >
   inline bool isAliased( const T* alias ) const noexcept {
      return vector_.isAliased( &unview( *alias ) );
   }
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   CPE vector_;  //!< The sparse vector/sparse vector cross product expression.
   //**********************************************************************************************
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
