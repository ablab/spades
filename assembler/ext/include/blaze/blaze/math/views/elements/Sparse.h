//=================================================================================================
/*!
//  \file blaze/math/views/elements/Sparse.h
//  \brief Elements specialization for sparse vectors
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_SPARSE_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/Elements.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Subvector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/typetraits/IsComputation.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/elements/BaseTemplate.h>
#include <blaze/math/views/elements/ElementsData.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/util/Assert.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE VECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Elements for sparse vectors.
// \ingroup elements
//
// This specialization of Elements adapts the class template to the requirements of sparse vectors.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
class Elements<VT,TF,false,CEAs...>
   : public View< SparseVector< Elements<VT,TF,false,CEAs...>, TF > >
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
   using This = Elements<VT,TF,false,CEAs...>;

   //! Base type of this Elements instance.
   using BaseType = View< SparseVector<This,TF> >;

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
   //**********************************************************************************************

   //**ElementsElement class definition************************************************************
   /*!\brief Access proxy for a specific element of the selection of elements.
   */
   template< typename IteratorType >  // Type of the sparse vector iterator
   class ElementsElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the ElementsElement class.
      //
      // \param pos Iterator to the current position within the element selection.
      // \param offset The offset within the according sparse vector.
      */
      inline ElementsElement( IteratorType pos, size_t index )
         : pos_  ( pos   )  // Iterator to the current position within the element selection
         , index_( index )  // Index within the according element selection
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the selected sparse element.
      //
      // \param v The new value of the selected sparse element.
      // \return Reference to the selected sparse element.
      */
      template< typename T > inline ElementsElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the selected sparse element.
      //
      // \param v The right-hand side value for the addition.
      // \return Reference to the selected sparse element.
      */
      template< typename T > inline ElementsElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the selected sparse element.
      //
      // \param v The right-hand side value for the subtraction.
      // \return Reference to the selected sparse element.
      */
      template< typename T > inline ElementsElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the selected sparse element.
      //
      // \param v The right-hand side value for the multiplication.
      // \return Reference to the selected sparse element.
      */
      template< typename T > inline ElementsElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the selected sparse element.
      //
      // \param v The right-hand side value for the division.
      // \return Reference to the selected sparse element.
      */
      template< typename T > inline ElementsElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the selected sparse element at the current iterator position.
      //
      // \return Reference to the selected sparse element at the current iterator position.
      */
      inline const ElementsElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse selected element.
      //
      // \return The current value of the sparse selected element.
      */
      inline decltype(auto) value() const {
         return pos_->value();
      }
      //*******************************************************************************************

      //**Index function***************************************************************************
      /*!\brief Access to the current index of the sparse element.
      //
      // \return The current index of the sparse element.
      */
      inline size_t index() const {
         return index_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the element selection.
      size_t index_;      //!< Index within the according element selection.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**ElementsIterator class definition***********************************************************
   /*!\brief Iterator over the elements of the selection of elements.
   */
   template< typename ET              // Type of the element selection
           , typename IteratorType >  // Type of the sparse vector iterator
   class ElementsIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;      //!< The iterator category.
      using ValueType        = ElementsElement<IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                      //!< Pointer return type.
      using ReferenceType    = ValueType;                      //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                      //!< Difference between two iterators.

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
         , index_   ( 0UL)       // Index within the according element selection
         , pos_     ()           // Iterator to the current position within the element selection
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ElementsIterator class.
      //
      // \param elements Pointer to the element selection.
      // \param index Index of the initial element.
      */
      inline ElementsIterator( ET* elements, size_t index )
         : elements_( elements )  // Pointer to the element selection
         , index_   ( index )     // Index within the according element selection
         , pos_     ()            // Iterator to the current position within the element selection
      {
         for( ; index_<elements_->size(); ++index_ ) {
            pos_ = elements_->operand().find( elements_->idx(index_) );
            if( pos_ != elements_->operand().end() ) break;
         }
      }
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the ElementsIterator class.
      //
      // \param elements Pointer to the element selection.
      // \param index Index of the initial element.
      // \param pos Initial position of the iterator.
      */
      inline ElementsIterator( ET* elements, size_t index, IteratorType pos )
         : elements_( elements )  // Pointer to the element selection
         , index_   ( index )     // Index within the according element selection
         , pos_     ( pos )       // Iterator to the current position within the element selection
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different ElementsIterator instances.
      //
      // \param it The elements iterator to be copied.
      */
      template< typename ET2, typename IteratorType2 >
      inline ElementsIterator( const ElementsIterator<ET2,IteratorType2>& it )
         : elements_( it.elements_ )  // Pointer to the element selection
         , index_   ( it.index_ )     // Index within the according element selection
         , pos_     ( it.pos_ )       // Iterator to the current position within the element selection
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline ElementsIterator& operator++() {
         ++index_;
         for( ; index_<elements_->size(); ++index_ ) {
            pos_ = elements_->operand().find( elements_->idx(index_) );
            if( pos_ != elements_->operand().end() ) break;
         }
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

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse subvector element.
      //
      // \return Reference to the sparse subvector element.
      */
      inline ReferenceType operator*() const {
         return ReferenceType( pos_, index_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse subvector element.
      //
      // \return Pointer to the sparse subvector element.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, index_ );
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

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const ElementsIterator& rhs ) const {
         size_t counter( 0UL );
         for( size_t j=rhs.index_; j<index_; ++j ) {
            if( elements_->find( j ) != elements_->end() )
               ++counter;
         }
         return counter;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      ET*          elements_;  //!< Pointer to the element selection.
      size_t       index_;     //!< Index within the according element selection.
      IteratorType pos_;       //!< Iterator to the current element of the element selection.
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      template< typename VT2, bool TF2, bool DF2, typename... CEAs2 > friend class Elements;
      template< typename ET2, typename IteratorType2 > friend class ElementsIterator;
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
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = false;

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

   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   reserve( size_t n );
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const ElementType& value );
   inline Iterator insert( size_t index, const ElementType& value );
   inline void     append( size_t index, const ElementType& value, bool check=false );
   //@}
   //**********************************************************************************************

   //**Erase functions*****************************************************************************
   /*!\name Erase functions */
   //@{
   inline void     erase( size_t index );
   inline Iterator erase( Iterator pos );
   inline Iterator erase( Iterator first, Iterator last );

   template< typename Pred, typename = DisableIf_t< IsIntegral_v<Pred> > >
   inline void erase( Pred predicate );

   template< typename Pred >
   inline void erase( Iterator first, Iterator last, Pred predicate );
   //@}
   //**********************************************************************************************

   //**Lookup functions****************************************************************************
   /*!\name Lookup functions */
   //@{
   inline Iterator      find      ( size_t index );
   inline ConstIterator find      ( size_t index ) const;
   inline Iterator      lowerBound( size_t index );
   inline ConstIterator lowerBound( size_t index ) const;
   inline Iterator      upperBound( size_t index );
   inline ConstIterator upperBound( size_t index ) const;
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
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename VT2 >   inline void assign   ( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void assign   ( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 >   inline void addAssign( const Vector<VT2,TF>& rhs );
   template< typename VT2 >   inline void subAssign( const Vector<VT2,TF>& rhs );
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
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE  ( VT );
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
/*!\brief Constructor for elements on sparse vectors.
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename... REAs >  // Optional arguments
inline Elements<VT,TF,false,CEAs...>::Elements( VT& vector, REAs... args )
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Reference
   Elements<VT,TF,false,CEAs...>::operator[]( size_t index )
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstReference
   Elements<VT,TF,false,CEAs...>::operator[]( size_t index ) const
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Reference
   Elements<VT,TF,false,CEAs...>::at( size_t index )
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstReference
   Elements<VT,TF,false,CEAs...>::at( size_t index ) const
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
/*!\brief Returns an iterator to the first element of the element selection.
//
// \return Iterator to the first element of the element selection.
//
// This function returns an iterator to the first element of the element selection.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::begin()
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::begin() const
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::cbegin() const
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::end()
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::end() const
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::cend() const
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

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

   decltype(auto) left( derestrict( *this ) );

   assign( left, tmp );

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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator=( const Elements& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
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
      assign( left, tmp );
   }
   else {
      assign( left, rhs );
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

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

   if( IsReference_v<Right> || right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      assign( left, tmp );
   }
   else {
      assign( left, right );
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
// \param rhs The right-hand side vector to be added to the sparse element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator+=( const Vector<VT2,TF>& rhs )
{
   using blaze::addAssign;

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
      addAssign( left, tmp );
   }
   else {
      addAssign( left, right );
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
// \param rhs The right-hand side vector to be subtracted from the sparse element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator-=( const Vector<VT2,TF>& rhs )
{
   using blaze::subAssign;

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
      subAssign( left, tmp );
   }
   else {
      subAssign( left, right );
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
// \param rhs The right-hand side vector to be multiplied with the sparse element selection.
// \return Reference to the assigned element selection.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator*=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   using MultType = MultTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( MultType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( MultType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const MultType tmp( *this * (*rhs) );

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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator/=( const DenseVector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType_t<VT2> );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   using DivType = DivTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( DivType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( DivType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( DivType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const DivType tmp( *this / (*rhs) );

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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::operator%=( const Vector<VT2,TF>& rhs )
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline VT& Elements<VT,TF,false,CEAs...>::operand() noexcept
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline const VT& Elements<VT,TF,false,CEAs...>::operand() const noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the element selection.
//
// \return The maximum capacity of the element selection.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline size_t Elements<VT,TF,false,CEAs...>::capacity() const noexcept
{
   return nonZeros() + vector_.capacity() - vector_.nonZeros();
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline size_t Elements<VT,TF,false,CEAs...>::nonZeros() const
{
   size_t counter( 0UL );
   for( ConstIterator element=begin(); element!=end(); ++element ) {
      ++counter;
   }
   return counter;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline void Elements<VT,TF,false,CEAs...>::reset()
{
   for( size_t i=0UL; i<size(); ++i )
      vector_.erase( idx(i) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the element selection.
//
// \param n The new minimum capacity of the element selection.
// \return void
//
// This function increases the capacity of the element selection to at least \a n elements. The
// current values of the elements are preserved.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
void Elements<VT,TF,false,CEAs...>::reserve( size_t n )
{
   const size_t current( capacity() );

   if( n > current ) {
      vector_.reserve( vector_.capacity() + n - current );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting an element of the element selection.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the element selection. In case the element
// selection already contains an element with index \a index its value is modified, else a new
// element with the given \a value is inserted.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::set( size_t index, const ElementType& value )
{
   return Iterator( this, index, vector_.set( idx(index), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the element selection.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid element access index.
//
// This function inserts a new element into the element selection. However, duplicate elements
// are not allowed. In case the element selection already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::insert( size_t index, const ElementType& value )
{
   return Iterator( this, index, vector_.insert( idx(index), value ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the element selection.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill an element selection with elements. It
// appends a new element to the end of the element selection without any memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the element selection
//  - the current number of non-zero elements must be smaller than the capacity of the selection
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline void Elements<VT,TF,false,CEAs...>::append( size_t index, const ElementType& value, bool check )
{
   if( !check || !isDefault<strict>( value ) )
      vector_.insert( idx(index), value );
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the element selection.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the element selection.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline void Elements<VT,TF,false,CEAs...>::erase( size_t index )
{
   vector_.erase( idx(index) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the element selection.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the element selection.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::erase( Iterator pos )
{
   const size_t index( pos.index_ );

   if( index == size() )
      return pos;

   vector_.erase( pos.pos_ );
   return Iterator( this, index+1UL );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the element selection.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the element selection.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::erase( Iterator first, Iterator last )
{
   for( ; first!=last; ++first ) {
      vector_.erase( first.pos_ );
   }

   return Iterator( this, last.index_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from the sparse subvector.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the sparse subvector. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove
// all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedVector<double,blaze::rowVector> a;
   // ... Resizing and initialization

   auto sv = subvector( a, 5UL, 7UL );
   sv.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Elements<VT,TF,false,CEAs...>::erase( Pred predicate )
{
   erase( begin(), end(), predicate );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing specific elements from a range of the sparse subvector.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the sparse subvector.
// The elements are selected by the given unary predicate \a predicate, which is expected
// to accept a single argument of the type of the elements and to be pure. The following
// example demonstrates how to remove all elements that are smaller than a certain threshold
// value:

   \code
   blaze::CompressedVector<double,blaze::rowVector> a;
   // ... Resizing and initialization

   auto sv = subvector( a, 5UL, 7UL );
   sv.erase( sv.begin(), sv.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Pred >   // Type of the unary predicate
inline void Elements<VT,TF,false,CEAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   for( ; first!=last; ++first ) {
      if( predicate( first->value() ) )
         vector_.erase( first.pos_ );
   }
}
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific element in the element selection.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the element
// selection. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the element selection (the end() iterator) is returned. Note that
// the returned iterator is subject to invalidation due to inserting operations via the subscript
// operator, the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::find( size_t index )
{
   const Iterator_t<VT> pos( vector_.find( idx(index) ) );

   if( pos != vector_.end() )
      return Iterator( this, index, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific element in the element selection.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the element
// selection. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the element selection (the end() iterator) is returned. Note that
// the returned iterator is subject to invalidation due to inserting operations via the subscript
// operator, the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::find( size_t index ) const
{
   const ConstIterator_t<VT> pos( vector_.find( idx(index) ) );

   if( pos != vector_.end() )
      return ConstIterator( this, index, pos );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::lowerBound( size_t index )
{
   for( ; index<size(); ++index ) {
      const Iterator_t<VT> pos( vector_.find( idx(index) ) );
      if( pos != vector_.end() )
         return Iterator( this, index, pos );
   }
   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::lowerBound( size_t index ) const
{
   for( ; index<size(); ++index ) {
      const ConstIterator_t<VT> pos( vector_.find( idx(index) ) );
      if( pos != vector_.end() )
         return ConstIterator( this, index, pos );
   }
   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::Iterator
   Elements<VT,TF,false,CEAs...>::upperBound( size_t index )
{
   while( (++index) < size() ) {
      const Iterator_t<VT> pos( vector_.find( idx(index) ) );
      if( pos != vector_.end() )
         return Iterator( this, index, pos );
   }
   return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned sparse subvector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
inline typename Elements<VT,TF,false,CEAs...>::ConstIterator
   Elements<VT,TF,false,CEAs...>::upperBound( size_t index ) const
{
   while( (++index) < size() ) {
      const ConstIterator_t<VT> pos( vector_.find( idx(index) ) );
      if( pos != vector_.end() )
         return ConstIterator( this, index, pos );
   }
   return end();
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the scalar value
inline Elements<VT,TF,false,CEAs...>&
   Elements<VT,TF,false,CEAs...>::scale( const Other& scalar )
{
   for( Iterator element=begin(); element!=end(); ++element )
      element->value() *= scalar;
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the foreign expression
inline bool Elements<VT,TF,false,CEAs...>::canAlias( const Other* alias ) const noexcept
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename Other >  // Data type of the foreign expression
inline bool Elements<VT,TF,false,CEAs...>::isAliased( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Elements<VT,TF,false,CEAs...>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   using RT = If_t< IsComputation_v<VT2>, ElementType_t<VT>, const ElementType_t<VT2>& >;

   reserve( (*rhs).size() );

   for( size_t i=0UL; i<size(); ++i ) {
      RT value( (*rhs)[i] );
      if( !isDefault<strict>( value ) )
         vector_.set( idx(i), std::move( value ) );
      else vector_.erase( idx(i) );
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
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Elements<VT,TF,false,CEAs...>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   using RT = If_t< IsComputation_v<VT2>, ElementType_t<VT>, const ElementType_t<VT2>& >;

   size_t i( 0UL );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      for( ; i<element->index(); ++i )
         vector_.erase( idx(i) );
      RT value( element->value() );
      if( !isDefault<strict>( value ) )
         vector_.set( idx(i), std::move( value ) );
      else vector_.erase( idx(i) );
      ++i;
   }
   for( ; i<size(); ++i ) {
      vector_.erase( idx(i) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a vector.
//
// \param rhs The right-hand side vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline void Elements<VT,TF,false,CEAs...>::addAssign( const Vector<VT2,TF>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a vector.
//
// \param rhs The right-hand side vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT         // Type of the sparse vector
        , bool TF             // Transpose flag
        , typename... CEAs >  // Compile time element arguments
template< typename VT2 >    // Type of the right-hand side vector
inline void Elements<VT,TF,false,CEAs...>::subAssign( const Vector<VT2,TF>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
