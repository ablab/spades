//=================================================================================================
/*!
//  \file blaze/math/views/subvector/Sparse.h
//  \brief Subvector specialization for sparse vectors
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_SPARSE_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_SPARSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SparseVector.h>
#include <blaze/math/constraints/Subvector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/subvector/BaseTemplate.h>
#include <blaze/math/views/subvector/SubvectorData.h>
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
//  CLASS TEMPLATE SPECIALIZATION FOR SPARSE SUBVECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Subvector for sparse subvectors.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the requirements of sparse
// subvectors.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector<VT,AF,TF,false,CSAs...>
   : public View< SparseVector< Subvector<VT,AF,TF,false,CSAs...>, TF > >
   , private SubvectorData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = SubvectorData<CSAs...>;               //!< The type of the SubvectorData base class.
   using Operand  = If_t< IsExpression_v<VT>, VT, VT& >;  //!< Composite data type of the vector expression.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Subvector instance.
   using This = Subvector<VT,AF,TF,false,CSAs...>;

   using BaseType      = View< SparseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = VT;                             //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<VT,CSAs...>;   //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;    //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;              //!< Type of the subvector elements.
   using ReturnType    = ReturnType_t<VT>;               //!< Return type for expression template evaluations
   using CompositeType = const Subvector&;               //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   using ConstReference = ConstReference_t<VT>;

   //! Reference to a non-constant subvector value.
   using Reference = If_t< IsConst_v<VT>, ConstReference, Reference_t<VT> >;
   //**********************************************************************************************

   //**SubvectorElement class definition***********************************************************
   /*!\brief Access proxy for a specific element of the sparse subvector.
   */
   template< typename VectorType      // Type of the sparse vector
           , typename IteratorType >  // Type of the sparse vector iterator
   class SubvectorElement
      : private SparseElement
   {
    public:
      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubvectorElement class.
      //
      // \param pos Iterator to the current position within the sparse subvector.
      // \param offset The offset within the according sparse vector.
      */
      inline SubvectorElement( IteratorType pos, size_t offset )
         : pos_   ( pos    )  // Iterator to the current position within the sparse subvector
         , offset_( offset )  // Offset within the according sparse vector
      {}
      //*******************************************************************************************

      //**Assignment operator**********************************************************************
      /*!\brief Assignment to the accessed sparse subvector element.
      //
      // \param v The new value of the sparse subvector element.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator=( const T& v ) {
         *pos_ = v;
         return *this;
      }
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the addition.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator+=( const T& v ) {
         *pos_ += v;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the subtraction.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator-=( const T& v ) {
         *pos_ -= v;
         return *this;
      }
      //*******************************************************************************************

      //**Multiplication assignment operator*******************************************************
      /*!\brief Multiplication assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the multiplication.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator*=( const T& v ) {
         *pos_ *= v;
         return *this;
      }
      //*******************************************************************************************

      //**Division assignment operator*************************************************************
      /*!\brief Division assignment to the accessed sparse subvector element.
      //
      // \param v The right-hand side value for the division.
      // \return Reference to the sparse subvector element.
      */
      template< typename T > inline SubvectorElement& operator/=( const T& v ) {
         *pos_ /= v;
         return *this;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the sparse subvector element at the current iterator position.
      //
      // \return Reference to the sparse subvector element at the current iterator position.
      */
      inline const SubvectorElement* operator->() const {
         return this;
      }
      //*******************************************************************************************

      //**Value function***************************************************************************
      /*!\brief Access to the current value of the sparse subvector element.
      //
      // \return The current value of the sparse subvector element.
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
         return pos_->index() - offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;  //!< Iterator to the current position within the sparse subvector.
      size_t offset_;     //!< Offset within the according sparse vector.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**SubvectorIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse subvector.
   */
   template< typename VectorType      // Type of the sparse vector
           , typename IteratorType >  // Type of the sparse vector iterator
   class SubvectorIterator
   {
    public:
      //**Type definitions*************************************************************************
      using IteratorCategory = std::forward_iterator_tag;                  //!< The iterator category.
      using ValueType        = SubvectorElement<VectorType,IteratorType>;  //!< Type of the underlying elements.
      using PointerType      = ValueType;                                  //!< Pointer return type.
      using ReferenceType    = ValueType;                                  //!< Reference return type.
      using DifferenceType   = ptrdiff_t;                                  //!< Difference between two iterators.

      // STL iterator requirements
      using iterator_category = IteratorCategory;  //!< The iterator category.
      using value_type        = ValueType;         //!< Type of the underlying elements.
      using pointer           = PointerType;       //!< Pointer return type.
      using reference         = ReferenceType;     //!< Reference return type.
      using difference_type   = DifferenceType;    //!< Difference between two iterators.
      //*******************************************************************************************

      //**Default constructor**********************************************************************
      /*!\brief Default constructor for the SubvectorIterator class.
      */
      inline SubvectorIterator()
         : pos_   ()  // Iterator to the current sparse element
         , offset_()  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor for the SubvectorIterator class.
      //
      // \param iterator Iterator to the current sparse element.
      // \param index The starting index of the subvector within the sparse vector.
      */
      inline SubvectorIterator( IteratorType iterator, size_t index )
         : pos_   ( iterator )  // Iterator to the current sparse element
         , offset_( index    )  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubvectorIterator instances.
      //
      // \param it The subvector iterator to be copied.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline SubvectorIterator( const SubvectorIterator<VectorType2,IteratorType2>& it )
         : pos_   ( it.base()   )  // Iterator to the current sparse element.
         , offset_( it.offset() )  // The offset of the subvector within the sparse vector
      {}
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubvectorIterator& operator++() {
         ++pos_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubvectorIterator operator++( int ) {
         const SubvectorIterator tmp( *this );
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
         return ReferenceType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the current sparse subvector element.
      //
      // \return Pointer to the sparse subvector element.
      */
      inline PointerType operator->() const {
         return PointerType( pos_, offset_ );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline bool operator==( const SubvectorIterator<VectorType2,IteratorType2>& rhs ) const {
         return base() == rhs.base();
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      template< typename VectorType2, typename IteratorType2 >
      inline bool operator!=( const SubvectorIterator<VectorType2,IteratorType2>& rhs ) const {
         return !( *this == rhs );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two subvector iterators.
      //
      // \param rhs The right-hand side subvector iterator.
      // \return The number of elements between the two subvector iterators.
      */
      inline DifferenceType operator-( const SubvectorIterator& rhs ) const {
         return pos_ - rhs.pos_;
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the subvector iterator.
      //
      // \return The current position of the subvector iterator.
      */
      inline IteratorType base() const {
         return pos_;
      }
      //*******************************************************************************************

      //**Offset function**************************************************************************
      /*!\brief Access to the offset of the subvector iterator.
      //
      // \return The offset of the subvector iterator.
      */
      inline size_t offset() const noexcept {
         return offset_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType pos_;     //!< Iterator to the current sparse element.
      size_t       offset_;  //!< The offset of the subvector within the sparse vector.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = SubvectorIterator< const VT, ConstIterator_t<VT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<VT>, ConstIterator, SubvectorIterator< VT, Iterator_t<VT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = VT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RSAs >
   explicit inline Subvector( VT& vector, RSAs... args );

   Subvector( const Subvector& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Subvector() = default;
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
                            inline Subvector& operator= ( initializer_list<ElementType> list );
                            inline Subvector& operator= ( const Subvector& rhs );
   template< typename VT2 > inline Subvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator*=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator/=( const DenseVector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator%=( const Vector<VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::offset;
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
   template< typename Other > inline Subvector& scale( const Other& scalar );
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
   template< typename VT2 >   inline void addAssign( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void addAssign( const SparseVector<VT2,TF>& rhs );
   template< typename VT2 >   inline void subAssign( const DenseVector <VT2,TF>& rhs );
   template< typename VT2 >   inline void subAssign( const SparseVector<VT2,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand vector_;  //!< The vector containing the subvector.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE  ( VT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBVECTOR_TYPE  ( VT );
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
/*!\brief Constructor for sparse subvectors.
//
// \param vector The sparse vector containing the subvector.
// \param args The runtime subvector arguments.
// \exception std::invalid_argument Invalid subvector specification.
//
// By default, the provided subvector arguments are checked at runtime. In case the subvector is
// not properly specified (i.e. if the specified offset is greater than the size of the given
// vector or the subvector is specified beyond the size of the vector) a \a std::invalid_argument
// exception is thrown. The checks can be skipped by providing the optional \a blaze::unchecked
// argument.
*/
template< typename VT         // Type of the sparse vector
        , AlignmentFlag AF    // Alignment flag
        , bool TF             // Transpose flag
        , size_t... CSAs >    // Compile time subvector arguments
template< typename... RSAs >  // Runtime subvector arguments
inline Subvector<VT,AF,TF,false,CSAs...>::Subvector( VT& vector, RSAs... args )
   : DataType( args... )  // Base class initialization
   , vector_ ( vector  )  // The vector containing the subvector
{
   if( isChecked( args... ) ) {
      if( offset() + size() > vector.size() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( offset() + size() <= vector.size(), "Invalid subvector specification" );
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
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Reference
   Subvector<VT,AF,TF,false,CSAs...>::operator[]( size_t index )
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return vector_[offset()+index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subscript operator for the direct access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access index.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstReference
   Subvector<VT,AF,TF,false,CSAs...>::operator[]( size_t index ) const
{
   BLAZE_USER_ASSERT( index < size(), "Invalid subvector access index" );
   return const_cast<const VT&>( vector_ )[offset()+index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Reference
   Subvector<VT,AF,TF,false,CSAs...>::at( size_t index )
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the subvector elements.
//
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstReference
   Subvector<VT,AF,TF,false,CSAs...>::at( size_t index ) const
{
   if( index >= size() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid subvector access index" );
   }
   return (*this)[index];
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::begin()
{
   if( offset() == 0UL )
      return Iterator( vector_.begin(), offset() );
   else
      return Iterator( vector_.lowerBound( offset() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::begin() const
{
   if( offset() == 0UL )
      return ConstIterator( vector_.cbegin(), offset() );
   else
      return ConstIterator( vector_.lowerBound( offset() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first element of the subvector.
//
// \return Iterator to the first element of the subvector.
//
// This function returns an iterator to the first element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::cbegin() const
{
   if( offset() == 0UL )
      return ConstIterator( vector_.cbegin(), offset() );
   else
      return ConstIterator( vector_.lowerBound( offset() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::end()
{
   if( offset() + size() == vector_.size() )
      return Iterator( vector_.end(), offset() );
   else
      return Iterator( vector_.lowerBound( offset() + size() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::end() const
{
   if( offset() + size() == vector_.size() )
      return ConstIterator( vector_.cend(), offset() );
   else
      return ConstIterator( vector_.lowerBound( offset() + size() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last element of the subvector.
//
// \return Iterator just past the last element of the subvector.
//
// This function returns an iterator just past the last element of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::cend() const
{
   if( offset() + size() == vector_.size() )
      return ConstIterator( vector_.cend(), offset() );
   else
      return ConstIterator( vector_.lowerBound( offset() + size() ), offset() );
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
/*!\brief List assignment to all subvector elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to subvector.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// This assignment operator offers the option to directly assign to all elements of the subvector
// by means of an initializer list. The subvector elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the subvector, a \a std::invalid_argument exception
// is thrown. Also, if the underlying vector \a VT is restricted and the assignment would violate
// an invariant of the vector, a \a std::invalid_argument exception is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator=( initializer_list<ElementType> list )
{
   using blaze::assign;

   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to subvector" );
   }

   const InitializerVector<ElementType,TF> tmp( list, size() );

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Subvector.
//
// \param rhs Sparse subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator=( const Subvector& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &vector_ == &rhs.vector_ && offset() == rhs.offset() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      reset();
      assign( left, tmp );
   }
   else {
      reset();
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> || right.canAlias( this ) ) {
      const ResultType_t<VT2> tmp( right );
      reset();
      assign( left, tmp );
   }
   else {
      reset();
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
// \param rhs The right-hand side vector to be added to the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator+=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator-=( const Vector<VT2,TF>& rhs )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
   assign( left, tmp );

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
// \param rhs The right-hand side vector to be multiplied with the sparse subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator*=( const Vector<VT2,TF>& rhs )
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

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator/=( const DenseVector<VT2,TF>& rhs )
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

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::operator%=( const Vector<VT2,TF>& rhs )
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

   if( !tryAssign( vector_, tmp, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   left.reset();
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
/*!\brief Returns the vector containing the subvector.
//
// \return The vector containing the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline VT& Subvector<VT,AF,TF,false,CSAs...>::operand() noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the vector containing the subvector.
//
// \return The vector containing the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline const VT& Subvector<VT,AF,TF,false,CSAs...>::operand() const noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the sparse subvector.
//
// \return The capacity of the sparse subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,AF,TF,false,CSAs...>::capacity() const noexcept
{
   return nonZeros() + vector_.capacity() - vector_.nonZeros();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always smaller than the size of the subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,AF,TF,false,CSAs...>::nonZeros() const
{
   return end() - begin();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline void Subvector<VT,AF,TF,false,CSAs...>::reset()
{
   vector_.erase( vector_.lowerBound( offset() ), vector_.lowerBound( offset() + size() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Setting the minimum capacity of the sparse subvector.
//
// \param n The new minimum capacity of the sparse subvector.
// \return void
//
// This function increases the capacity of the sparse subvector to at least \a n elements. The
// current values of the subvector elements are preserved.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
void Subvector<VT,AF,TF,false,CSAs...>::reserve( size_t n )
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
/*!\brief Setting an element of the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the sparse subvector. In case the sparse subvector
// already contains an element with index \a index its value is modified, else a new element with
// the given \a value is inserted.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::set( size_t index, const ElementType& value )
{
   return Iterator( vector_.set( offset() + index, value ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Inserting an element into the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid sparse subvector access index.
//
// This function inserts a new element into the sparse subvector. However, duplicate elements
// are not allowed. In case the sparse subvector already contains an element at index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::insert( size_t index, const ElementType& value )
{
   return Iterator( vector_.insert( offset() + index, value ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Appending an element to the sparse subvector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a sparse subvector with elements. It
// appends a new element to the end of the sparse subvector without any memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the sparse subvector
//  - the current number of non-zero elements must be smaller than the capacity of the subvector
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline void Subvector<VT,AF,TF,false,CSAs...>::append( size_t index, const ElementType& value, bool check )
{
   if( offset() + size() == vector_.size() )
      vector_.append( offset() + index, value, check );
   else if( !check || !isDefault<strict>( value ) )
      vector_.insert( offset() + index, value );
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
/*!\brief Erasing an element from the sparse subvector.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the sparse subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline void Subvector<VT,AF,TF,false,CSAs...>::erase( size_t index )
{
   vector_.erase( offset() + index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing an element from the sparse subvector.
//
// \param pos Iterator to the element to be erased.
// \return void
//
// This function erases an element from the sparse subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::erase( Iterator pos )
{
   return Iterator( vector_.erase( pos.base() ), offset() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Erasing a range of elements from the sparse subvector.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the sparse subvector.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::erase( Iterator first, Iterator last )
{
   return Iterator( vector_.erase( first.base(), last.base() ), offset() );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Pred     // Type of the unary predicate
        , typename >        // Type restriction on the unary predicate
inline void Subvector<VT,AF,TF,false,CSAs...>::erase( Pred predicate )
{
   vector_.erase( begin().base(), end().base(), predicate );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Pred >   // Type of the unary predicate
inline void Subvector<VT,AF,TF,false,CSAs...>::erase( Iterator first, Iterator last, Pred predicate )
{
   vector_.erase( first.base(), last.base(), predicate );
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
/*!\brief Searches for a specific subvector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// subvector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse subvector (the end() iterator) is returned. Note that
// the returned sparse subvector iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::find( size_t index )
{
   const Iterator_t<VT> pos( vector_.find( offset() + index ) );

   if( pos != vector_.end() )
      return Iterator( pos, offset() );
   else
      return end();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Searches for a specific subvector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// subvector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the sparse subvector (the end() iterator) is returned. Note that
// the returned sparse subvector iterator is subject to invalidation due to inserting operations
// via the subscript operator, the set() function or the insert() function!
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::find( size_t index ) const
{
   const ConstIterator_t<VT> pos( vector_.find( offset() + index ) );

   if( pos != vector_.end() )
      return Iterator( pos, offset() );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::lowerBound( size_t index )
{
   return Iterator( vector_.lowerBound( offset() + index ), offset() );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::lowerBound( size_t index ) const
{
   return ConstIterator( vector_.lowerBound( offset() + index ), offset() );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::Iterator
   Subvector<VT,AF,TF,false,CSAs...>::upperBound( size_t index )
{
   return Iterator( vector_.upperBound( offset() + index ), offset() );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,AF,TF,false,CSAs...>::ConstIterator
   Subvector<VT,AF,TF,false,CSAs...>::upperBound( size_t index ) const
{
   return ConstIterator( vector_.upperBound( offset() + index ), offset() );
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
/*!\brief Scaling of the sparse subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the sparse subvector.
//
// This function scales the subvector by applying the given scalar value \a scalar to each
// element of the subvector. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the scalar value
inline Subvector<VT,AF,TF,false,CSAs...>&
   Subvector<VT,AF,TF,false,CSAs...>::scale( const Other& scalar )
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
/*!\brief Returns whether the sparse subvector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse subvector, \a false if not.
//
// This function returns whether the given address can alias with the sparse subvector. In
// contrast to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,AF,TF,false,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the sparse subvector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this sparse subvector, \a false if not.
//
// This function returns whether the given address is aliased with the sparse subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,AF,TF,false,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the subvector can be used in SMP assignments.
//
// \return \a true in case the subvector can be used in SMP assignments, \a false if not.
//
// This function returns whether the subvector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// vector).
*/
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool Subvector<VT,AF,TF,false,CSAs...>::canSMPAssign() const noexcept
{
   return false;
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Subvector<VT,AF,TF,false,CSAs...>::assign( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   reserve( (*rhs).size() );

   for( size_t i=0UL; i<size(); ++i ) {
      append( i, (*rhs)[i], true );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,AF,TF,false,CSAs...>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   reserve( (*rhs).nonZeros() );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      append( element->index(), element->value(), true );
   }
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Subvector<VT,AF,TF,false,CSAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   reset();
   assign( tmp );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,AF,TF,false,CSAs...>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   using AddType = AddTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   reset();
   assign( tmp );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline void Subvector<VT,AF,TF,false,CSAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   reset();
   assign( tmp );
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
template< typename VT       // Type of the sparse vector
        , AlignmentFlag AF  // Alignment flag
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,AF,TF,false,CSAs...>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   using SubType = SubTrait_t< ResultType, ResultType_t<VT2> >;

   BLAZE_CONSTRAINT_MUST_BE_SPARSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   reset();
   assign( tmp );
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
