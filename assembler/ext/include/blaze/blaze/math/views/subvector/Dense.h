//=================================================================================================
/*!
//  \file blaze/math/views/subvector/Dense.h
//  \brief Subvector specialization for dense vectors
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

#ifndef _BLAZE_MATH_VIEWS_SUBVECTOR_DENSE_H_
#define _BLAZE_MATH_VIEWS_SUBVECTOR_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/Subvector.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/dense/InitializerVector.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/Computation.h>
#include <blaze/math/expressions/CrossExpr.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/Forward.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDDiv.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/subvector/BaseTemplate.h>
#include <blaze/math/views/subvector/SubvectorData.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/HostDevice.h>
#include <blaze/system/Inline.h>
#include <blaze/system/MacroDisable.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR UNALIGNED DENSE SUBVECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Subvector for unaligned dense subvectors.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the requirements of unaligned
// dense subvectors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector<VT,unaligned,TF,true,CSAs...>
   : public View< DenseVector< Subvector<VT,unaligned,TF,true,CSAs...>, TF > >
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
   using This = Subvector<VT,unaligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = VT;                            //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<VT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;             //!< Type of the subvector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the subvector elements.
   using ReturnType    = ReturnType_t<VT>;              //!< Return type for expression template evaluations
   using CompositeType = const Subvector&;              //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   using ConstReference = ConstReference_t<VT>;

   //! Reference to a non-constant subvector value.
   using Reference = If_t< IsConst_v<VT>, ConstReference, Reference_t<VT> >;

   //! Pointer to a constant subvector value.
   using ConstPointer = ConstPointer_t<VT>;

   //! Pointer to a non-constant subvector value.
   using Pointer = If_t< IsConst_v<VT> || !HasMutableDataAccess_v<VT>, ConstPointer, Pointer_t<VT> >;
   //**********************************************************************************************

   //**SubvectorIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the dense subvector.
   */
   template< typename IteratorType >  // Type of the dense vector iterator
   class SubvectorIterator
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
      /*!\brief Default constructor of the SubvectorIterator class.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator()
         : iterator_ (       )  // Iterator to the current subvector element
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubvectorIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator( IteratorType iterator, bool isMemoryAligned )
         : iterator_ ( iterator        )  // Iterator to the current subvector element
         , isAligned_( isMemoryAligned )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubvectorIterator instances.
      //
      // \param it The subvector iterator to be copied
      */
      template< typename IteratorType2 >
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator( const SubvectorIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current subvector element
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator& operator+=( size_t inc ) {
         iterator_ += inc;
         return *this;
      }
      //*******************************************************************************************

      //**Subtraction assignment operator**********************************************************
      /*!\brief Subtraction assignment operator.
      //
      // \param dec The decrement of the iterator.
      // \return The decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const SubvectorIterator operator++( int ) {
         return SubvectorIterator( iterator_++, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline BLAZE_DEVICE_CALLABLE SubvectorIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline BLAZE_DEVICE_CALLABLE const SubvectorIterator operator--( int ) {
         return SubvectorIterator( iterator_--, isAligned_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return Reference to the current value.
      */
      inline BLAZE_DEVICE_CALLABLE ReferenceType operator*() const {
         return *iterator_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return Pointer to the element at the current iterator position.
      */
      inline BLAZE_DEVICE_CALLABLE IteratorType operator->() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Load of a SIMD element of the dense subvector.
      //
      // \return The loaded SIMD element.
      //
      // This function performs a load of the current SIMD element of the subvector iterator.
      // This function must \b NOT be called explicitly! It is used internally for the performance
      // optimized evaluation of expression templates. Calling this function explicitly might
      // result in erroneous results and/or in compilation errors.
      */
      inline SIMDType load() const {
         if( isAligned_ ) {
            return loada();
         }
         else {
            return loadu();
         }
      }
      //*******************************************************************************************

      //**Loada function***************************************************************************
      /*!\brief Aligned load of a SIMD element of the dense subvector.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an aligned load of the current SIMD element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loada() const {
         return iterator_.loada();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of a SIMD element of the dense subvector.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an unaligned load of the current SIMD element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loadu() const {
         return iterator_.loadu();
      }
      //*******************************************************************************************

      //**Store function***************************************************************************
      /*!\brief Store of a SIMD element of the dense subvector.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs a store of the current SIMD element of the subvector iterator.
      // This function must \b NOT be called explicitly! It is used internally for the performance
      // optimized evaluation of expression templates. Calling this function explicitly might
      // result in erroneous results and/or in compilation errors.
      */
      inline void store( const SIMDType& value ) const {
         if( isAligned_ ) {
            storea( value );
         }
         else {
            storeu( value );
         }
      }
      //*******************************************************************************************

      //**Storea function**************************************************************************
      /*!\brief Aligned store of a SIMD element of the dense subvector.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned store of the current SIMD element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storea( const SIMDType& value ) const {
         iterator_.storea( value );
      }
      //*******************************************************************************************

      //**Storeu function**************************************************************************
      /*!\brief Unaligned store of a SIMD element of the dense subvector.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an unaligned store of the current SIMD element of the subvector
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storeu( const SIMDType& value ) const {
         iterator_.storeu( value );
      }
      //*******************************************************************************************

      //**Stream function**************************************************************************
      /*!\brief Aligned, non-temporal store of a SIMD element of the dense subvector.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned, non-temporal store of the current SIMD element of the
      // subvector iterator. This function must \b NOT be called explicitly! It is used internally
      // for the performance optimized evaluation of expression templates. Calling this function
      // explicitly might result in erroneous results and/or in compilation errors.
      */
      inline void stream( const SIMDType& value ) const {
         iterator_.stream( value );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator==( const SubvectorIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator!=( const SubvectorIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<( const SubvectorIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>( const SubvectorIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator<=( const SubvectorIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubvectorIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline BLAZE_DEVICE_CALLABLE bool operator>=( const SubvectorIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline BLAZE_DEVICE_CALLABLE DifferenceType operator-( const SubvectorIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubvectorIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const SubvectorIterator operator+( const SubvectorIterator& it, size_t inc ) {
         return SubvectorIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubvectorIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const SubvectorIterator operator+( size_t inc, const SubvectorIterator& it ) {
         return SubvectorIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubvectorIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline BLAZE_DEVICE_CALLABLE const SubvectorIterator operator-( const SubvectorIterator& it, size_t dec ) {
         return SubvectorIterator( it.iterator_ - dec, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the subvector iterator.
      //
      // \return The current position of the subvector iterator.
      */
      inline BLAZE_DEVICE_CALLABLE IteratorType base() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**IsAligned function***********************************************************************
      /*!\brief Access to the iterator's memory alignment flag.
      //
      // \return \a true in case the iterator is aligned, \a false if it is not.
      */
      inline bool isAligned() const {
         return isAligned_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;   //!< Iterator to the current subvector element.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = SubvectorIterator< ConstIterator_t<VT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<VT>, ConstIterator, SubvectorIterator< Iterator_t<VT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = VT::simdEnabled;

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
                            inline Subvector& operator= ( const ElementType& rhs );
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

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   template< typename Other > inline Subvector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDDiv_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename VT2, AlignmentFlag AF2, bool TF2, size_t... CSAs2 >
   inline bool canAlias( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT2, AlignmentFlag AF2, bool TF2, size_t... CSAs2 >
   inline bool isAliased( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT2 >
   inline auto assign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT2> >;

   template< typename VT2 >
   inline auto assign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT2> >;

   template< typename VT2 > inline void assign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto addAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT2> >;

   template< typename VT2 >
   inline auto addAssign ( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT2> >;

   template< typename VT2 > inline void addAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto subAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT2> >;

   template< typename VT2 >
   inline auto subAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT2> >;

   template< typename VT2 > inline void subAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto multAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT2> >;

   template< typename VT2 >
   inline auto multAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT2> >;

   template< typename VT2 > inline void multAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto divAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT2> >;

   template< typename VT2 >
   inline auto divAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT2> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand vector_;        //!< The vector containing the subvector.
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the subvector is fully aligned
                                with respect to the given element type and the available instruction
                                set. In case the subvector is fully aligned it is possible to use
                                aligned loads and stores instead of unaligned loads and stores. In
                                order to be aligned, the first element of the subvector must be
                                aligned. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename VT2, AlignmentFlag AF2, bool TF2, bool DF2, size_t... CSAs2 > friend class Subvector;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
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
/*!\brief Constructor for unaligned dense subvectors.
//
// \param vector The dense vector containing the subvector.
// \param args The runtime subvector arguments.
// \exception std::invalid_argument Invalid subvector specification.
//
// By default, the provided subvector arguments are checked at runtime. In case the subvector is
// not properly specified (i.e. if the specified offset is greater than the size of the given
// vector or the subvector is specified beyond the size of the vector) a \a std::invalid_argument
// exception is thrown. The checks can be skipped by providing the optional \a blaze::unchecked
// argument.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , size_t... CSAs >    // Compile time subvector arguments
template< typename... RSAs >  // Runtime subvector arguments
inline Subvector<VT,unaligned,TF,true,CSAs...>::Subvector( VT& vector, RSAs... args )
   : DataType  ( args... )  // Base class initialization
   , vector_   ( vector  )  // The vector containing the subvector
   , isAligned_( simdEnabled && IsContiguous_v<VT> &&
                 vector.data() != nullptr && checkAlignment( data() ) )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::Reference
   Subvector<VT,unaligned,TF,true,CSAs...>::operator[]( size_t index )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstReference
   Subvector<VT,unaligned,TF,true,CSAs...>::operator[]( size_t index ) const
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
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::Reference
   Subvector<VT,unaligned,TF,true,CSAs...>::at( size_t index )
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
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid subvector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access index.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstReference
   Subvector<VT,unaligned,TF,true,CSAs...>::at( size_t index ) const
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
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::Pointer
   Subvector<VT,unaligned,TF,true,CSAs...>::data() noexcept
{
   return vector_.data() + offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstPointer
   Subvector<VT,unaligned,TF,true,CSAs...>::data() const noexcept
{
   return vector_.data() + offset();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::Iterator
   Subvector<VT,unaligned,TF,true,CSAs...>::begin()
{
   return Iterator( vector_.begin() + offset(), isAligned_ );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,unaligned,TF,true,CSAs...>::begin() const
{
   return ConstIterator( vector_.cbegin() + offset(), isAligned_ );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,unaligned,TF,true,CSAs...>::cbegin() const
{
   return ConstIterator( vector_.cbegin() + offset(), isAligned_ );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::Iterator
   Subvector<VT,unaligned,TF,true,CSAs...>::end()
{
   return Iterator( vector_.begin() + offset() + size(), isAligned_ );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,unaligned,TF,true,CSAs...>::end() const
{
   return ConstIterator( vector_.cbegin() + offset() + size(), isAligned_ );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,unaligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,unaligned,TF,true,CSAs...>::cend() const
{
   return ConstIterator( vector_.cbegin() + offset() + size(), isAligned_ );
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
/*!\brief Homogenous assignment to all subvector elements.
//
// \param rhs Scalar value to be assigned to all subvector elements.
// \return Reference to the assigned subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t iend( offset() + size() );
   decltype(auto) left( derestrict( vector_ ) );

   for( size_t i=offset(); i<iend; ++i ) {
      if( !IsRestricted_v<VT> || trySet( vector_, i, rhs ) )
         left[i] = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all subvector elements.
//
// \param list The initializer list.
// \return Reference to the assigned subvector.
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to subvector" );
   }

   if( IsRestricted_v<VT> ) {
      const InitializerVector<ElementType,TF> tmp( list, size() );
      if( !tryAssign( vector_, tmp, offset() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *this ) );

   std::fill( std::copy( list.begin(), list.end(), left.begin() ), left.end(), ElementType() );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Subvector.
//
// \param rhs Dense subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator=( const Subvector& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && offset() == rhs.offset() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Subvector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator=( const Vector<VT2,TF>& rhs )
{
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
// \param rhs The right-hand side vector to be added to the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator+=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryAddAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \param rhs The right-hand side vector to be subtracted from the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator-=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !trySubAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \param rhs The right-hand side vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator*=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryMultAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator/=( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryDivAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::operator%=( const Vector<VT2,TF>& rhs )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline VT& Subvector<VT,unaligned,TF,true,CSAs...>::operand() noexcept
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline const VT& Subvector<VT,unaligned,TF,true,CSAs...>::operand() const noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the dense subvector.
//
// \return The minimum capacity of the dense subvector.
//
// This function returns the minimum capacity of the dense subvector, which corresponds to the
// current size plus padding.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,unaligned,TF,true,CSAs...>::spacing() const noexcept
{
   return vector_.spacing() - offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense subvector.
//
// \return The maximum capacity of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,unaligned,TF,true,CSAs...>::capacity() const noexcept
{
   return vector_.capacity() - offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,unaligned,TF,true,CSAs...>::nonZeros() const
{
   size_t nonzeros( 0 );

   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i ) {
      if( !isDefault( vector_[i] ) )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline void Subvector<VT,unaligned,TF,true,CSAs...>::reset()
{
   using blaze::clear;

   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i )
      clear( vector_[i] );
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
/*!\brief Scaling of the dense subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the dense subvector.
//
// This function scales the subvector by applying the given scalar value \a scalar to each
// element of the subvector. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the scalar value
inline Subvector<VT,unaligned,TF,true,CSAs...>&
   Subvector<VT,unaligned,TF,true,CSAs...>::scale( const Other& scalar )
{
   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i )
      vector_[i] *= scalar;
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
/*!\brief Returns whether the dense subvector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,unaligned,TF,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector can alias with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT        // Type of the dense vector
        , bool TF            // Transpose flag
        , size_t... CSAs >   // Compile time subvector arguments
template< typename VT2       // Data type of the foreign dense subvector
        , AlignmentFlag AF2  // Alignment flag of the foreign dense subvector
        , bool TF2           // Transpose flag of the foreign dense subvector
        , size_t... CSAs2 >  // Compile time subvector arguments of the foreign dense subvector
inline bool
   Subvector<VT,unaligned,TF,true,CSAs...>::canAlias( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset() + size() > alias->offset() ) &&
            ( offset() < alias->offset() + alias->size() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,unaligned,TF,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector is aliased with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT        // Type of the dense vector
        , bool TF            // Transpose flag
        , size_t... CSAs >   // Compile time subvector arguments
template< typename VT2       // Data type of the foreign dense subvector
        , AlignmentFlag AF2  // Alignment flag of the foreign dense subvector
        , bool TF2           // Transpose flag of the foreign dense subvector
        , size_t... CSAs2 >  // Compile time subvector arguments of the foreign dense subvector
inline bool
   Subvector<VT,unaligned,TF,true,CSAs...>::isAliased( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset() + size() > alias->offset() ) &&
            ( offset() < alias->offset() + alias->size() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the subvector is properly aligned in memory.
//
// \return \a true in case the subvector is aligned, \a false if not.
//
// This function returns whether the subvector is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the subvector are guaranteed to conform to the
// alignment restrictions of the underlying element type.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool Subvector<VT,unaligned,TF,true,CSAs...>::isAligned() const noexcept
{
   return isAligned_;
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
// subvector).
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool Subvector<VT,unaligned,TF,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,unaligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,unaligned,TF,true,CSAs...>::load( size_t index ) const noexcept
{
   if( isAligned_ )
      return loada( index );
   else
      return loadu( index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,unaligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,unaligned,TF,true,CSAs...>::loada( size_t index ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   return vector_.loada( offset()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense
// subvector. The index must be smaller than the number of subvector elements and it must be
// a multiple of the number of values inside the SIMD element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,unaligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,unaligned,TF,true,CSAs...>::loadu( size_t index ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );

   return vector_.loadu( offset()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,unaligned,TF,true,CSAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   if( isAligned_ )
      storea( index, value );
   else
      storeu( index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,unaligned,TF,true,CSAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   vector_.storea( offset()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,unaligned,TF,true,CSAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );

   vector_.storeu( offset()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the
// dense subvector. The index must be smaller than the number of subvector elements and it
// must be a multiple of the number of values inside the SIMD element. This function
// must \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result in
// erroneous results and/or in compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,unaligned,TF,true,CSAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   if( isAligned_ )
      vector_.stream( offset()+index, value );
   else
      vector_.storeu( offset()+index, value );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::assign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] = (*rhs)[i    ];
      vector_[offset()+i+1UL] = (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] = (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::assign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   if( useStreaming && isAligned_ &&
       ( size() > ( cacheSize/( sizeof(ElementType) * 3UL ) ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<size(); ++i ) {
         *left = *right;
      }
   }
   else
   {
      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<size(); ++i ) {
         *left = *right; ++left; ++right;
      }
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,unaligned,TF,true,CSAs...>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] = element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] += (*rhs)[i    ];
      vector_[offset()+i+1UL] += (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] += (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left += *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,unaligned,TF,true,CSAs...>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] += element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] -= (*rhs)[i    ];
      vector_[offset()+i+1UL] -= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] -= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left -= *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,unaligned,TF,true,CSAs...>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] -= element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::multAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] *= (*rhs)[i    ];
      vector_[offset()+i+1UL] *= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] *= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::multAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left *= *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,unaligned,TF,true,CSAs...>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( vector_[offset()+i] );
      vector_[offset()+i] *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( vector_[offset()+i] );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::divAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] /= (*rhs)[i    ];
      vector_[offset()+i+1UL] /= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] /= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,unaligned,TF,true,CSAs...>::divAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left /= *right; ++left; ++right;
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED DENSE SUBVECTORS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Subvector for aligned dense subvectors.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the requirements of aligned
// dense subvectors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector<VT,aligned,TF,true,CSAs...>
   : public View< DenseVector< Subvector<VT,aligned,TF,true,CSAs...>, TF > >
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
   using This = Subvector<VT,aligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = VT;                            //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<VT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<VT>;             //!< Type of the subvector elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the subvector elements.
   using ReturnType    = ReturnType_t<VT>;              //!< Return type for expression template evaluations
   using CompositeType = const Subvector&;              //!< Data type for composite expression templates.

   //! Reference to a constant subvector value.
   using ConstReference = ConstReference_t<VT>;

   //! Reference to a non-constant subvector value.
   using Reference = If_t< IsConst_v<VT>, ConstReference, Reference_t<VT> >;

   //! Pointer to a constant subvector value.
   using ConstPointer = ConstPointer_t<VT>;

   //! Pointer to a non-constant subvector value.
   using Pointer = If_t< IsConst_v<VT> || !HasMutableDataAccess_v<VT>, ConstPointer, Pointer_t<VT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<VT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<VT>, ConstIterator, Iterator_t<VT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = VT::simdEnabled;

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
                            inline Subvector& operator= ( const ElementType& rhs );
                            inline Subvector& operator= ( initializer_list<ElementType> list );
                            inline Subvector& operator= ( const Subvector& rhs );
   template< typename VT2 > inline Subvector& operator= ( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator+=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator-=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator*=( const Vector<VT2,TF>& rhs );
   template< typename VT2 > inline Subvector& operator/=( const DenseVector<VT2,TF>&  rhs );
   template< typename VT2 > inline Subvector& operator%=( const Vector<VT2,TF>&  rhs );
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::offset;
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
   template< typename Other > inline Subvector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && VT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedMultAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper structure for the explicit application of the SFINAE principle.
   template< typename VT2 >
   static constexpr bool VectorizedDivAssign_v =
      ( VectorizedAssign_v<VT2> &&
        HasSIMDDiv_v< ElementType, ElementType_t<VT2> > );
   //**********************************************************************************************

   //**SIMD properties*****************************************************************************
   //! The number of elements packed within a single SIMD element.
   static constexpr size_t SIMDSIZE = SIMDTrait<ElementType>::size;
   //**********************************************************************************************

 public:
   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other >
   inline bool canAlias( const Other* alias ) const noexcept;

   template< typename VT2, AlignmentFlag AF2, bool TF2, size_t... CSAs2 >
   inline bool canAlias( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename VT2, AlignmentFlag AF2, bool TF2, size_t... CSAs2 >
   inline bool isAliased( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t index ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t index ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t index, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t index, const SIMDType& value ) noexcept;

   template< typename VT2 >
   inline auto assign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedAssign_v<VT2> >;

   template< typename VT2 >
   inline auto assign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedAssign_v<VT2> >;

   template< typename VT2 > inline void assign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto addAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<VT2> >;

   template< typename VT2 >
   inline auto addAssign ( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<VT2> >;

   template< typename VT2 > inline void addAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto subAssign ( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<VT2> >;

   template< typename VT2 >
   inline auto subAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<VT2> >;

   template< typename VT2 > inline void subAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto multAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedMultAssign_v<VT2> >;

   template< typename VT2 >
   inline auto multAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedMultAssign_v<VT2> >;

   template< typename VT2 > inline void multAssign( const SparseVector<VT2,TF>& rhs );

   template< typename VT2 >
   inline auto divAssign( const DenseVector <VT2,TF>& rhs ) -> DisableIf_t< VectorizedDivAssign_v<VT2> >;

   template< typename VT2 >
   inline auto divAssign( const DenseVector <VT2,TF>& rhs ) -> EnableIf_t< VectorizedDivAssign_v<VT2> >;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand vector_;  //!< The vector containing the subvector.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename VT2, AlignmentFlag AF2, bool TF2, bool DF2, size_t... CSAs2 > friend class Subvector;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE   ( VT );
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
/*!\brief Constructor for aligned dense subvectors.
//
// \param vector The dense vector containing the subvector.
// \param args The runtime subvector arguments.
// \exception std::invalid_argument Invalid subvector specification.
//
// By default, the provided subvector arguments are checked at runtime. In case the subvector is
// not properly specified (i.e. if the specified offset is greater than the size of the given
// vector or the subvector is specified beyond the size of the vector) a \a std::invalid_argument
// exception is thrown. The checks can be skipped by providing the optional \a blaze::unchecked
// argument.
*/
template< typename VT         // Type of the dense vector
        , bool TF             // Transpose flag
        , size_t... CSAs >    // Compile time subvector arguments
template< typename... RSAs >  // Runtime subvector arguments
inline Subvector<VT,aligned,TF,true,CSAs...>::Subvector( VT& vector, RSAs... args )
   : DataType( args... )  // Base class initialization
   , vector_ ( vector  )  // The vector containing the subvector
{
   if( isChecked( args... ) )
   {
      if( offset() + size() > vector.size() ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector specification" );
      }

      if( simdEnabled && IsContiguous_v<VT> && !checkAlignment( data() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid subvector alignment" );
      }
   }
   else
   {
      BLAZE_USER_ASSERT( offset() + size() <= vector.size(), "Invalid subvector specification" );
      BLAZE_USER_ASSERT( !simdEnabled || !IsContiguous_v<VT> || checkAlignment( data() ), "Invalid subvector alignment" );
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
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::Reference
   Subvector<VT,aligned,TF,true,CSAs...>::operator[]( size_t index )
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
// \param index Access index. The index must be smaller than the number of subvector columns.
// \return Reference to the accessed value.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstReference
   Subvector<VT,aligned,TF,true,CSAs...>::operator[]( size_t index ) const
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::Reference
   Subvector<VT,aligned,TF,true,CSAs...>::at( size_t index )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstReference
   Subvector<VT,aligned,TF,true,CSAs...>::at( size_t index ) const
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
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::Pointer
   Subvector<VT,aligned,TF,true,CSAs...>::data() noexcept
{
   return vector_.data() + offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the subvector elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstPointer
   Subvector<VT,aligned,TF,true,CSAs...>::data() const noexcept
{
   return vector_.data() + offset();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::Iterator
   Subvector<VT,aligned,TF,true,CSAs...>::begin()
{
   return ( vector_.begin() + offset() );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,aligned,TF,true,CSAs...>::begin() const
{
   return ( vector_.cbegin() + offset() );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,aligned,TF,true,CSAs...>::cbegin() const
{
   return ( vector_.cbegin() + offset() );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::Iterator
   Subvector<VT,aligned,TF,true,CSAs...>::end()
{
   return ( vector_.begin() + offset() + size() );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,aligned,TF,true,CSAs...>::end() const
{
   return ( vector_.cbegin() + offset() + size() );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline typename Subvector<VT,aligned,TF,true,CSAs...>::ConstIterator
   Subvector<VT,aligned,TF,true,CSAs...>::cend() const
{
   return ( vector_.cbegin() + offset() + size() );
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
/*!\brief Homogenous assignment to all subvector elements.
//
// \param rhs Scalar value to be assigned to all subvector elements.
// \return Reference to the assigned subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t iend( offset() + size() );
   decltype(auto) left( derestrict( vector_ ) );

   for( size_t i=offset(); i<iend; ++i ) {
      if( !IsRestricted_v<VT> || trySet( vector_, i, rhs ) )
         left[i] = rhs;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all subvector elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to subvector.
//
// This assignment operator offers the option to directly assign to all elements of the subvector
// by means of an initializer list. The subvector elements are assigned the values from the given
// initializer list. Missing values are reset to their default state. Note that in case the size
// of the initializer list exceeds the size of the subvector, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator=( initializer_list<ElementType> list )
{
   if( list.size() > size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to subvector" );
   }

   if( IsRestricted_v<VT> ) {
      const InitializerVector<ElementType,TF> tmp( list, size() );
      if( !tryAssign( vector_, tmp, offset() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
      }
   }

   decltype(auto) left( derestrict( *this ) );

   std::fill( std::copy( list.begin(), list.end(), begin() ), end(), ElementType() );

   BLAZE_INTERNAL_ASSERT( isIntact( vector_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Subvector.
//
// \param rhs Dense subvector to be copied.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Subvector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two subvectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator=( const Subvector& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( &rhs == this || ( &vector_ == &rhs.vector_ && offset() == rhs.offset() ) )
      return *this;

   if( size() != rhs.size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Subvector sizes do not match" );
   }

   if( !tryAssign( vector_, rhs, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( *rhs );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator=( const Vector<VT2,TF>& rhs )
{
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
// \param rhs The right-hand side vector to be added to the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator+=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryAddAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \param rhs The right-hand side vector to be subtracted from the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator-=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !trySubAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \param rhs The right-hand side vector to be multiplied with the dense subvector.
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator*=( const Vector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryMultAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Vector sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator/=( const DenseVector<VT2,TF>& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT2>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT2> );

   if( size() != (*rhs).size() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<VT>, CompositeType_t<VT2>, const VT2& >;
   Right right( *rhs );

   if( !tryDivAssign( vector_, right, offset() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted vector" );
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
// \return Reference to the assigned subvector.
// \exception std::invalid_argument Invalid vector size for cross product.
// \exception std::invalid_argument Invalid assignment to restricted vector.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side vector
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::operator%=( const Vector<VT2,TF>& rhs )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline VT& Subvector<VT,aligned,TF,true,CSAs...>::operand() noexcept
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline const VT& Subvector<VT,aligned,TF,true,CSAs...>::operand() const noexcept
{
   return vector_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the minimum capacity of the dense subvector.
//
// \return The minimum capacity of the dense subvector.
//
// This function returns the minimum capacity of the dense subvector, which corresponds to the
// current size plus padding.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,aligned,TF,true,CSAs...>::spacing() const noexcept
{
   return vector_.spacing() - offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense subvector.
//
// \return The maximum capacity of the dense subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,aligned,TF,true,CSAs...>::capacity() const noexcept
{
   return vector_.capacity() - offset();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the subvector.
//
// \return The number of non-zero elements in the subvector.
//
// Note that the number of non-zero elements is always less than or equal to the current size
// of the subvector.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline size_t Subvector<VT,aligned,TF,true,CSAs...>::nonZeros() const
{
   size_t nonzeros( 0 );

   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i ) {
      if( !isDefault( vector_[i] ) )
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline void Subvector<VT,aligned,TF,true,CSAs...>::reset()
{
   using blaze::clear;

   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i )
      clear( vector_[i] );
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
/*!\brief Scaling of the dense subvector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the subvector scaling.
// \return Reference to the dense subvector.
//
// This function scales the subvector by applying the given scalar value \a scalar to each
// element of the subvector. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the scalar value
inline Subvector<VT,aligned,TF,true,CSAs...>&
   Subvector<VT,aligned,TF,true,CSAs...>::scale( const Other& scalar )
{
   const size_t iend( offset() + size() );
   for( size_t i=offset(); i<iend; ++i )
      vector_[i] *= scalar;
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
/*!\brief Returns whether the dense subvector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,aligned,TF,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector can alias with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address can alias with the dense subvector.
// In contrast to the isAliased() function this function is allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT        // Type of the dense vector
        , bool TF            // Transpose flag
        , size_t... CSAs >   // Compile time subvector arguments
template< typename VT2       // Data type of the foreign dense subvector
        , AlignmentFlag AF2  // Alignment flag of the foreign dense subvector
        , bool TF2           // Transpose flag of the foreign dense subvector
        , size_t... CSAs2 >  // Compile time subvector arguments of the foreign dense subvector
inline bool
   Subvector<VT,aligned,TF,true,CSAs...>::canAlias( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset() + size() > alias->offset() ) &&
            ( offset() < alias->offset() + alias->size() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename Other >  // Data type of the foreign expression
inline bool Subvector<VT,aligned,TF,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return vector_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the dense subvector is aliased with the given dense subvector \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this dense subvector, \a false if not.
//
// This function returns whether the given address is aliased with the dense subvector.
// In contrast to the canAlias() function this function is not allowed to use compile time
// expressions to optimize the evaluation.
*/
template< typename VT        // Type of the dense vector
        , bool TF            // Transpose flag
        , size_t... CSAs >   // Compile time subvector arguments
template< typename VT2       // Data type of the foreign dense subvector
        , AlignmentFlag AF2  // Alignment flag of the foreign dense subvector
        , bool TF2           // Transpose flag of the foreign dense subvector
        , size_t... CSAs2 >  // Compile time subvector arguments of the foreign dense subvector
inline bool
   Subvector<VT,aligned,TF,true,CSAs...>::isAliased( const Subvector<VT2,AF2,TF2,true,CSAs2...>* alias ) const noexcept
{
   return ( vector_.isAliased( &alias->vector_ ) &&
            ( offset() + size() > alias->offset() ) &&
            ( offset() < alias->offset() + alias->size() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the subvector is properly aligned in memory.
//
// \return \a true in case the subvector is aligned, \a false if not.
//
// This function returns whether the subvector is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of the subvector are guaranteed to conform to the
// alignment restrictions of the underlying element type.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool Subvector<VT,aligned,TF,true,CSAs...>::isAligned() const noexcept
{
   return true;
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
// subvector).
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
inline bool Subvector<VT,aligned,TF,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( size() > SMP_DVECASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense subvector. The index
// must be smaller than the number of subvector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,aligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,aligned,TF,true,CSAs...>::load( size_t index ) const noexcept
{
   return loada( index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,aligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,aligned,TF,true,CSAs...>::loada( size_t index ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   return vector_.loada( offset()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the dense subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE typename Subvector<VT,aligned,TF,true,CSAs...>::SIMDType
   Subvector<VT,aligned,TF,true,CSAs...>::loadu( size_t index ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );

   return vector_.loadu( offset()+index );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store a specific SIMD element of the dense subvector. The index
// must be smaller than the number of subvector elements and it must be a multiple of the
// number of values inside the SIMD element. This function must \b NOT be called explicitly!
// It is used internally for the performance optimized evaluation of expression templates.
// Calling this function explicitly might result in erroneous results and/or in compilation
// errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,aligned,TF,true,CSAs...>::store( size_t index, const SIMDType& value ) noexcept
{
   storea( index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,aligned,TF,true,CSAs...>::storea( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   vector_.storea( offset()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store a specific SIMD element of the dense subvector.
// The index must be smaller than the number of subvector elements and it must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,aligned,TF,true,CSAs...>::storeu( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );

   vector_.storeu( offset()+index, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the subvector.
//
// \param index Access index. The index must be smaller than the number of subvector elements.
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store a specific SIMD element of the
// dense subvector. The index must be smaller than the number of subvector elements and it
// must be a multiple of the number of values inside the SIMD element. This function
// must \b NOT be called explicitly! It is used internally for the performance optimized
// evaluation of expression templates. Calling this function explicitly might result in
// erroneous results and/or in compilation errors.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
BLAZE_ALWAYS_INLINE void
   Subvector<VT,aligned,TF,true,CSAs...>::stream( size_t index, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( index < size()            , "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index + SIMDSIZE <= size(), "Invalid subvector access index" );
   BLAZE_INTERNAL_ASSERT( index % SIMDSIZE == 0UL   , "Invalid subvector access index" );

   vector_.stream( offset()+index, value );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >   // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::assign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] = (*rhs)[i    ];
      vector_[offset()+i+1UL] = (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] = (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::assign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   if( useStreaming && size() > ( cacheSize/( sizeof(ElementType) * 3UL ) ) && !(*rhs).isAliased( this ) )
   {
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<size(); ++i ) {
         *left = *right; ++left; ++right;
      }
   }
   else
   {
      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<size(); ++i ) {
         *left = *right; ++left; ++right;
      }
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,aligned,TF,true,CSAs...>::assign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] = element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] += (*rhs)[i    ];
      vector_[offset()+i+1UL] += (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] += (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::addAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left += *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,aligned,TF,true,CSAs...>::addAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] += element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] -= (*rhs)[i    ];
      vector_[offset()+i+1UL] -= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] -= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::subAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left -= *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,aligned,TF,true,CSAs...>::subAssign( const SparseVector<VT2,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      vector_[offset()+element->index()] -= element->value();
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::multAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedMultAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] *= (*rhs)[i    ];
      vector_[offset()+i+1UL] *= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] *= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the multiplication assignment of a dense vector.
//
// \param rhs The right-hand side dense vector to be multiplied.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::multAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedMultAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left *= *right; ++left; ++right;
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side sparse vector
inline void Subvector<VT,aligned,TF,true,CSAs...>::multAssign( const SparseVector<VT2,TF>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   size_t i( 0UL );

   for( ConstIterator_t<VT2> element=(*rhs).begin(); element!=(*rhs).end(); ++element ) {
      const size_t index( element->index() );
      for( ; i<index; ++i )
         reset( vector_[offset()+i] );
      vector_[offset()+i] *= element->value();
      ++i;
   }

   for( ; i<size(); ++i ) {
      reset( vector_[offset()+i] );
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
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::divAssign( const DenseVector<VT2,TF>& rhs )
   -> DisableIf_t< VectorizedDivAssign_v<VT2> >
{
   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   for( size_t i=0UL; i<ipos; i+=2UL ) {
      vector_[offset()+i    ] /= (*rhs)[i    ];
      vector_[offset()+i+1UL] /= (*rhs)[i+1UL];
   }
   if( ipos < size() ) {
      vector_[offset()+ipos] /= (*rhs)[ipos];
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the division assignment of a dense vector.
//
// \param rhs The right-hand side dense vector divisor.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename VT       // Type of the dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
template< typename VT2 >    // Type of the right-hand side dense vector
inline auto Subvector<VT,aligned,TF,true,CSAs...>::divAssign( const DenseVector<VT2,TF>& rhs )
   -> EnableIf_t< VectorizedDivAssign_v<VT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( size() == (*rhs).size(), "Invalid vector sizes" );

   const size_t ipos( prevMultiple( size(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= size(), "Invalid end calculation" );

   size_t i( 0UL );
   Iterator left( begin() );
   ConstIterator_t<VT2> right( (*rhs).begin() );

   for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<ipos; i+=SIMDSIZE ) {
      left.store( left.load() / right.load() ); left += SIMDSIZE; right += SIMDSIZE;
   }
   for( ; i<size(); ++i ) {
      *left /= *right; ++left; ++right;
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
/*!\brief Specialization of Subvector for dense vector/dense vector cross products.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the special case of dense
// vector/dense vector cross products.
*/
template< typename VT1      // Type of the left-hand side dense vector
        , typename VT2      // Type of the right-hand side dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector< DVecDVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >
   : public View< DenseVector< Subvector< DVecDVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >, TF > >
   , private SubvectorData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = DVecDVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = SubvectorData<CSAs...>;         //!< The type of the SubvectorData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Subvector instance.
   using This = Subvector<CPE,unaligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = CPE;                           //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<RT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;            //!< Type of the subvector elements.
   using ReturnType    = ReturnType_t<CPE>;             //!< Return type for expression template evaluations
   using CompositeType = const ResultType;              //!< Data type for composite expression templates.
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
   /*!\brief Constructor for dense vector/dense vector cross product subvectors.
   //
   // \param vector The dense vector/dense vector cross product expression.
   // \param args The runtime subvector arguments.
   // \exception std::invalid_argument Invalid subvector specification.
   */
   template< typename... RSAs >  // Optional subvector arguments
   explicit inline Subvector( const CPE& vector, RSAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The dense vector/dense vector cross product expression
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
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset()+index];
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
   using DataType::offset;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the subvector.
   //
   // \return The cross product expression containing the subvector.
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
/*!\brief Specialization of Subvector for dense vector/sparse vector cross products.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the special case of dense
// vector/sparse vector cross products.
*/
template< typename VT1      // Type of the left-hand side dense vector
        , typename VT2      // Type of the right-hand side sparse vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector< DVecSVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >
   : public View< DenseVector< Subvector< DVecSVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >, TF > >
   , private SubvectorData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = DVecSVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = SubvectorData<CSAs...>;         //!< The type of the SubvectorData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Subvector instance.
   using This = Subvector<CPE,unaligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = CPE;                           //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<RT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;            //!< Type of the subvector elements.
   using ReturnType    = ReturnType_t<CPE>;             //!< Return type for expression template evaluations
   using CompositeType = const ResultType;              //!< Data type for composite expression templates.
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
   /*!\brief Constructor for dense vector/sparse vector cross product subvectors.
   //
   // \param vector The dense vector/sparse vector cross product expression.
   // \param args The runtime subvector arguments.
   // \exception std::invalid_argument Invalid subvector specification.
   */
   template< typename... RSAs >  // Optional subvector arguments
   explicit inline Subvector( const CPE& vector, RSAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The dense vector/sparse vector cross product expression
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
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset()+index];
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
   using DataType::offset;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the subvector.
   //
   // \return The cross product expression containing the subvector.
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
/*!\brief Specialization of Subvector for sparse vector/dense vector cross products.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the special case of sparse
// vector/dense vector cross products.
*/
template< typename VT1      // Type of the left-hand side sparse vector
        , typename VT2      // Type of the right-hand side dense vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector< SVecDVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >
   : public View< DenseVector< Subvector< SVecDVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >, TF > >
   , private SubvectorData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = SVecDVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = SubvectorData<CSAs...>;         //!< The type of the SubvectorData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Subvector instance.
   using This = Subvector<CPE,unaligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = CPE;                           //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<RT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;            //!< Type of the subvector elements.
   using ReturnType    = ReturnType_t<CPE>;             //!< Return type for expression template evaluations
   using CompositeType = const ResultType;              //!< Data type for composite expression templates.
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
   /*!\brief Constructor for sparse vector/dense vector cross product subvectors.
   //
   // \param vector The sparse vector/dense vector cross product expression.
   // \param args The runtime subvector arguments.
   // \exception std::invalid_argument Invalid subvector specification.
   */
   template< typename... RSAs >  // Optional subvector arguments
   explicit inline Subvector( const CPE& vector, RSAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The sparse vector/dense vector cross product expression
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
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset()+index];
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
   using DataType::offset;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the subvector.
   //
   // \return The cross product expression containing the subvector.
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
/*!\brief Specialization of Subvector for sparse vector/sparse vector cross products.
// \ingroup subvector
//
// This specialization of Subvector adapts the class template to the special case of sparse
// vector/sparse vector cross products.
*/
template< typename VT1      // Type of the left-hand side sparse vector
        , typename VT2      // Type of the right-hand side sparse vector
        , bool TF           // Transpose flag
        , size_t... CSAs >  // Compile time subvector arguments
class Subvector< SVecSVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >
   : public View< DenseVector< Subvector< SVecSVecCrossExpr<VT1,VT2,TF>, unaligned, TF, true, CSAs... >, TF > >
   , private SubvectorData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using CPE      = SVecSVecCrossExpr<VT1,VT2,TF>;  //!< Type of the cross product expression.
   using RT       = ResultType_t<CPE>;              //!< Result type of the cross product expression.
   using DataType = SubvectorData<CSAs...>;         //!< The type of the SubvectorData base class.
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Subvector instance.
   using This = Subvector<CPE,unaligned,TF,true,CSAs...>;

   using BaseType      = View< DenseVector<This,TF> >;  //!< Base type of this Subvector instance.
   using ViewedType    = CPE;                           //!< The type viewed by this Subvector instance.
   using ResultType    = SubvectorTrait_t<RT,CSAs...>;  //!< Result type for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<CPE>;            //!< Type of the subvector elements.
   using ReturnType    = ReturnType_t<CPE>;             //!< Return type for expression template evaluations
   using CompositeType = const ResultType;              //!< Data type for composite expression templates.
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
   /*!\brief Constructor for sparse vector/sparse vector cross product subvectors.
   //
   // \param vector The sparse vector/sparse vector cross product expression.
   // \param args The runtime subvector arguments.
   // \exception std::invalid_argument Invalid subvector specification.
   */
   template< typename... RSAs >  // Optional subvector arguments
   explicit inline Subvector( const CPE& vector, RSAs... args )
      : DataType( args... )  // Base class initialization
      , vector_ ( vector  )  // The sparse vector/sparse vector cross product expression
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
   //**********************************************************************************************

   //**Subscript operator**************************************************************************
   /*!\brief Subscript operator for the direct access to the vector elements.
   //
   // \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
   // \return The resulting value.
   */
   inline ReturnType operator[]( size_t index ) const {
      BLAZE_INTERNAL_ASSERT( index < size(), "Invalid vector access index" );
      return vector_[offset()+index];
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
   using DataType::offset;
   using DataType::size;
   //**********************************************************************************************

   //**Operand access******************************************************************************
   /*!\brief Returns the cross product expression containing the subvector.
   //
   // \return The cross product expression containing the subvector.
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
   inline bool isAliased( const T* alias ) const {
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
