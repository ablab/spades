//=================================================================================================
/*!
//  \file blaze/math/views/submatrix/Dense.h
//  \brief Submatrix specialization for dense matrices
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

#ifndef _BLAZE_MATH_VIEWS_SUBMATRIX_DENSE_H_
#define _BLAZE_MATH_VIEWS_SUBMATRIX_DENSE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <iterator>
#include <blaze/math/Aliases.h>
#include <blaze/math/AlignmentFlag.h>
#include <blaze/math/constraints/ColumnMajorMatrix.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/DenseMatrix.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/RowMajorMatrix.h>
#include <blaze/math/constraints/Submatrix.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/TransExpr.h>
#include <blaze/math/constraints/UniTriangular.h>
#include <blaze/math/dense/InitializerMatrix.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseMatrix.h>
#include <blaze/math/expressions/View.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/PrevMultiple.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/math/SIMD.h>
#include <blaze/math/StorageOrder.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/SchurTrait.h>
#include <blaze/math/traits/SubmatrixTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/typetraits/HasMutableDataAccess.h>
#include <blaze/math/typetraits/HasSIMDAdd.h>
#include <blaze/math/typetraits/HasSIMDMult.h>
#include <blaze/math/typetraits/HasSIMDSub.h>
#include <blaze/math/typetraits/IsContiguous.h>
#include <blaze/math/typetraits/IsDiagonal.h>
#include <blaze/math/typetraits/IsExpression.h>
#include <blaze/math/typetraits/IsHermitian.h>
#include <blaze/math/typetraits/IsLower.h>
#include <blaze/math/typetraits/IsRestricted.h>
#include <blaze/math/typetraits/IsSIMDCombinable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsStrictlyLower.h>
#include <blaze/math/typetraits/IsStrictlyUpper.h>
#include <blaze/math/typetraits/IsSymmetric.h>
#include <blaze/math/typetraits/IsTriangular.h>
#include <blaze/math/typetraits/IsUniLower.h>
#include <blaze/math/typetraits/IsUniUpper.h>
#include <blaze/math/typetraits/IsUpper.h>
#include <blaze/math/typetraits/RequiresEvaluation.h>
#include <blaze/math/views/Check.h>
#include <blaze/math/views/submatrix/BaseTemplate.h>
#include <blaze/math/views/submatrix/SubmatrixData.h>
#include <blaze/system/Blocking.h>
#include <blaze/system/CacheSize.h>
#include <blaze/system/Inline.h>
#include <blaze/system/Optimizations.h>
#include <blaze/system/Thresholds.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/AlignmentCheck.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Vectorizable.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/MaybeUnused.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsConst.h>
#include <blaze/util/typetraits/IsReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR UNALIGNED ROW-MAJOR DENSE SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Submatrix for unaligned row-major dense submatrices.
// \ingroup submatrix
//
// This Specialization of Submatrix adapts the class template to the requirements of unaligned
// row-major dense submatrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
class Submatrix<MT,unaligned,false,true,CSAs...>
   : public View< DenseMatrix< Submatrix<MT,unaligned,false,true,CSAs...>, false > >
   , private SubmatrixData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = SubmatrixData<CSAs...>;               //!< The type of the SubmatrixData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Submatrix instance.
   using This = Submatrix<MT,unaligned,false,true,CSAs...>;

   //! Base type of this Submatrix instance.
   using BaseType = View< DenseMatrix<This,false> >;

   using ViewedType    = MT;                            //!< The type viewed by this Submatrix instance.
   using ResultType    = SubmatrixTrait_t<MT,CSAs...>;  //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;    //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;             //!< Type of the submatrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the submatrix elements.
   using ReturnType    = ReturnType_t<MT>;              //!< Return type for expression template evaluations
   using CompositeType = const Submatrix&;              //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant submatrix value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant submatrix value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant submatrix value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename IteratorType >  // Type of the dense matrix iterator
   class SubmatrixIterator
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
      /*!\brief Default constructor of the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : iterator_ (       )  // Iterator to the current submatrix element
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubmatrixIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline SubmatrixIterator( IteratorType iterator, bool isMemoryAligned )
         : iterator_ ( iterator        )  // Iterator to the current submatrix element
         , isAligned_( isMemoryAligned )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current submatrix element
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline SubmatrixIterator& operator+=( size_t inc ) {
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
      inline SubmatrixIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         return SubmatrixIterator( iterator_++, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline SubmatrixIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator--( int ) {
         return SubmatrixIterator( iterator_--, isAligned_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReferenceType operator*() const {
         return *iterator_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return Pointer to the element at the current iterator position.
      */
      inline IteratorType operator->() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs a load of the current SIMD element of the submatrix iterator.
      // This function must \b NOT be called explicitly! It is used internally for the performance
      // optimized evaluation of expression templates. Calling this function explicitly might
      // result in erroneous results and/or in compilation errors.
      */
      inline SIMDType load() const noexcept {
         if( isAligned_ )
            return loada();
         else
            return loadu();
      }
      //*******************************************************************************************

      //**Loada function***************************************************************************
      /*!\brief Aligned load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an aligned load of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for
      // the performance optimized evaluation of expression templates. Calling this function
      // explicitly might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loada() const noexcept {
         return iterator_.loada();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an unaligned load of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loadu() const noexcept {
         return iterator_.loadu();
      }
      //*******************************************************************************************

      //**Store function***************************************************************************
      /*!\brief Store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs a store of the current SIMD element of the submatrix iterator.
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
      /*!\brief Aligned store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned store of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storea( const SIMDType& value ) const {
         iterator_.storea( value );
      }
      //*******************************************************************************************

      //**Storeu function**************************************************************************
      /*!\brief Unaligned store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an unaligned store of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storeu( const SIMDType& value ) const {
         iterator_.storeu( value );
      }
      //*******************************************************************************************

      //**Stream function**************************************************************************
      /*!\brief Aligned, non-temporal store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned, non-temporal store of the current SIMD element of the
      // submatrix iterator. This function must \b NOT be called explicitly! It is used internally
      // for the performance optimized evaluation of expression templates. Calling this function
      // explicitly might result in erroneous results and/or in compilation errors.
      */
      inline void stream( const SIMDType& value ) const {
         iterator_.stream( value );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const SubmatrixIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const SubmatrixIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const SubmatrixIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const SubmatrixIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const SubmatrixIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const SubmatrixIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( const SubmatrixIterator& it, size_t inc ) {
         return SubmatrixIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubmatrixIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( size_t inc, const SubmatrixIterator& it ) {
         return SubmatrixIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param dec The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const SubmatrixIterator operator-( const SubmatrixIterator& it, size_t dec ) {
         return SubmatrixIterator( it.iterator_ - dec, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**IsAligned function***********************************************************************
      /*!\brief Access to the iterator's memory alignment flag.
      //
      // \return \a true in case the iterator is aligned, \a false if it is not.
      */
      inline bool isAligned() const noexcept {
         return isAligned_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;   //!< Iterator to the current submatrix element.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = SubmatrixIterator< ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, SubmatrixIterator< Iterator_t<MT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RSAs >
   explicit inline Submatrix( MT& matrix, RSAs... args );

   Submatrix( const Submatrix& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Submatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Submatrix& operator=( const ElementType& rhs );
   inline Submatrix& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Submatrix& operator=( const Submatrix& rhs );

   template< typename MT2, bool SO2 >
   inline Submatrix& operator=( const Matrix<MT2,SO2>& rhs );

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO2 >
   inline auto operator+=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO2 >
   inline auto operator-=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO2 >
   inline auto operator%=( const Matrix<MT2,SO2>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;
   using DataType::column;
   using DataType::rows;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Submatrix& transpose();
   inline Submatrix& ctranspose();

   template< typename Other > inline Submatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<MT2> > );
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

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;        //!< The matrix containing the submatrix.
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the submatrix is fully aligned
                                with respect to the given element type and the available instruction
                                set. In case the submatrix is fully aligned it is possible to use
                                aligned loads and stores instead of unaligned loads and stores. In
                                order to be aligned, the first element of each row/column must be
                                aligned. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, AlignmentFlag AF2, bool SO2, bool DF2, size_t... CSAs2 > friend class Submatrix;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
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
/*!\brief Constructor for unaligned row-major dense submatrices.
//
// \param matrix The dense matrix containing the submatrix.
// \param args The runtime submatrix arguments.
// \exception std::invalid_argument Invalid submatrix specification.
//
// By default, the provided submatrix arguments are checked at runtime. In case the submatrix is
// not properly specified (i.e. if the specified submatrix is not contained in the given dense
// matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by providing
// the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CSAs >    // Compile time submatrix arguments
template< typename... RSAs >  // Runtime submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>::Submatrix( MT& matrix, RSAs... args )
   : DataType  ( args... )  // Base class initialization
   , matrix_   ( matrix  )  // The matrix containing the submatrix
   , isAligned_( simdEnabled && IsContiguous_v<MT> &&
                 matrix.data() != nullptr && checkAlignment( data() ) &&
                 ( rows() < 2UL || ( matrix.spacing() % SIMDSIZE ) == 0UL ) )
{
   if( isChecked( args... ) ) {
      if( ( row() + rows() > matrix_.rows() ) || ( column() + columns() > matrix_.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row()    + rows()    <= matrix_.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column() + columns() <= matrix_.columns(), "Invalid submatrix specification" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Reference
   Submatrix<MT,unaligned,false,true,CSAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstReference
   Submatrix<MT,unaligned,false,true,CSAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Reference
   Submatrix<MT,unaligned,false,true,CSAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstReference
   Submatrix<MT,unaligned,false,true,CSAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Pointer
   Submatrix<MT,unaligned,false,true,CSAs...>::data() noexcept
{
   return matrix_.data() + row()*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstPointer
   Submatrix<MT,unaligned,false,true,CSAs...>::data() const noexcept
{
   return matrix_.data() + row()*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of row/column \a i.
//
// \param i The row/column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row/column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Pointer
   Submatrix<MT,unaligned,false,true,CSAs...>::data( size_t i ) noexcept
{
   return matrix_.data() + (row()+i)*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of row/column \a i.
//
// \param i The row/column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in row/column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstPointer
   Submatrix<MT,unaligned,false,true,CSAs...>::data( size_t i ) const noexcept
{
   return matrix_.data() + (row()+i)*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Iterator
   Submatrix<MT,unaligned,false,true,CSAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return Iterator( matrix_.begin( row() + i ) + column(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,false,true,CSAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ConstIterator( matrix_.cbegin( row() + i ) + column(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,false,true,CSAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ConstIterator( matrix_.cbegin( row() + i ) + column(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::Iterator
   Submatrix<MT,unaligned,false,true,CSAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return Iterator( matrix_.begin( row() + i ) + column() + columns(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,false,true,CSAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ConstIterator( matrix_.cbegin( row() + i ) + column() + columns(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,false,true,CSAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ConstIterator( matrix_.cbegin( row() + i ) + column() + columns(), isAligned_ );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t iend( row() + rows() );
   decltype(auto) left( derestrict( matrix_ ) );

   for( size_t i=row(); i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j ) {
         if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, j, rhs ) )
            left(i,j) = rhs;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all submatrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid initializer list dimension.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the submatrix
// by means of an initializer list. The submatrix elements are assigned the values from the given
// initializer list. Missing values are initialized as default. Note that in case the size of the
// top-level initializer list does not match the number of rows of the submatrix or the size of
// any nested list exceeds the number of columns, a \a std::invalid_argument exception is thrown.
// Also, if the underlying matrix \a MT is restricted and the assignment would violate an
// invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to submatrix" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      if( !tryAssign( matrix_, tmp, row(), column() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), left.begin(i) ), left.end(i), ElementType() );
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Submatrix.
//
// \param rhs Dense submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::operator=( const Submatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() && column() == rhs.column() ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::operator=( const Matrix<MT2,SO2>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( right ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !tryAddAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpAddAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator+=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySubAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSubAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator-=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySchurAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSchurAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO2 >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::operator%=( const Matrix<MT2,SO2>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

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
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline MT& Submatrix<MT,unaligned,false,true,CSAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline const MT& Submatrix<MT,unaligned,false,true,CSAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,false,true,CSAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,false,true,CSAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row/column.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,false,true,CSAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,false,true,CSAs...>::nonZeros() const
{
   const size_t iend( row() + rows() );
   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t i=row(); i<iend; ++i )
      for( size_t j=column(); j<jend; ++j )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,false,true,CSAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t j=column(); j<jend; ++j )
      if( !isDefault( matrix_(row()+i,j) ) )
         ++nonzeros;

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
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,unaligned,false,true,CSAs...>::reset()
{
   using blaze::clear;

   for( size_t i=row(); i<row()+rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j )
         clear( matrix_(i,j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// This function resets the values in the specified row/column to their default value. In case
// the storage order is set to \a rowMajor the function resets the values in row \a i, in case
// the storage order is set to \a columnMajor the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,unaligned,false,true,CSAs...>::reset( size_t i )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( max( i+1UL, column() ) )
                           :( max( i, column() ) ) )
                        :( column() ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( min( i, column()+columns() ) )
                           :( min( i+1UL, column()+columns() ) ) )
                        :( column()+columns() ) );

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row()+i,j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,false,true,CSAs...>::hasOverlap() const noexcept
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric_v<MT> || IsHermitian_v<MT>, "Invalid matrix detected" );

   if( ( row() + rows() <= column() ) || ( column() + columns() <= row() ) )
      return false;
   else return true;
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
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::transpose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::ctranspose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, ctrans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( ctrans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales the submatrix by applying the given scalar value \a scalar to each
// element of the submatrix. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used
// to scale a submatrix on a lower or upper unitriangular matrix. The attempt to scale
// such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the scalar value
inline Submatrix<MT,unaligned,false,true,CSAs...>&
   Submatrix<MT,unaligned,false,true,CSAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t iend( row() + rows() );

   for( size_t i=row(); i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) *= scalar;
   }

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
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,unaligned,false,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,unaligned,false,true,CSAs...>::canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,unaligned,false,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,unaligned,false,true,CSAs...>::isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each row/column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,false,true,CSAs...>::isAligned() const noexcept
{
   return isAligned_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the submatrix).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,false,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the column index (in case of a row-major matrix) or
// the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,false,true,CSAs...>::load( size_t i, size_t j ) const noexcept
{
   if( isAligned_ )
      return loada( i, j );
   else
      return loadu( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,false,true,CSAs...>::loada( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   return matrix_.loada( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,false,true,CSAs...>::loadu( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );

   return matrix_.loadu( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than the
// number of columns. Additionally, the column index (in case of a row-major matrix) or the row
// index (in case of a column-major matrix) must be a multiple of the number of values inside
// the SIMD element. This function must \b NOT be called explicitly! It is used internally
// for the performance optimized evaluation of expression templates. Calling this function
// explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,false,true,CSAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( isAligned_ )
      storea( i, j, value );
   else
      storeu( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,false,true,CSAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   matrix_.storea( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,false,true,CSAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );

   matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must be
// smaller than the number of columns. Additionally, the column index (in case of a row-major
// matrix) or the row index (in case of a column-major matrix) must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,false,true,CSAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   if( isAligned_ )
      matrix_.stream( row()+i, column()+j, value );
   else
      matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row()+i,column()+j    ) = (*rhs)(i,j    );
         matrix_(row()+i,column()+j+1UL) = (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(row()+i,column()+jpos) = (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   if( useStreaming && isAligned_ &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; j<jpos; j+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<jpos; j+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) = (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+i,column()+i) += (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(row()+i,column()+j    ) += (*rhs)(i,j    );
            matrix_(row()+i,column()+j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(row()+i,column()+jpos) += (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) += (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+i,column()+i) -= (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(row()+i,column()+j    ) -= (*rhs)(i,j    );
            matrix_(row()+i,column()+j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(row()+i,column()+jpos) -= (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) -= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row()+i,column()+j    ) *= (*rhs)(i,j    );
         matrix_(row()+i,column()+j+1UL) *= (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(row()+i,column()+jpos) *= (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

      size_t j( 0UL );
      Iterator left( begin(i) );
      ConstIterator_t<MT2> right( (*rhs).begin(i) );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<columns(); ++j ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) *= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,false,true,CSAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR UNALIGNED COLUMN-MAJOR DENSE SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Submatrix for unaligned column-major dense submatrices.
// \ingroup submatrix
//
// This Specialization of Submatrix adapts the class template to the requirements of unaligned
// column-major dense submatrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
class Submatrix<MT,unaligned,true,true,CSAs...>
   : public View< DenseMatrix< Submatrix<MT,unaligned,true,true,CSAs...>, true > >
   , private SubmatrixData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = SubmatrixData<CSAs...>;               //!< The type of the SubmatrixData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Submatrix instance.
   using This = Submatrix<MT,unaligned,true,true,CSAs...>;

   //! Base type of this Submatrix instance.
   using BaseType = View< DenseMatrix<This,true> >;

   using ViewedType    = MT;                            //!< The type viewed by this Submatrix instance.
   using ResultType    = SubmatrixTrait_t<MT,CSAs...>;  //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;    //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;             //!< Type of the submatrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the submatrix elements.
   using ReturnType    = ReturnType_t<MT>;              //!< Return type for expression template evaluations
   using CompositeType = const Submatrix&;              //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant submatrix value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant submatrix value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant submatrix value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;
   //**********************************************************************************************

   //**SubmatrixIterator class definition**********************************************************
   /*!\brief Iterator over the elements of the sparse submatrix.
   */
   template< typename IteratorType >  // Type of the dense matrix iterator
   class SubmatrixIterator
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
      /*!\brief Default constructor of the SubmatrixIterator class.
      */
      inline SubmatrixIterator()
         : iterator_ (       )  // Iterator to the current submatrix element
         , isAligned_( false )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Constructor of the SubmatrixIterator class.
      //
      // \param iterator Iterator to the initial element.
      // \param finalIterator The final iterator for SIMD operations.
      // \param remainingElements The number of remaining elements beyond the final iterator.
      // \param isMemoryAligned Memory alignment flag.
      */
      inline SubmatrixIterator( IteratorType iterator, bool isMemoryAligned )
         : iterator_ ( iterator        )  // Iterator to the current submatrix element
         , isAligned_( isMemoryAligned )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Constructor******************************************************************************
      /*!\brief Conversion constructor from different SubmatrixIterator instances.
      //
      // \param it The submatrix iterator to be copied.
      */
      template< typename IteratorType2 >
      inline SubmatrixIterator( const SubmatrixIterator<IteratorType2>& it )
         : iterator_ ( it.base()      )  // Iterator to the current submatrix element
         , isAligned_( it.isAligned() )  // Memory alignment flag
      {}
      //*******************************************************************************************

      //**Addition assignment operator*************************************************************
      /*!\brief Addition assignment operator.
      //
      // \param inc The increment of the iterator.
      // \return The incremented iterator.
      */
      inline SubmatrixIterator& operator+=( size_t inc ) {
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
      inline SubmatrixIterator& operator-=( size_t dec ) {
         iterator_ -= dec;
         return *this;
      }
      //*******************************************************************************************

      //**Prefix increment operator****************************************************************
      /*!\brief Pre-increment operator.
      //
      // \return Reference to the incremented iterator.
      */
      inline SubmatrixIterator& operator++() {
         ++iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix increment operator***************************************************************
      /*!\brief Post-increment operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator++( int ) {
         return SubmatrixIterator( iterator_++, isAligned_ );
      }
      //*******************************************************************************************

      //**Prefix decrement operator****************************************************************
      /*!\brief Pre-decrement operator.
      //
      // \return Reference to the decremented iterator.
      */
      inline SubmatrixIterator& operator--() {
         --iterator_;
         return *this;
      }
      //*******************************************************************************************

      //**Postfix decrement operator***************************************************************
      /*!\brief Post-decrement operator.
      //
      // \return The previous position of the iterator.
      */
      inline const SubmatrixIterator operator--( int ) {
         return SubmatrixIterator( iterator_--, isAligned_ );
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return The resulting value.
      */
      inline ReferenceType operator*() const {
         return *iterator_;
      }
      //*******************************************************************************************

      //**Element access operator******************************************************************
      /*!\brief Direct access to the element at the current iterator position.
      //
      // \return Pointer to the element at the current iterator position.
      */
      inline IteratorType operator->() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**Load function****************************************************************************
      /*!\brief Load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs a load of the current SIMD element of the submatrix iterator.
      // This function must \b NOT be called explicitly! It is used internally for the performance
      // optimized evaluation of expression templates. Calling this function explicitly might
      // result in erroneous results and/or in compilation errors.
      */
      inline SIMDType load() const noexcept {
         if( isAligned_ )
            return loada();
         else
            return loadu();
      }
      //*******************************************************************************************

      //**Loada function***************************************************************************
      /*!\brief Aligned load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an aligned load of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loada() const noexcept {
         return iterator_.loada();
      }
      //*******************************************************************************************

      //**Loadu function***************************************************************************
      /*!\brief Unaligned load of a SIMD element of the dense submatrix.
      //
      // \return The loaded SIMD element.
      //
      // This function performs an unaligned load of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline SIMDType loadu() const noexcept {
         return iterator_.loadu();
      }
      //*******************************************************************************************

      //**Store function***************************************************************************
      /*!\brief Store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs a store of the current SIMD element of the submatrix iterator.
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
      /*!\brief Aligned store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned store of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storea( const SIMDType& value ) const {
         iterator_.storea( value );
      }
      //*******************************************************************************************

      //**Storeu function**************************************************************************
      /*!\brief Unaligned store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an unaligned store of the current SIMD element of the submatrix
      // iterator. This function must \b NOT be called explicitly! It is used internally for the
      // performance optimized evaluation of expression templates. Calling this function explicitly
      // might result in erroneous results and/or in compilation errors.
      */
      inline void storeu( const SIMDType& value ) const {
         iterator_.storeu( value );
      }
      //*******************************************************************************************

      //**Stream function**************************************************************************
      /*!\brief Aligned, non-temporal store of a SIMD element of the dense submatrix.
      //
      // \param value The SIMD element to be stored.
      // \return void
      //
      // This function performs an aligned, non-temporal store of the current SIMD element of the
      // submatrix iterator. This function must \b NOT be called explicitly! It is used internally
      // for the performance optimized evaluation of expression templates. Calling this function
      // explicitly might result in erroneous results and/or in compilation errors.
      */
      inline void stream( const SIMDType& value ) const {
         iterator_.stream( value );
      }
      //*******************************************************************************************

      //**Equality operator************************************************************************
      /*!\brief Equality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators refer to the same element, \a false if not.
      */
      inline bool operator==( const SubmatrixIterator& rhs ) const {
         return iterator_ == rhs.iterator_;
      }
      //*******************************************************************************************

      //**Inequality operator**********************************************************************
      /*!\brief Inequality comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the iterators don't refer to the same element, \a false if they do.
      */
      inline bool operator!=( const SubmatrixIterator& rhs ) const {
         return iterator_ != rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-than operator***********************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller, \a false if not.
      */
      inline bool operator<( const SubmatrixIterator& rhs ) const {
         return iterator_ < rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-than operator********************************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater, \a false if not.
      */
      inline bool operator>( const SubmatrixIterator& rhs ) const {
         return iterator_ > rhs.iterator_;
      }
      //*******************************************************************************************

      //**Less-or-equal-than operator**************************************************************
      /*!\brief Less-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is smaller or equal, \a false if not.
      */
      inline bool operator<=( const SubmatrixIterator& rhs ) const {
         return iterator_ <= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Greater-or-equal-than operator***********************************************************
      /*!\brief Greater-than comparison between two SubmatrixIterator objects.
      //
      // \param rhs The right-hand side iterator.
      // \return \a true if the left-hand side iterator is greater or equal, \a false if not.
      */
      inline bool operator>=( const SubmatrixIterator& rhs ) const {
         return iterator_ >= rhs.iterator_;
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Calculating the number of elements between two iterators.
      //
      // \param rhs The right-hand side iterator.
      // \return The number of elements between the two iterators.
      */
      inline DifferenceType operator-( const SubmatrixIterator& rhs ) const {
         return iterator_ - rhs.iterator_;
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be incremented.
      // \param inc The number of elements the iterator is incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( const SubmatrixIterator& it, size_t inc ) {
         return SubmatrixIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Addition operator************************************************************************
      /*!\brief Addition between an integral value and a SubmatrixIterator.
      //
      // \param inc The number of elements the iterator is incremented.
      // \param it The iterator to be incremented.
      // \return The incremented iterator.
      */
      friend inline const SubmatrixIterator operator+( size_t inc, const SubmatrixIterator& it ) {
         return SubmatrixIterator( it.iterator_ + inc, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Subtraction operator*********************************************************************
      /*!\brief Subtraction between a SubmatrixIterator and an integral value.
      //
      // \param it The iterator to be decremented.
      // \param inc The number of elements the iterator is decremented.
      // \return The decremented iterator.
      */
      friend inline const SubmatrixIterator operator-( const SubmatrixIterator& it, size_t dec ) {
         return SubmatrixIterator( it.iterator_ - dec, it.isAligned_ );
      }
      //*******************************************************************************************

      //**Base function****************************************************************************
      /*!\brief Access to the current position of the submatrix iterator.
      //
      // \return The current position of the submatrix iterator.
      */
      inline IteratorType base() const {
         return iterator_;
      }
      //*******************************************************************************************

      //**IsAligned function***********************************************************************
      /*!\brief Access to the iterator's memory alignment flag.
      //
      // \return \a true in case the iterator is aligned, \a false if it is not.
      */
      inline bool isAligned() const noexcept {
         return isAligned_;
      }
      //*******************************************************************************************

    private:
      //**Member variables*************************************************************************
      IteratorType iterator_;   //!< Iterator to the current submatrix element.
      bool         isAligned_;  //!< Memory alignment flag.
      //*******************************************************************************************
   };
   //**********************************************************************************************

   //**Type definitions****************************************************************************
   //! Iterator over constant elements.
   using ConstIterator = SubmatrixIterator< ConstIterator_t<MT> >;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, SubmatrixIterator< Iterator_t<MT> > >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RSAs >
   explicit inline Submatrix( MT& matrix, RSAs... args );

   Submatrix( const Submatrix& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Submatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Submatrix& operator=( const ElementType& rhs );
   inline Submatrix& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Submatrix& operator=( const Submatrix& rhs );

   template< typename MT2, bool SO >
   inline Submatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;
   using DataType::column;
   using DataType::rows;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Submatrix& transpose();
   inline Submatrix& ctranspose();

   template< typename Other > inline Submatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<MT2> > );
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

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;        //!< The matrix containing the submatrix.
   const bool isAligned_;  //!< Memory alignment flag.
                           /*!< The alignment flag indicates whether the submatrix is fully aligned
                                with respect to the given element type and the available instruction
                                set. In case the submatrix is fully aligned it is possible to use
                                aligned loads and stores instead of unaligned loads and stores. In
                                order to be aligned, the first element of each row/column must be
                                aligned. */
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, AlignmentFlag AF2, bool SO2, bool DF2, size_t... CSAs2 > friend class Submatrix;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
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
/*!\brief The constructor for unaligned column-major dense submatrices.
//
// \param matrix The dense matrix containing the submatrix.
// \param args The runtime submatrix arguments.
// \exception std::invalid_argument Invalid submatrix specification.
//
// By default, the provided submatrix arguments are checked at runtime. In case the submatrix is
// not properly specified (i.e. if the specified submatrix is not contained in the given dense
// matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by providing
// the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CSAs >    // Compile time submatrix arguments
template< typename... RSAs >  // Runtime submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>::Submatrix( MT& matrix, RSAs... args )
   : DataType  ( args... )  // Base class initialization
   , matrix_   ( matrix  )  // The matrix containing the submatrix
   , isAligned_( simdEnabled && IsContiguous_v<MT> &&
                 matrix.data() != nullptr && checkAlignment( data() ) &&
                 ( columns() < 2UL || ( matrix.spacing() % SIMDSIZE ) == 0UL ) )
{
   if( isChecked( args... ) ) {
      if( ( row() + rows() > matrix_.rows() ) || ( column() + columns() > matrix_.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }
   }
   else {
      BLAZE_USER_ASSERT( row()    + rows()    <= matrix_.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column() + columns() <= matrix_.columns(), "Invalid submatrix specification" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Reference
   Submatrix<MT,unaligned,true,true,CSAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstReference
   Submatrix<MT,unaligned,true,true,CSAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Reference
   Submatrix<MT,unaligned,true,true,CSAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstReference
   Submatrix<MT,unaligned,true,true,CSAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Pointer
   Submatrix<MT,unaligned,true,true,CSAs...>::data() noexcept
{
   return matrix_.data() + row() + column()*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstPointer
   Submatrix<MT,unaligned,true,true,CSAs...>::data() const noexcept
{
   return matrix_.data() + row() + column()*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Pointer
   Submatrix<MT,unaligned,true,true,CSAs...>::data( size_t j ) noexcept
{
   return matrix_.data() + row() + (column()+j)*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstPointer
   Submatrix<MT,unaligned,true,true,CSAs...>::data( size_t j ) const noexcept
{
   return matrix_.data() + row() + (column()+j)*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Iterator
   Submatrix<MT,unaligned,true,true,CSAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return Iterator( matrix_.begin( column() + j ) + row(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,true,true,CSAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ConstIterator( matrix_.cbegin( column() + j ) + row(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,true,true,CSAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ConstIterator( matrix_.cbegin( column() + j ) + row(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::Iterator
   Submatrix<MT,unaligned,true,true,CSAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return Iterator( matrix_.begin( column() + j ) + row() + rows(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,true,true,CSAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ConstIterator( matrix_.cbegin( column() + j ) + row() + rows(), isAligned_ );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,unaligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,unaligned,true,true,CSAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ConstIterator( matrix_.cbegin( column() + j ) + row() + rows(), isAligned_ );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t jend( column() + columns() );
   decltype(auto) left( derestrict( matrix_ ) );

   for( size_t j=column(); j<jend; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i ) {
         if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, j, rhs ) )
            left(i,j) = rhs;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all submatrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to submatrix.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the submatrix
// by means of an initializer list. The submatrix elements are assigned the values from the given
// initializer list. Missing values are initialized as default. Note that in case the size of the
// top-level initializer list does not match the number of rows of the submatrix or the size of
// any nested list exceeds the number of columns, a \a std::invalid_argument exception is thrown.
// Also, if the underlying matrix \a MT is restricted and the assignment would violate an
// invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   using blaze::reset;

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to submatrix" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      if( !tryAssign( matrix_, tmp, row(), column() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      size_t j( 0UL );
      for( const auto& element : rowList ) {
         left(i,j) = element;
         ++j;
      }
      for( ; j<columns(); ++j ) {
         reset( left(i,j) );
      }
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Submatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::operator=( const Submatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() && column() == rhs.column() ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( right ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !tryAddAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpAddAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySubAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSubAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySchurAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSchurAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

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
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline MT& Submatrix<MT,unaligned,true,true,CSAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline const MT& Submatrix<MT,unaligned,true,true,CSAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total
// number of elements of a column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,true,true,CSAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,true,true,CSAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,true,true,CSAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,true,true,CSAs...>::nonZeros() const
{
   const size_t iend( row() + rows() );
   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t j=column(); j<jend; ++j )
      for( size_t i=row(); i<iend; ++i )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,unaligned,true,true,CSAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( row() + rows() );
   size_t nonzeros( 0UL );

   for( size_t i=row(); i<iend; ++i )
      if( !isDefault( matrix_(i,column()+j) ) )
         ++nonzeros;

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
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,unaligned,true,true,CSAs...>::reset()
{
   using blaze::clear;

   for( size_t j=column(); j<column()+columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i )
         clear( matrix_(i,j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,unaligned,true,true,CSAs...>::reset( size_t j )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( max( j+1UL, row() ) )
                           :( max( j, row() ) ) )
                        :( row() ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( min( j, row()+rows() ) )
                           :( min( j+1UL, row()+rows() ) ) )
                        :( row()+rows() ) );

   for( size_t i=ibegin; i<iend; ++i )
      clear( matrix_(i,column()+j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,true,true,CSAs...>::hasOverlap() const noexcept
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric_v<MT> || IsHermitian_v<MT>, "Invalid matrix detected" );

   if( ( row() + rows() <= column() ) || ( column() + columns() <= row() ) )
      return false;
   else return true;
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
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::transpose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::ctranspose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, ctrans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( ctrans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales the submatrix by applying the given scalar value \a scalar to each
// element of the submatrix. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used
// to scale a submatrix on a lower or upper unitriangular matrix. The attempt to scale
// such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the scalar value
inline Submatrix<MT,unaligned,true,true,CSAs...>&
   Submatrix<MT,unaligned,true,true,CSAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jend( column() + columns() );

   for( size_t j=column(); j<jend; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) *= scalar;
   }

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
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,unaligned,true,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,unaligned,true,true,CSAs...>::canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,unaligned,true,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,unaligned,true,true,CSAs...>::isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,true,true,CSAs...>::isAligned() const noexcept
{
   return isAligned_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the submatrix).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,unaligned,true,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,true,true,CSAs...>::load( size_t i, size_t j ) const noexcept
{
   if( isAligned_ )
      return loada( i, j );
   else
      return loadu( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,true,true,CSAs...>::loada( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.loada( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,unaligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,unaligned,true,true,CSAs...>::loadu( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.loadu( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense submatrix. The
// row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,true,true,CSAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   if( isAligned_ )
      storea( i, j, value );
   else
      storeu( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,true,true,CSAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.storea( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,true,true,CSAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the row index must be a
// multiple of the number of values inside the SIMD element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,unaligned,true,true,CSAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   if( isAligned_ )
      matrix_.stream( row()+i, column()+j, value );
   else
      matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row()+i    ,column()+j) = (*rhs)(i    ,j);
         matrix_(row()+i+1UL,column()+j) = (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(row()+ipos,column()+j) = (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   if( useStreaming && isAligned_ &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; i<ipos; i+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) = (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+j,column()+j) += (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(row()+i    ,column()+j) += (*rhs)(i    ,j);
            matrix_(row()+i+1UL,column()+j) += (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(row()+ipos,column()+j) += (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) += (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+j,column()+j) -= (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(row()+i    ,column()+j) -= (*rhs)(i    ,j);
            matrix_(row()+i+1UL,column()+j) -= (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(row()+ipos,column()+j) -= (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) -= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row()+i    ,column()+j) *= (*rhs)(i    ,j);
         matrix_(row()+i+1UL,column()+j) *= (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(row()+ipos,column()+j) *= (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a column-major dense
//        matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,unaligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

      size_t i( 0UL );
      Iterator left( begin(j) );
      ConstIterator_t<MT2> right( (*rhs).begin(j) );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<rows(); ++i ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) *= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,unaligned,true,true,CSAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED ROW-MAJOR DENSE SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Submatrix for aligned row-major dense submatrices.
// \ingroup submatrix
//
// This Specialization of Submatrix adapts the class template to the requirements of aligned
// row-major dense submatrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
class Submatrix<MT,aligned,false,true,CSAs...>
   : public View< DenseMatrix< Submatrix<MT,aligned,false,true,CSAs...>, false > >
   , private SubmatrixData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = SubmatrixData<CSAs...>;               //!< The type of the SubmatrixData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Submatrix instance.
   using This = Submatrix<MT,aligned,false,true,CSAs...>;

   //! Base type of this Submatrix instance.
   using BaseType = View< DenseMatrix<This,false> >;

   using ViewedType    = MT;                            //!< The type viewed by this Submatrix instance.
   using ResultType    = SubmatrixTrait_t<MT,CSAs...>;  //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;    //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;             //!< Type of the submatrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the submatrix elements.
   using ReturnType    = ReturnType_t<MT>;              //!< Return type for expression template evaluations
   using CompositeType = const Submatrix&;              //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant submatrix value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant submatrix value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant submatrix value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RSAs >
   explicit inline Submatrix( MT& matrix, RSAs... args );

   Submatrix( const Submatrix& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Submatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t i ) noexcept;
   inline ConstPointer   data  ( size_t i ) const noexcept;
   inline Iterator       begin ( size_t i );
   inline ConstIterator  begin ( size_t i ) const;
   inline ConstIterator  cbegin( size_t i ) const;
   inline Iterator       end   ( size_t i );
   inline ConstIterator  end   ( size_t i ) const;
   inline ConstIterator  cend  ( size_t i ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Submatrix& operator=( const ElementType& rhs );
   inline Submatrix& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Submatrix& operator=( const Submatrix& rhs );

   template< typename MT2, bool SO >
   inline Submatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;
   using DataType::column;
   using DataType::rows;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Submatrix& transpose();
   inline Submatrix& ctranspose();

   template< typename Other > inline Submatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<MT2> > );
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

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>& rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,false>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, AlignmentFlag AF2, bool SO2, bool DF2, size_t... CSAs2 > friend class Submatrix;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE   ( MT );
   BLAZE_CONSTRAINT_MUST_BE_ROW_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE   ( MT );
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
/*!\brief Constructor for aligned row-major dense submatrices.
//
// \param matrix The dense matrix containing the submatrix.
// \param args The runtime submatrix arguments.
// \exception std::invalid_argument Invalid submatrix specification.
//
// By default, the provided submatrix arguments are checked at runtime. In case the submatrix is
// not properly specified (i.e. if the specified submatrix is not contained in the given dense
// matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by providing
// the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CSAs >    // Compile time submatrix arguments
template< typename... RSAs >  // Runtime submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>::Submatrix( MT& matrix, RSAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the submatrix
{
   if( isChecked( args... ) )
   {
      if( ( row() + rows() > matrix_.rows() ) || ( column() + columns() > matrix_.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }

      if( simdEnabled && IsContiguous_v<MT> &&
          ( !checkAlignment( data() ) ||
            ( rows() > 1UL && matrix_.spacing() % SIMDSIZE != 0UL ) ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix alignment" );
      }
   }
   else
   {
      BLAZE_USER_ASSERT( row()    + rows()    <= matrix_.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column() + columns() <= matrix_.columns(), "Invalid submatrix specification" );

      BLAZE_USER_ASSERT( !simdEnabled || !IsContiguous_v<MT> || checkAlignment( data() ), "Invalid submatrix alignment" );
      BLAZE_USER_ASSERT( !simdEnabled || !IsContiguous_v<MT> || rows() <= 1UL || matrix_.spacing() % SIMDSIZE == 0UL, "Invalid submatrix alignment" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Reference
   Submatrix<MT,aligned,false,true,CSAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstReference
   Submatrix<MT,aligned,false,true,CSAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Reference
   Submatrix<MT,aligned,false,true,CSAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstReference
   Submatrix<MT,aligned,false,true,CSAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Pointer
   Submatrix<MT,aligned,false,true,CSAs...>::data() noexcept
{
   return matrix_.data() + row()*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstPointer
   Submatrix<MT,aligned,false,true,CSAs...>::data() const noexcept
{
   return matrix_.data() + row()*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix in row \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Pointer
   Submatrix<MT,aligned,false,true,CSAs...>::data( size_t i ) noexcept
{
   return matrix_.data() + (row()+i)*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of row \a i.
//
// \param i The row index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix in row \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstPointer
   Submatrix<MT,aligned,false,true,CSAs...>::data( size_t i ) const noexcept
{
   return matrix_.data() + (row()+i)*spacing() + column();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Iterator
   Submatrix<MT,aligned,false,true,CSAs...>::begin( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.begin( row() + i ) + column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,false,true,CSAs...>::begin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row() + i ) + column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator to the first non-zero element of row/column \a i.
//
// This function returns a row/column iterator to the first non-zero element of row/column \a i.
// In case the storage order is set to \a rowMajor the function returns an iterator to the first
// non-zero element of row \a i, in case the storage flag is set to \a columnMajor the function
// returns an iterator to the first non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,false,true,CSAs...>::cbegin( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row() + i ) + column() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::Iterator
   Submatrix<MT,aligned,false,true,CSAs...>::end( size_t i )
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.begin( row() + i ) + column() + columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,false,true,CSAs...>::end( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row() + i ) + column() + columns() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of row/column \a i.
//
// \param i The row/column index.
// \return Iterator just past the last non-zero element of row/column \a i.
//
// This function returns an row/column iterator just past the last non-zero element of row/column
// \a i. In case the storage order is set to \a rowMajor the function returns an iterator just
// past the last non-zero element of row \a i, in case the storage flag is set to \a columnMajor
// the function returns an iterator just past the last non-zero element of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,false,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,false,true,CSAs...>::cend( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid dense submatrix row access index" );
   return ( matrix_.cbegin( row() + i ) + column() + columns() );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t iend( row() + rows() );
   decltype(auto) left( derestrict( matrix_ ) );

   for( size_t i=row(); i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j ) {
         if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, j, rhs ) )
            left(i,j) = rhs;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all submatrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to submatrix.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the submatrix
// by means of an initializer list. The submatrix elements are assigned the values from the given
// initializer list. Missing values are initialized as default. Note that in case the size of the
// top-level initializer list does not match the number of rows of the submatrix or the size of
// any nested list exceeds the number of columns, a \a std::invalid_argument exception is thrown.
// Also, if the underlying matrix \a MT is restricted and the assignment would violate an
// invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to submatrix" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      if( !tryAssign( matrix_, tmp, row(), column() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      std::fill( std::copy( rowList.begin(), rowList.end(), left.begin(i) ), left.end(i), ElementType() );
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Submatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::operator=( const Submatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() && column() == rhs.column() ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given matrix. In case the current sizes
// of the two matrices don't match, a \a std::invalid_argument exception is thrown. Also, if
// the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( right ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !tryAddAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpAddAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySubAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSubAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySchurAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSchurAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

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
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline MT& Submatrix<MT,aligned,false,true,CSAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline const MT& Submatrix<MT,aligned,false,true,CSAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two rows/columns.
//
// \return The spacing between the beginning of two rows/columns.
//
// This function returns the spacing between the beginning of two rows/columns, i.e. the
// total number of elements of a row/column. In case the storage order is set to \a rowMajor
// the function returns the spacing between two rows, in case the storage flag is set to
// \a columnMajor the function returns the spacing between two columns.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,false,true,CSAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,false,true,CSAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified row/column.
//
// \param i The index of the row/column.
// \return The current capacity of row/column \a i.
//
// This function returns the current capacity of the specified row/column. In case the
// storage order is set to \a rowMajor the function returns the capacity of row \a i,
// in case the storage flag is set to \a columnMajor the function returns the capacity
// of column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,false,true,CSAs...>::capacity( size_t i ) const noexcept
{
   MAYBE_UNUSED( i );

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   return columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,false,true,CSAs...>::nonZeros() const
{
   const size_t iend( row() + rows() );
   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t i=row(); i<iend; ++i )
      for( size_t j=column(); j<jend; ++j )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified row/column.
//
// \param i The index of the row/column.
// \return The number of non-zero elements of row/column \a i.
//
// This function returns the current number of non-zero elements in the specified row/column.
// In case the storage order is set to \a rowMajor the function returns the number of non-zero
// elements in row \a i, in case the storage flag is set to \a columnMajor the function returns
// the number of non-zero elements in column \a i.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,false,true,CSAs...>::nonZeros( size_t i ) const
{
   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t j=column(); j<jend; ++j )
      if( !isDefault( matrix_(row()+i,j) ) )
         ++nonzeros;

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
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,aligned,false,true,CSAs...>::reset()
{
   using blaze::clear;

   for( size_t i=row(); i<row()+rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j )
         clear( matrix_(i,j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified row/column to the default initial values.
//
// \param i The index of the row/column.
// \return void
//
// This function resets the values in the specified row/column to their default value. In case
// the storage order is set to \a rowMajor the function resets the values in row \a i, in case
// the storage order is set to \a columnMajor the function resets the values in column \a i.
// Note that the capacity of the row/column remains unchanged.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,aligned,false,true,CSAs...>::reset( size_t i )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( i < rows(), "Invalid row access index" );

   const size_t jbegin( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( max( i+1UL, column() ) )
                           :( max( i, column() ) ) )
                        :( column() ) );
   const size_t jend  ( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( min( i, column()+columns() ) )
                           :( min( i+1UL, column()+columns() ) ) )
                        :( column()+columns() ) );

   for( size_t j=jbegin; j<jend; ++j )
      clear( matrix_(row()+i,j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,false,true,CSAs...>::hasOverlap() const noexcept
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric_v<MT> || IsHermitian_v<MT>, "Invalid matrix detected" );

   if( ( row() + rows() <= column() ) || ( column() + columns() <= row() ) )
      return false;
   else return true;
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
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::transpose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::ctranspose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, ctrans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( ctrans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales the submatrix by applying the given scalar value \a scalar to each
// element of the submatrix. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used
// to scale a submatrix on a lower or upper unitriangular matrix. The attempt to scale
// such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the scalar value
inline Submatrix<MT,aligned,false,true,CSAs...>&
   Submatrix<MT,aligned,false,true,CSAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t iend( row() + rows() );

   for( size_t i=row(); i<iend; ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT> )
                           ?( ( IsStrictlyUpper_v<MT> )
                              ?( max( i+1UL, column() ) )
                              :( max( i, column() ) ) )
                           :( column() ) );
      const size_t jend  ( ( IsLower_v<MT> )
                           ?( ( IsStrictlyLower_v<MT> )
                              ?( min( i, column()+columns() ) )
                              :( min( i+1UL, column()+columns() ) ) )
                           :( column()+columns() ) );

      for( size_t j=jbegin; j<jend; ++j )
         matrix_(i,j) *= scalar;
   }

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
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,aligned,false,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,aligned,false,true,CSAs...>::canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,aligned,false,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,aligned,false,true,CSAs...>::isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each row of the submatrix are guaranteed to conform
// to the alignment restrictions of the underlying element type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,false,true,CSAs...>::isAligned() const noexcept
{
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the submatrix).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,false,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the column index (in case of a row-major matrix) or
// the row index (in case of a column-major matrix) must be a multiple of the number of values
// inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,false,true,CSAs...>::load( size_t i, size_t j ) const noexcept
{
   return loada( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,false,true,CSAs...>::loada( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   return matrix_.loada( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,false,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,false,true,CSAs...>::loadu( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );

   return matrix_.loadu( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense submatrix. The
// row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the column index (in case of a row-major matrix)
// or the row index (in case of a column-major matrix) must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,false,true,CSAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   return storea( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the column index (in case of a row-major matrix) or the
// row index (in case of a column-major matrix) must be a multiple of the number of values inside
// the SIMD element. This function must \b NOT be called explicitly! It is used internally for
// the performance optimized evaluation of expression templates. Calling this function explicitly
// might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,false,true,CSAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   return matrix_.storea( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense
// submatrix. The row index must be smaller than the number of rows and the column index must
// be smaller than the number of columns. Additionally, the column index (in case of a row-major
// matrix) or the row index (in case of a column-major matrix) must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,false,true,CSAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );

   matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense submatrix. The row index must be smaller than the number of rows and the column index
// must be smaller than the number of columns. Additionally, the column index (in case of a
// row-major matrix) or the row index (in case of a column-major matrix) must be a multiple
// of the number of values inside the SIMD element. This function must \b NOT be called
// explicitly! It is used internally for the performance optimized evaluation of expression
// templates. Calling this function explicitly might result in erroneous results and/or in
// compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,false,true,CSAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j + SIMDSIZE <= columns(), "Invalid column access index" );
   BLAZE_INTERNAL_ASSERT( j % SIMDSIZE == 0UL, "Invalid column access index" );

   matrix_.stream( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row()+i,column()+j    ) = (*rhs)(i,j    );
         matrix_(row()+i,column()+j+1UL) = (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(row()+i,column()+jpos) = (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   if( useStreaming &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; j<jpos; j+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t i=0UL; i<rows(); ++i )
      {
         size_t j( 0UL );
         Iterator left( begin(i) );
         ConstIterator_t<MT2> right( (*rhs).begin(i) );

         for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<jpos; j+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; j<columns(); ++j ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) = (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+i,column()+i) += (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(row()+i,column()+j    ) += (*rhs)(i,j    );
            matrix_(row()+i,column()+j+1UL) += (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(row()+i,column()+jpos) += (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) += (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+i,column()+i) -= (*rhs)(i,i);
      }
      else {
         for( size_t j=0UL; j<jpos; j+=2UL ) {
            matrix_(row()+i,column()+j    ) -= (*rhs)(i,j    );
            matrix_(row()+i,column()+j+1UL) -= (*rhs)(i,j+1UL);
         }
         if( jpos < columns() ) {
            matrix_(row()+i,column()+jpos) -= (*rhs)(i,jpos);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jbegin( ( IsUpper_v<MT2> )
                           ?( prevMultiple( ( IsStrictlyUpper_v<MT2> ? i+1UL : i ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t jend  ( ( IsLower_v<MT2> )
                           ?( IsStrictlyLower_v<MT2> ? i : i+1UL )
                           :( columns() ) );
      BLAZE_INTERNAL_ASSERT( jbegin <= jend, "Invalid loop indices detected" );

      const size_t jpos( prevMultiple( jend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= jend, "Invalid end calculation" );

      size_t j( jbegin );
      Iterator left( begin(i) + jbegin );
      ConstIterator_t<MT2> right( (*rhs).begin(i) + jbegin );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jend; ++j ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) -= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t jpos( prevMultiple( columns(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

   for( size_t i=0UL; i<rows(); ++i ) {
      for( size_t j=0UL; j<jpos; j+=2UL ) {
         matrix_(row()+i,column()+j    ) *= (*rhs)(i,j    );
         matrix_(row()+i,column()+j+1UL) *= (*rhs)(i,j+1UL);
      }
      if( jpos < columns() ) {
         matrix_(row()+i,column()+jpos) *= (*rhs)(i,jpos);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      const size_t jpos( prevMultiple( columns(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( jpos <= columns(), "Invalid end calculation" );

      size_t j( 0UL );
      Iterator left( begin(i) );
      ConstIterator_t<MT2> right( (*rhs).begin(i) );

      for( ; (j+SIMDSIZE*3UL) < jpos; j+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<jpos; j+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; j<columns(); ++j ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t ii=0UL; ii<rows(); ii+=block ) {
      const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
      for( size_t jj=0UL; jj<columns(); jj+=block ) {
         const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
         for( size_t i=ii; i<iend; ++i ) {
            for( size_t j=jj; j<jend; ++j ) {
               matrix_(row()+i,column()+j) *= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,false,true,CSAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+element->index(),column()+j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************








//=================================================================================================
//
//  CLASS TEMPLATE SPECIALIZATION FOR ALIGNED COLUMN-MAJOR DENSE SUBMATRICES
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Specialization of Submatrix for aligned column-major dense submatrices.
// \ingroup submatrix
//
// This Specialization of Submatrix adapts the class template to the requirements of aligned
// column-major dense submatrices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
class Submatrix<MT,aligned,true,true,CSAs...>
   : public View< DenseMatrix< Submatrix<MT,aligned,true,true,CSAs...>, true > >
   , private SubmatrixData<CSAs...>
{
 private:
   //**Type definitions****************************************************************************
   using DataType = SubmatrixData<CSAs...>;               //!< The type of the SubmatrixData base class.
   using Operand  = If_t< IsExpression_v<MT>, MT, MT& >;  //!< Composite data type of the matrix expression.
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT1, typename MT2 >
   static constexpr bool EnforceEvaluation_v =
      ( IsRestricted_v<MT1> && RequiresEvaluation_v<MT2> );
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   //! Type of this Submatrix instance.
   using This = Submatrix<MT,aligned,true,true,CSAs...>;

   //! Base type of this Submatrix instance.
   using BaseType = View< DenseMatrix<This,true> >;

   using ViewedType    = MT;                            //!< The type viewed by this Submatrix instance.
   using ResultType    = SubmatrixTrait_t<MT,CSAs...>;  //!< Result type for expression template evaluations.
   using OppositeType  = OppositeType_t<ResultType>;    //!< Result type with opposite storage order for expression template evaluations.
   using TransposeType = TransposeType_t<ResultType>;   //!< Transpose type for expression template evaluations.
   using ElementType   = ElementType_t<MT>;             //!< Type of the submatrix elements.
   using SIMDType      = SIMDTrait_t<ElementType>;      //!< SIMD type of the submatrix elements.
   using ReturnType    = ReturnType_t<MT>;              //!< Return type for expression template evaluations
   using CompositeType = const Submatrix&;              //!< Data type for composite expression templates.

   //! Reference to a constant submatrix value.
   using ConstReference = ConstReference_t<MT>;

   //! Reference to a non-constant submatrix value.
   using Reference = If_t< IsConst_v<MT>, ConstReference, Reference_t<MT> >;

   //! Pointer to a constant submatrix value.
   using ConstPointer = ConstPointer_t<MT>;

   //! Pointer to a non-constant submatrix value.
   using Pointer = If_t< IsConst_v<MT> || !HasMutableDataAccess_v<MT>, ConstPointer, Pointer_t<MT> >;

   //! Iterator over constant elements.
   using ConstIterator = ConstIterator_t<MT>;

   //! Iterator over non-constant elements.
   using Iterator = If_t< IsConst_v<MT>, ConstIterator, Iterator_t<MT> >;
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool simdEnabled = MT::simdEnabled;

   //! Compilation switch for the expression template assignment strategy.
   static constexpr bool smpAssignable = MT::smpAssignable;

   //! Compilation switch for the expression template evaluation strategy.
   static constexpr bool compileTimeArgs = DataType::compileTimeArgs;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   template< typename... RSAs >
   explicit inline Submatrix( MT& matrix, RSAs... args );

   Submatrix( const Submatrix& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~Submatrix() = default;
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator()( size_t i, size_t j );
   inline ConstReference operator()( size_t i, size_t j ) const;
   inline Reference      at( size_t i, size_t j );
   inline ConstReference at( size_t i, size_t j ) const;
   inline Pointer        data  () noexcept;
   inline ConstPointer   data  () const noexcept;
   inline Pointer        data  ( size_t j ) noexcept;
   inline ConstPointer   data  ( size_t j ) const noexcept;
   inline Iterator       begin ( size_t j );
   inline ConstIterator  begin ( size_t j ) const;
   inline ConstIterator  cbegin( size_t j ) const;
   inline Iterator       end   ( size_t j );
   inline ConstIterator  end   ( size_t j ) const;
   inline ConstIterator  cend  ( size_t j ) const;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline Submatrix& operator=( const ElementType& rhs );
   inline Submatrix& operator=( initializer_list< initializer_list<ElementType> > list );
   inline Submatrix& operator=( const Submatrix& rhs );

   template< typename MT2, bool SO >
   inline Submatrix& operator=( const Matrix<MT2,SO>& rhs );

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator+=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator-=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;

   template< typename MT2, bool SO >
   inline auto operator%=( const Matrix<MT2,SO>& rhs )
      -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   using DataType::row;
   using DataType::column;
   using DataType::rows;
   using DataType::columns;

   inline MT&       operand() noexcept;
   inline const MT& operand() const noexcept;

   inline size_t spacing() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t capacity( size_t i ) const noexcept;
   inline size_t nonZeros() const;
   inline size_t nonZeros( size_t i ) const;
   inline void   reset();
   inline void   reset( size_t i );
   //@}
   //**********************************************************************************************

   //**Numeric functions***************************************************************************
   /*!\name Numeric functions */
   //@{
   inline Submatrix& transpose();
   inline Submatrix& ctranspose();

   template< typename Other > inline Submatrix& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

 private:
   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAssign_v =
      ( useOptimizedKernels &&
        simdEnabled && MT2::simdEnabled &&
        IsSIMDCombinable_v< ElementType, ElementType_t<MT2> > );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedAddAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDAdd_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSubAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDSub_v< ElementType, ElementType_t<MT2> > &&
        !IsDiagonal_v<MT2> );
   //**********************************************************************************************

   //**********************************************************************************************
   //! Helper variable template for the explicit application of the SFINAE principle.
   template< typename MT2 >
   static constexpr bool VectorizedSchurAssign_v =
      ( VectorizedAssign_v<MT2> &&
        HasSIMDMult_v< ElementType, ElementType_t<MT2> > );
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

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   template< typename Other >
   inline bool isAliased( const Other* alias ) const noexcept;

   template< typename MT2, AlignmentFlag AF2, bool SO2, size_t... CSAs2 >
   inline bool isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept;

   inline bool isAligned   () const noexcept;
   inline bool canSMPAssign() const noexcept;

   BLAZE_ALWAYS_INLINE SIMDType load ( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loada( size_t i, size_t j ) const noexcept;
   BLAZE_ALWAYS_INLINE SIMDType loadu( size_t i, size_t j ) const noexcept;

   BLAZE_ALWAYS_INLINE void store ( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storea( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void storeu( size_t i, size_t j, const SIMDType& value ) noexcept;
   BLAZE_ALWAYS_INLINE void stream( size_t i, size_t j, const SIMDType& value ) noexcept;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 >
   inline auto assign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAssign_v<MT2> >;

   template< typename MT2 > inline void assign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void assign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 >
   inline auto addAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedAddAssign_v<MT2> >;

   template< typename MT2 > inline void addAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void addAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 >
   inline auto subAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSubAssign_v<MT2> >;

   template< typename MT2 > inline void subAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void subAssign( const SparseMatrix<MT2,false>& rhs );

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> DisableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 >
   inline auto schurAssign( const DenseMatrix<MT2,true>& rhs ) -> EnableIf_t< VectorizedSchurAssign_v<MT2> >;

   template< typename MT2 > inline void schurAssign( const DenseMatrix<MT2,false>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,true>&  rhs );
   template< typename MT2 > inline void schurAssign( const SparseMatrix<MT2,false>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline bool hasOverlap() const noexcept;
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Operand matrix_;  //!< The matrix containing the submatrix.
   //@}
   //**********************************************************************************************

   //**Friend declarations*************************************************************************
   template< typename MT2, AlignmentFlag AF2, bool SO2, bool DF2, size_t... CSAs2 > friend class Submatrix;
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSEXPR_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SUBMATRIX_TYPE      ( MT );
   BLAZE_CONSTRAINT_MUST_BE_COLUMN_MAJOR_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE        ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE      ( MT );
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
/*!\brief Constructor for aligned column-major dense submatrices.
//
// \param matrix The dense matrix containing the submatrix.
// \param args The runtime submatrix arguments.
// \exception std::invalid_argument Invalid submatrix specification.
//
// By default, the provided submatrix arguments are checked at runtime. In case the submatrix is
// not properly specified (i.e. if the specified submatrix is not contained in the given dense
// matrix) a \a std::invalid_argument exception is thrown. The checks can be skipped by providing
// the optional \a blaze::unchecked argument.
*/
template< typename MT         // Type of the dense matrix
        , size_t... CSAs >    // Compile time submatrix arguments
template< typename... RSAs >  // Runtime submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>::Submatrix( MT& matrix, RSAs... args )
   : DataType( args... )  // Base class initialization
   , matrix_ ( matrix  )  // The matrix containing the submatrix
{
   if( isChecked( args... ) )
   {
      if( ( row() + rows() > matrix_.rows() ) || ( column() + columns() > matrix_.columns() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix specification" );
      }

      if( simdEnabled && IsContiguous_v<MT> &&
          ( !checkAlignment( data() ) ||
            ( columns() > 1UL && matrix_.spacing() % SIMDSIZE != 0UL ) ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid submatrix alignment" );
      }
   }
   else
   {
      BLAZE_USER_ASSERT( row()    + rows()    <= matrix_.rows()   , "Invalid submatrix specification" );
      BLAZE_USER_ASSERT( column() + columns() <= matrix_.columns(), "Invalid submatrix specification" );

      BLAZE_USER_ASSERT( !simdEnabled || !IsContiguous_v<MT> || checkAlignment( data() ), "Invalid submatrix alignment" );
      BLAZE_USER_ASSERT( !simdEnabled || !IsContiguous_v<MT> || columns() <= 1UL || matrix_.spacing() % SIMDSIZE == 0UL, "Invalid submatrix alignment" );
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
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Reference
   Submatrix<MT,aligned,true,true,CSAs...>::operator()( size_t i, size_t j )
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief 2D-access to the dense submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function only performs an index check in case BLAZE_USER_ASSERT() is active. In contrast,
// the at() function is guaranteed to perform a check of the given access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstReference
   Submatrix<MT,aligned,true,true,CSAs...>::operator()( size_t i, size_t j ) const
{
   BLAZE_USER_ASSERT( i < rows()   , "Invalid row access index"    );
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return const_cast<const MT&>( matrix_ )(row()+i,column()+j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Reference
   Submatrix<MT,aligned,true,true,CSAs...>::at( size_t i, size_t j )
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checked access to the submatrix elements.
//
// \param i Access index for the row. The index has to be in the range \f$[0..M-1]\f$.
// \param j Access index for the column. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid matrix access index.
//
// In contrast to the function call operator this function always performs a check of the given
// access indices.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstReference
   Submatrix<MT,aligned,true,true,CSAs...>::at( size_t i, size_t j ) const
{
   if( i >= rows() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid row access index" );
   }
   if( j >= columns() ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid column access index" );
   }
   return (*this)(i,j);
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Pointer
   Submatrix<MT,aligned,true,true,CSAs...>::data() noexcept
{
   return matrix_.data() + row() + column()*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements.
//
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage of the dense submatrix. Note that
// you can NOT assume that all matrix elements lie adjacent to each other! The dense submatrix
// may use techniques such as padding to improve the alignment of the data.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstPointer
   Submatrix<MT,aligned,true,true,CSAs...>::data() const noexcept
{
   return matrix_.data() + row() + column()*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Pointer
   Submatrix<MT,aligned,true,true,CSAs...>::data( size_t j ) noexcept
{
   return matrix_.data() + row() + (column()+j)*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Low-level data access to the submatrix elements of column \a j.
//
// \param j The column index.
// \return Pointer to the internal element storage.
//
// This function returns a pointer to the internal storage for the elements in column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstPointer
   Submatrix<MT,aligned,true,true,CSAs...>::data( size_t j ) const noexcept
{
   return matrix_.data() + row() + (column()+j)*spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Iterator
   Submatrix<MT,aligned,true,true,CSAs...>::begin( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.begin( column() + j ) + row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,true,true,CSAs...>::begin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column() + j ) + row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator to the first non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator to the first non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,true,true,CSAs...>::cbegin( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column() + j ) + row() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::Iterator
   Submatrix<MT,aligned,true,true,CSAs...>::end( size_t j )
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.begin( column() + j ) + row() + rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,true,true,CSAs...>::end( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column() + j ) + row() + rows() );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns an iterator just past the last non-zero element of column \a j.
//
// \param j The column index.
// \return Iterator just past the last non-zero element of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline typename Submatrix<MT,aligned,true,true,CSAs...>::ConstIterator
   Submatrix<MT,aligned,true,true,CSAs...>::cend( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid dense submatrix column access index" );
   return ( matrix_.cbegin( column() + j ) + row() + rows() );
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
/*!\brief Homogenous assignment to all submatrix elements.
//
// \param rhs Scalar value to be assigned to all submatrix elements.
// \return Reference to the assigned submatrix.
//
// This function homogeneously assigns the given value to all dense matrix elements. Note that in
// case the underlying dense matrix is a lower/upper matrix only lower/upper and diagonal elements
// of the underlying matrix are modified.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::operator=( const ElementType& rhs )
{
   const size_t jend( column() + columns() );
   decltype(auto) left( derestrict( matrix_ ) );

   for( size_t j=column(); j<jend; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i ) {
         if( !IsRestricted_v<MT> || IsTriangular_v<MT> || trySet( matrix_, i, j, rhs ) )
            left(i,j) = rhs;
      }
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief List assignment to all submatrix elements.
//
// \param list The initializer list.
// \exception std::invalid_argument Invalid assignment to submatrix.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// This assignment operator offers the option to directly assign to all elements of the submatrix
// by means of an initializer list. The submatrix elements are assigned the values from the given
// initializer list. Missing values are initialized as default. Note that in case the size of the
// top-level initializer list does not match the number of rows of the submatrix or the size of
// any nested list exceeds the number of columns, a \a std::invalid_argument exception is thrown.
// Also, if the underlying matrix \a MT is restricted and the assignment would violate an
// invariant of the matrix, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::operator=( initializer_list< initializer_list<ElementType> > list )
{
   using blaze::reset;

   if( list.size() != rows() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to submatrix" );
   }

   if( IsRestricted_v<MT> ) {
      const InitializerMatrix<ElementType> tmp( list, columns() );
      if( !tryAssign( matrix_, tmp, row(), column() ) ) {
         BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
      }
   }

   decltype(auto) left( derestrict( *this ) );
   size_t i( 0UL );

   for( const auto& rowList : list ) {
      size_t j( 0UL );
      for( const auto& element : rowList ) {
         left(i,j) = element;
         ++j;
      }
      for( ; j<columns(); ++j ) {
         reset( left(i,j) );
      }
      ++i;
   }

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Copy assignment operator for Submatrix.
//
// \param rhs Sparse submatrix to be copied.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Submatrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the current
// sizes of the two submatrices don't match, a \a std::invalid_argument exception is thrown. Also,
// if the underlying matrix \a MT is a lower triangular, upper triangular, or symmetric matrix
// and the assignment would violate its lower, upper, or symmetry property, respectively, a
// \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::operator=( const Submatrix& rhs )
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );

   if( this == &rhs || ( &matrix_ == &rhs.matrix_ && row() == rhs.row() && column() == rhs.column() ) )
      return *this;

   if( rows() != rhs.rows() || columns() != rhs.columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Submatrix sizes do not match" );
   }

   if( !tryAssign( matrix_, rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( rhs.canAlias( this ) ) {
      const ResultType tmp( rhs );
      smpAssign( left, tmp );
   }
   else {
      smpAssign( left, rhs );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Assignment operator for different matrices.
//
// \param rhs Matrix to be assigned.
// \return Reference to the assigned submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// The dense submatrix is initialized as a copy of the given dense submatrix. In case the
// current sizes of the two matrices don't match, a \a std::invalid_argument exception is
// thrown. Also, if the underlying matrix \a MT is a symmetric matrix and the assignment
// would violate its symmetry, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::operator=( const Matrix<MT2,SO>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   using Right = If_t< IsRestricted_v<MT>, CompositeType_t<MT2>, const MT2& >;
   Right right( *rhs );

   if( !tryAssign( matrix_, right, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsReference_v<Right> && right.canAlias( this ) ) {
      const ResultType_t<MT2> tmp( right );
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      if( IsSparseMatrix_v<MT2> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( right ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !tryAddAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const AddType tmp( *this + (*rhs) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpAddAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Addition assignment operator for the addition of a matrix (\f$ A+=B \f$).
//
// \param rhs The right-hand side matrix to be added to the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator+=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using AddType = AddTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( AddType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( AddType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const AddType tmp( *this + (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySubAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SubType tmp( *this - (*rhs ) );
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSubAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Subtraction assignment operator for the subtraction of a matrix (\f$ A-=B \f$).
//
// \param rhs The right-hand side matrix to be subtracted from the submatrix.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO >         // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator-=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SubType = SubTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( SubType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SubType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SubType tmp( *this - (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> DisableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   if( !trySchurAssign( matrix_, *rhs, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( ( ( IsSymmetric_v<MT> || IsHermitian_v<MT> ) && hasOverlap() ) || (*rhs).canAlias( this ) ) {
      const SchurType tmp( *this % (*rhs) );
      if( IsSparseMatrix_v<SchurType> )
         reset();
      smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );
   }
   else {
      smpSchurAssign( left, transIf< IsSymmetric_v<This> >( *rhs ) );
   }

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Schur product assignment operator for the multiplication of a matrix (\f$ A=B \f$).
//
// \param rhs The right-hand side matrix for the Schur product.
// \return Reference to the dense submatrix.
// \exception std::invalid_argument Matrix sizes do not match.
// \exception std::invalid_argument Invalid assignment to restricted matrix.
//
// In case the current sizes of the two matrices don't match, a \a std::invalid_argument exception
// is thrown. Also, if the underlying matrix \a MT is a lower triangular, upper triangular, or
// symmetric matrix and the assignment would violate its lower, upper, or symmetry property,
// respectively, a \a std::invalid_argument exception is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2      // Type of the right-hand side matrix
        , bool SO  >        // Storage order of the right-hand side matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::operator%=( const Matrix<MT2,SO>& rhs )
   -> EnableIf_t< EnforceEvaluation_v<MT,MT2>, Submatrix& >
{
   BLAZE_CONSTRAINT_MUST_BE_DENSE_MATRIX_TYPE  ( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<MT2> );

   using SchurType = SchurTrait_t< ResultType, ResultType_t<MT2> >;

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( SchurType );

   if( rows() != (*rhs).rows() || columns() != (*rhs).columns() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Matrix sizes do not match" );
   }

   const SchurType tmp( *this % (*rhs) );

   if( !tryAssign( matrix_, tmp, row(), column() ) ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid assignment to restricted matrix" );
   }

   decltype(auto) left( derestrict( *this ) );

   if( IsSparseMatrix_v<SchurType> ) {
      reset();
   }

   smpAssign( left, transIf< IsSymmetric_v<This> >( tmp ) );

   BLAZE_INTERNAL_ASSERT( isIntact( matrix_ ), "Invariant violation detected" );

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
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline MT& Submatrix<MT,aligned,true,true,CSAs...>::operand() noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the matrix containing the submatrix.
//
// \return The matrix containing the submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline const MT& Submatrix<MT,aligned,true,true,CSAs...>::operand() const noexcept
{
   return matrix_;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the spacing between the beginning of two columns.
//
// \return The spacing between the beginning of two columns.
//
// This function returns the spacing between the beginning of two columns, i.e. the total
// number of elements of a column.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,true,true,CSAs...>::spacing() const noexcept
{
   return matrix_.spacing();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the maximum capacity of the dense submatrix.
//
// \return The capacity of the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,true,true,CSAs...>::capacity() const noexcept
{
   return rows() * columns();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the current capacity of the specified column.
//
// \param j The index of the column.
// \return The current capacity of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,true,true,CSAs...>::capacity( size_t j ) const noexcept
{
   MAYBE_UNUSED( j );

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   return rows();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the dense submatrix
//
// \return The number of non-zero elements in the dense submatrix.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,true,true,CSAs...>::nonZeros() const
{
   const size_t iend( row() + rows() );
   const size_t jend( column() + columns() );
   size_t nonzeros( 0UL );

   for( size_t j=column(); j<jend; ++j )
      for( size_t i=row(); i<iend; ++i )
         if( !isDefault( matrix_(i,j) ) )
            ++nonzeros;

   return nonzeros;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns the number of non-zero elements in the specified column.
//
// \param j The index of the column.
// \return The number of non-zero elements of column \a j.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline size_t Submatrix<MT,aligned,true,true,CSAs...>::nonZeros( size_t j ) const
{
   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t iend( row() + rows() );
   size_t nonzeros( 0UL );

   for( size_t i=row(); i<iend; ++i )
      if( !isDefault( matrix_(i,column()+j) ) )
         ++nonzeros;

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
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,aligned,true,true,CSAs...>::reset()
{
   using blaze::clear;

   for( size_t j=column(); j<column()+columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i )
         clear( matrix_(i,j) );
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Reset the specified column to the default initial values.
//
// \param j The index of the column.
// \return void
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline void Submatrix<MT,aligned,true,true,CSAs...>::reset( size_t j )
{
   using blaze::clear;

   BLAZE_USER_ASSERT( j < columns(), "Invalid column access index" );

   const size_t ibegin( ( IsLower_v<MT> )
                        ?( ( IsUniLower_v<MT> || IsStrictlyLower_v<MT> )
                           ?( max( j+1UL, row() ) )
                           :( max( j, row() ) ) )
                        :( row() ) );
   const size_t iend  ( ( IsUpper_v<MT> )
                        ?( ( IsUniUpper_v<MT> || IsStrictlyUpper_v<MT> )
                           ?( min( j, row()+rows() ) )
                           :( min( j+1UL, row()+rows() ) ) )
                        :( row()+rows() ) );

   for( size_t i=ibegin; i<iend; ++i )
      clear( matrix_(i,column()+j) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Checking whether there exists an overlap in the context of a symmetric matrix.
//
// \return \a true in case an overlap exists, \a false if not.
//
// This function checks if in the context of a symmetric matrix the submatrix has an overlap with
// its counterpart. In case an overlap exists, the function return \a true, otherwise it returns
// \a false.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,true,true,CSAs...>::hasOverlap() const noexcept
{
   BLAZE_INTERNAL_ASSERT( IsSymmetric_v<MT> || IsHermitian_v<MT>, "Invalid matrix detected" );

   if( ( row() + rows() <= column() ) || ( column() + columns() <= row() ) )
      return false;
   else return true;
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
/*!\brief In-place transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::transpose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, trans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( trans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief In-place conjugate transpose of the submatrix.
//
// \return Reference to the transposed submatrix.
// \exception std::logic_error Invalid transpose of a non-quadratic submatrix.
// \exception std::logic_error Invalid transpose operation.
//
// This function transposes the dense submatrix in-place. Note that this function can only be used
// for quadratic submatrices, i.e. if the number of rows is equal to the number of columns. Also,
// the function fails if ...
//
//  - ... the submatrix contains elements from the upper part of the underlying lower matrix;
//  - ... the submatrix contains elements from the lower part of the underlying upper matrix;
//  - ... the result would be non-deterministic in case of a symmetric or Hermitian matrix.
//
// In all cases, a \a std::logic_error is thrown.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::ctranspose()
{
   if( rows() != columns() ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose of a non-quadratic submatrix" );
   }

   if( !tryAssign( matrix_, ctrans( *this ), row(), column() ) ) {
      BLAZE_THROW_LOGIC_ERROR( "Invalid transpose operation" );
   }

   decltype(auto) left( derestrict( *this ) );
   const ResultType tmp( ctrans( *this ) );

   smpAssign( left, tmp );

   return *this;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Scaling of the dense submatrix by the scalar value \a scalar (\f$ A=B*s \f$).
//
// \param scalar The scalar value for the submatrix scaling.
// \return Reference to the dense submatrix.
//
// This function scales the submatrix by applying the given scalar value \a scalar to each
// element of the submatrix. For built-in and \c complex data types it has the same effect
// as using the multiplication assignment operator. Note that the function cannot be used
// to scale a submatrix on a lower or upper unitriangular matrix. The attempt to scale
// such a submatrix results in a compile time error!
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the scalar value
inline Submatrix<MT,aligned,true,true,CSAs...>&
   Submatrix<MT,aligned,true,true,CSAs...>::scale( const Other& scalar )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_UNITRIANGULAR_MATRIX_TYPE( MT );

   const size_t jend( column() + columns() );

   for( size_t j=column(); j<jend; ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( ( IsStrictlyLower_v<MT> )
                              ?( max( j+1UL, row() ) )
                              :( max( j, row() ) ) )
                           :( row() ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( ( IsStrictlyUpper_v<MT> )
                              ?( min( j, row()+rows() ) )
                              :( min( j+1UL, row()+rows() ) ) )
                           :( row()+rows() ) );

      for( size_t i=ibegin; i<iend; ++i )
         matrix_(i,j) *= scalar;
   }

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
/*!\brief Returns whether the submatrix can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,aligned,true,true,CSAs...>::canAlias( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can alias with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address can alias with the submatrix. In contrast
// to the isAliased() function this function is allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,aligned,true,true,CSAs...>::canAlias( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename Other >  // Data type of the foreign expression
inline bool Submatrix<MT,aligned,true,true,CSAs...>::isAliased( const Other* alias ) const noexcept
{
   return matrix_.isAliased( &unview( *alias ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is aliased with the given dense submatrix \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this submatrix, \a false if not.
//
// This function returns whether the given address is aliased with the submatrix. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions to
// optimize the evaluation.
*/
template< typename MT        // Type of the dense matrix
        , size_t... CSAs >   // Compile time submatrix arguments
template< typename MT2       // Data type of the foreign dense submatrix
        , AlignmentFlag AF2  // Alignment flag of the foreign dense submatrix
        , bool SO2           // Storage order of the foreign dense submatrix
        , size_t... CSAs2 >  // Compile time submatrix arguments of the foreign dense submatrix
inline bool
   Submatrix<MT,aligned,true,true,CSAs...>::isAliased( const Submatrix<MT2,AF2,SO2,true,CSAs2...>* alias ) const noexcept
{
   return ( matrix_.isAliased( &alias->matrix_ ) &&
            ( row() + rows() > alias->row() ) &&
            ( row() < alias->row() + alias->rows() ) &&
            ( column() + columns() > alias->column() ) &&
            ( column() < alias->column() + alias->columns() ) );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix is properly aligned in memory.
//
// \return \a true in case the submatrix is aligned, \a false if not.
//
// This function returns whether the submatrix is guaranteed to be properly aligned in memory,
// i.e. whether the beginning and the end of each column of the submatrix are guaranteed to
// conform to the alignment restrictions of the underlying element type.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,true,true,CSAs...>::isAligned() const noexcept
{
   return true;
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Returns whether the submatrix can be used in SMP assignments.
//
// \return \a true in case the submatrix can be used in SMP assignments, \a false if not.
//
// This function returns whether the submatrix can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current number of
// rows and/or columns of the submatrix).
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
inline bool Submatrix<MT,aligned,true,true,CSAs...>::canSMPAssign() const noexcept
{
   return ( rows() * columns() >= SMP_DMATASSIGN_THRESHOLD );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs a load of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,true,true,CSAs...>::load( size_t i, size_t j ) const noexcept
{
   return loada( i, j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an aligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number
// of values inside the SIMD element. This function must \b NOT be called explicitly! It is
// used internally for the performance optimized evaluation of expression templates. Calling
// this function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,true,true,CSAs...>::loada( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.loada( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned load of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \return The loaded SIMD element.
//
// This function performs an unaligned load of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE typename Submatrix<MT,aligned,true,true,CSAs...>::SIMDType
   Submatrix<MT,aligned,true,true,CSAs...>::loadu( size_t i, size_t j ) const noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   return matrix_.loadu( row()+i, column()+j );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs a store of a specific SIMD element of the dense submatrix. The row
// index must be smaller than the number of rows and the column index must be smaller than
// the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,true,true,CSAs...>::store( size_t i, size_t j, const SIMDType& value ) noexcept
{
   storea( i, j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,true,true,CSAs...>::storea( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.storea( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Unaligned store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an unaligned store of a specific SIMD element of the dense submatrix.
// The row index must be smaller than the number of rows and the column index must be smaller
// than the number of columns. Additionally, the row index must be a multiple of the number of
// values inside the SIMD element. This function must \b NOT be called explicitly! It is used
// internally for the performance optimized evaluation of expression templates. Calling this
// function explicitly might result in erroneous results and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,true,true,CSAs...>::storeu( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.storeu( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Aligned, non-temporal store of a SIMD element of the submatrix.
//
// \param i Access index for the row. The index has to be in the range [0..M-1].
// \param j Access index for the column. The index has to be in the range [0..N-1].
// \param value The SIMD element to be stored.
// \return void
//
// This function performs an aligned, non-temporal store of a specific SIMD element of the
// dense submatrix. The row index must be smaller than the number of rows and the column
// index must be smaller than the number of columns. Additionally, the row index must be
// a multiple of the number of values inside the SIMD element. This function must \b NOT
// be called explicitly! It is used internally for the performance optimized evaluation of
// expression templates. Calling this function explicitly might result in erroneous results
// and/or in compilation errors.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
BLAZE_ALWAYS_INLINE void
   Submatrix<MT,aligned,true,true,CSAs...>::stream( size_t i, size_t j, const SIMDType& value ) noexcept
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( i < rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i + SIMDSIZE <= rows(), "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( i % SIMDSIZE == 0UL, "Invalid row access index" );
   BLAZE_INTERNAL_ASSERT( j < columns(), "Invalid column access index" );

   matrix_.stream( row()+i, column()+j, value );
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >   // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row()+i    ,column()+j) = (*rhs)(i    ,j);
         matrix_(row()+i+1UL,column()+j) = (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(row()+ipos,column()+j) = (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   if( useStreaming &&
       rows()*columns() > ( cacheSize / ( sizeof(ElementType) * 3UL ) ) &&
       !(*rhs).isAliased( this ) )
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; i<ipos; i+=SIMDSIZE ) {
            left.stream( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
   else
   {
      for( size_t j=0UL; j<columns(); ++j )
      {
         size_t i( 0UL );
         Iterator left( begin(j) );
         ConstIterator_t<MT2> right( (*rhs).begin(j) );

         for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<ipos; i+=SIMDSIZE ) {
            left.store( right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         }
         for( ; i<rows(); ++i ) {
            *left = *right; ++left; ++right;
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::assign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) = (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::assign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be assigned.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::assign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) = element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+j,column()+j) += (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(row()+i    ,column()+j) += (*rhs)(i    ,j);
            matrix_(row()+i+1UL,column()+j) += (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(row()+ipos,column()+j) += (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the addition assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedAddAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() + right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left += *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::addAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) += (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::addAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the addition assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be added.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::addAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) += element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      if( IsDiagonal_v<MT2> ) {
         matrix_(row()+j,column()+j) -= (*rhs)(j,j);
      }
      else {
         for( size_t i=0UL; i<ipos; i+=2UL ) {
            matrix_(row()+i    ,column()+j) -= (*rhs)(i    ,j);
            matrix_(row()+i+1UL,column()+j) -= (*rhs)(i+1UL,j);
         }
         if( ipos < rows() ) {
            matrix_(row()+ipos,column()+j) -= (*rhs)(ipos,j);
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the subtraction assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSubAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ibegin( ( IsLower_v<MT> )
                           ?( prevMultiple( ( IsStrictlyLower_v<MT> ? j+1UL : j ), SIMDSIZE ) )
                           :( 0UL ) );
      const size_t iend  ( ( IsUpper_v<MT> )
                           ?( IsStrictlyUpper_v<MT> ? j : j+1UL )
                           :( rows() ) );
      BLAZE_INTERNAL_ASSERT( ibegin <= iend, "Invalid loop indices detected" );

      const size_t ipos( prevMultiple( iend, SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= iend, "Invalid end calculation" );

      size_t i( ibegin );
      Iterator left( begin(j) + ibegin );
      ConstIterator_t<MT2> right( (*rhs).begin(j) + ibegin );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() - right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<iend; ++i ) {
         *left -= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::subAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) -= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::subAssign( const SparseMatrix<MT2,true>& rhs )
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element )
         matrix_(row()+element->index(),column()+j) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the subtraction assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix to be subtracted.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::subAssign( const SparseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element )
         matrix_(row()+i,column()+element->index()) -= element->value();
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> DisableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   const size_t ipos( prevMultiple( rows(), 2UL ) );
   BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

   for( size_t j=0UL; j<columns(); ++j ) {
      for( size_t i=0UL; i<ipos; i+=2UL ) {
         matrix_(row()+i    ,column()+j) *= (*rhs)(i    ,j);
         matrix_(row()+i+1UL,column()+j) *= (*rhs)(i+1UL,j);
      }
      if( ipos < rows() ) {
         matrix_(row()+ipos,column()+j) *= (*rhs)(ipos,j);
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief SIMD optimized implementation of the Schur product assignment of a column-major dense
//        matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline auto Submatrix<MT,aligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,true>& rhs )
   -> EnableIf_t< VectorizedSchurAssign_v<MT2> >
{
   BLAZE_CONSTRAINT_MUST_BE_VECTORIZABLE_TYPE( ElementType );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      const size_t ipos( prevMultiple( rows(), SIMDSIZE ) );
      BLAZE_INTERNAL_ASSERT( ipos <= rows(), "Invalid end calculation" );

      size_t i( 0UL );
      Iterator left( begin(j) );
      ConstIterator_t<MT2> right( (*rhs).begin(j) );

      for( ; (i+SIMDSIZE*3UL) < ipos; i+=SIMDSIZE*4UL ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<ipos; i+=SIMDSIZE ) {
         left.store( left.load() * right.load() ); left += SIMDSIZE; right += SIMDSIZE;
      }
      for( ; i<rows(); ++i ) {
         *left *= *right; ++left; ++right;
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major dense matrix.
//
// \param rhs The right-hand side dense matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side dense matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::schurAssign( const DenseMatrix<MT2,false>& rhs )
{
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   constexpr size_t block( BLOCK_SIZE );

   for( size_t jj=0UL; jj<columns(); jj+=block ) {
      const size_t jend( ( columns()<(jj+block) )?( columns() ):( jj+block ) );
      for( size_t ii=0UL; ii<rows(); ii+=block ) {
         const size_t iend( ( rows()<(ii+block) )?( rows() ):( ii+block ) );
         for( size_t j=jj; j<jend; ++j ) {
            for( size_t i=ii; i<iend; ++i ) {
               matrix_(row()+i,column()+j) *= (*rhs)(i,j);
            }
         }
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a column-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::schurAssign( const SparseMatrix<MT2,true>& rhs )
{
   using blaze::reset;

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t j=0UL; j<columns(); ++j )
   {
      size_t i( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(j); element!=(*rhs).end(j); ++element ) {
         for( ; i<element->index(); ++i )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++i;
      }

      for( ; i<rows(); ++i ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************


//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
/*!\brief Default implementation of the Schur product assignment of a row-major sparse matrix.
//
// \param rhs The right-hand side sparse matrix for the Schur product.
// \return void
//
// This function must \b NOT be called explicitly! It is used internally for the performance
// optimized evaluation of expression templates. Calling this function explicitly might result
// in erroneous results and/or in compilation errors. Instead of using this function use the
// assignment operator.
*/
template< typename MT       // Type of the dense matrix
        , size_t... CSAs >  // Compile time submatrix arguments
template< typename MT2 >    // Type of the right-hand side sparse matrix
inline void Submatrix<MT,aligned,true,true,CSAs...>::schurAssign( const SparseMatrix<MT2,false>& rhs )
{
   using blaze::reset;

   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT2 );

   BLAZE_INTERNAL_ASSERT( rows()    == (*rhs).rows()   , "Invalid number of rows"    );
   BLAZE_INTERNAL_ASSERT( columns() == (*rhs).columns(), "Invalid number of columns" );

   for( size_t i=0UL; i<rows(); ++i )
   {
      size_t j( 0UL );

      for( ConstIterator_t<MT2> element=(*rhs).begin(i); element!=(*rhs).end(i); ++element ) {
         for( ; j<element->index(); ++j )
            reset( matrix_(row()+i,column()+j) );
         matrix_(row()+i,column()+j) *= element->value();
         ++j;
      }

      for( ; j<columns(); ++j ) {
         reset( matrix_(row()+i,column()+j) );
      }
   }
}
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
