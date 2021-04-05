//=================================================================================================
/*!
//  \file blaze/math/sparse/CompressedVector.h
//  \brief Implementation of an arbitrarily sized compressed vector
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

#ifndef _BLAZE_MATH_SPARSE_COMPRESSEDVECTOR_H_
#define _BLAZE_MATH_SPARSE_COMPRESSEDVECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <algorithm>
#include <utility>
#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/DenseVector.h>
#include <blaze/math/constraints/RequiresEvaluation.h>
#include <blaze/math/constraints/SameTag.h>
#include <blaze/math/constraints/TransposeFlag.h>
#include <blaze/math/Exception.h>
#include <blaze/math/expressions/DenseVector.h>
#include <blaze/math/expressions/SparseVector.h>
#include <blaze/math/Forward.h>
#include <blaze/math/InitializerList.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/Serial.h>
#include <blaze/math/sparse/Forward.h>
#include <blaze/math/sparse/ValueIndexPair.h>
#include <blaze/math/sparse/VectorAccessProxy.h>
#include <blaze/math/traits/AddTrait.h>
#include <blaze/math/traits/BandTrait.h>
#include <blaze/math/traits/ColumnTrait.h>
#include <blaze/math/traits/CrossTrait.h>
#include <blaze/math/traits/DivTrait.h>
#include <blaze/math/traits/ElementsTrait.h>
#include <blaze/math/traits/EvaluateTrait.h>
#include <blaze/math/traits/KronTrait.h>
#include <blaze/math/traits/MapTrait.h>
#include <blaze/math/traits/MultTrait.h>
#include <blaze/math/traits/RepeatTrait.h>
#include <blaze/math/traits/RowTrait.h>
#include <blaze/math/traits/SubTrait.h>
#include <blaze/math/traits/SubvectorTrait.h>
#include <blaze/math/typetraits/HighType.h>
#include <blaze/math/typetraits/IsColumnVector.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsResizable.h>
#include <blaze/math/typetraits/IsRowVector.h>
#include <blaze/math/typetraits/IsScalar.h>
#include <blaze/math/typetraits/IsShrinkable.h>
#include <blaze/math/typetraits/IsSMPAssignable.h>
#include <blaze/math/typetraits/IsSparseMatrix.h>
#include <blaze/math/typetraits/IsSparseVector.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/math/typetraits/IsZero.h>
#include <blaze/math/typetraits/LowType.h>
#include <blaze/math/typetraits/TransposeFlag.h>
#include <blaze/system/Thresholds.h>
#include <blaze/system/TransposeFlag.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/algorithms/Transfer.h>
#include <blaze/util/Assert.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/SameSize.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/IntegralConstant.h>
#include <blaze/util/Memory.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsFloatingPoint.h>
#include <blaze/util/typetraits/IsIntegral.h>
#include <blaze/util/typetraits/RemoveConst.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup compressed_vector CompressedVector
// \ingroup sparse_vector
*/
/*!\brief Efficient implementation of an arbitrary sized sparse vector.
// \ingroup compressed_vector
//
// The CompressedVector class is the representation of an arbitrarily sized sparse vector,
// which stores only non-zero elements of arbitrary type. The type of the elements, the transpose
// flag, and the group tag of the vector can be specified via the three template parameters:

   \code
   namespace blaze {

   template< typename Type, bool TF, typename Tag >
   class CompressedVector;

   } // namespace blaze
   \endcode

//  - Type: specifies the type of the vector elements. CompressedVector can be used with any
//          non-cv-qualified, non-reference, non-pointer element type.
//  - TF  : specifies whether the vector is a row vector (\a blaze::rowVector) or a column
//          vector (\a blaze::columnVector). The default value is \a blaze::defaultTransposeFlag.
//  - Tag : optional type parameter to tag the vector. The default type is \a blaze::Group0.
//          See \ref grouping_tagging for details.
//
// Inserting/accessing elements in a compressed vector can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   // Creating a compressed column vector of size 100
   CompressedVector<double,columnVector> a( 100 );

   // The subscript operator provides access to all possible elements of the compressed vector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse vector, the element is inserted into the vector.
   a[42] = 2.0;

   // The second operation for inserting elements is the set() function. In case the element
   // is not contained in the vector it is inserted into the vector, if it is already contained
   // in the vector its value is modified.
   a.set( 45, -1.2 );

   // An alternative for inserting elements into the vector is the insert() function. However,
   // it inserts the element only in case the element is not already contained in the vector.
   a.insert( 50, 3.7 );

   // A very efficient way to add new elements to a sparse vector is the append() function.
   // Note that append() requires that the appended element's index is strictly larger than
   // the currently largest non-zero index of the vector and that the vector's capacity
   // is large enough to hold the new element.
   a.reserve( 10 );
   a.append( 51, -2.1 );

   // In order to traverse all non-zero elements currently stored in the vector, the begin()
   // and end() functions can be used. In the example, all non-zero elements of the vector are
   // traversed.
   for( CompressedVector<double,false>::Iterator it=a.begin(); it!=a.end(); ++it ) {
      ... = it->value();  // Access to the value of the non-zero element
      ... = it->index();  // Access to the index of the non-zero element
   }
   \endcode

// The use of CompressedVector is very natural and intuitive. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// vectors with fitting element types. The following example gives an impression of the use of
// CompressedVector:

   \code
   using blaze::CompressedVector;
   using blaze::DynamicVector;
   using blaze::CompressedMatrix;

   CompressedVector<double> a( 2 );  // Default constructed, non-initialized 2D vectors
   a[0] = 1.0;                       // Initialization of the first element
   a[1] = 2.0;                       // Initialization of the second element

   CompressedVector<double> b( 2 );        // Empty sparse vector
   DynamicVector<float>     c( 2, 2.0F );  // Directly, homogeneously initialized dense vector
   CompressedVector<double> d;             // Default constructed dynamic vector
   CompressedMatrix<double> A;             // Default constructed row-major matrix

   d = a + b;  // Vector addition between vectors of equal element type
   d = a - c;  // Vector subtraction between a sparse and dense vector with different element types
   d = a * b;  // Component-wise vector multiplication

   a *= 2.0;      // In-place scaling of vector
   d  = a * 2.0;  // Scaling of vector a
   d  = 2.0 * a;  // Scaling of vector a

   d += a - b;  // Addition assignment
   d -= a + c;  // Subtraction assignment
   d *= a * b;  // Multiplication assignment

   double scalar = trans( a ) * b;  // Scalar/dot/inner product between two vectors

   A = a * trans( b );  // Outer product between two vectors
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
class CompressedVector
   : public SparseVector< CompressedVector<Type,TF,Tag>, TF >
{
 private:
   //**Type definitions****************************************************************************
   using ElementBase  = ValueIndexPair<Type>;  //!< Base class for the compressed vector element.
   using IteratorBase = ElementBase*;          //!< Iterator over non-constant base elements.
   //**********************************************************************************************

   //**Private class Element***********************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Value-index-pair for the CompressedVector class.
   //
   // This struct grants access to the data members of the base class and adapts the copy and
   // move semantics of the value-index-pair.
   */
   struct Element
      : public ElementBase
   {
      //**Constructors*****************************************************************************
      Element() = default;
      Element( const Element& rhs ) = default;
      Element( Element&& rhs ) = default;
      //*******************************************************************************************

      //**Assignment operators*********************************************************************
      inline Element& operator=( const Element& rhs )
      {
         this->value_ = rhs.value_;
         return *this;
      }

      inline Element& operator=( Element&& rhs )
      {
         this->value_ = std::move( rhs.value_ );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& rhs )
         -> EnableIf_t< IsSparseElement_v<Other>, Element& >
      {
         this->value_ = rhs.value();
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& rhs )
         -> EnableIf_t< IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( rhs.value() );
         return *this;
      }

      template< typename Other >
      inline auto operator=( const Other& v )
         -> EnableIf_t< !IsSparseElement_v<Other>, Element& >
      {
         this->value_ = v;
         return *this;
      }

      template< typename Other >
      inline auto operator=( Other&& v )
         -> EnableIf_t< !IsSparseElement_v< RemoveReference_t<Other> > &&
                        IsRValueReference_v<Other&&>, Element& >
      {
         this->value_ = std::move( v );
         return *this;
      }
      //*******************************************************************************************

      //**Friend declarations**********************************************************************
      friend class CompressedVector;
      //*******************************************************************************************
   };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using This       = CompressedVector<Type,TF,Tag>;  //!< Type of this CompressedVector instance.
   using BaseType   = SparseVector<This,TF>;          //!< Base type of this CompressedVector instance.
   using ResultType = This;                           //!< Result type for expression template evaluations.

   //! Transpose type for expression template evaluations.
   using TransposeType = CompressedVector<Type,!TF,Tag>;

   using ElementType   = Type;                     //!< Type of the compressed vector elements.
   using TagType       = Tag;                      //!< Tag type of this CompressedVector instance.
   using ReturnType    = const Type&;              //!< Return type for expression template evaluations.
   using CompositeType = const CompressedVector&;  //!< Data type for composite expression templates.

   using Reference      = VectorAccessProxy<This>;  //!< Reference to a non-constant vector value.
   using ConstReference = const Type&;              //!< Reference to a constant vector value.
   using Iterator       = Element*;                 //!< Iterator over non-constant elements.
   using ConstIterator  = const Element*;           //!< Iterator over constant elements.
   //**********************************************************************************************

   //**Rebind struct definition********************************************************************
   /*!\brief Rebind mechanism to obtain a CompressedVector with different data/element type.
   */
   template< typename NewType >  // Data type of the other vector
   struct Rebind {
      using Other = CompressedVector<NewType,TF,Tag>;  //!< The type of the other CompressedVector.
   };
   //**********************************************************************************************

   //**Resize struct definition********************************************************************
   /*!\brief Resize mechanism to obtain a CompressedVector with a different fixed number of elements.
   */
   template< size_t NewN >  // Number of elements of the other vector
   struct Resize {
      using Other = CompressedVector<Type,TF,Tag>;  //!< The type of the other CompressedVector.
   };
   //**********************************************************************************************

   //**Compilation flags***************************************************************************
   //! Compilation flag for SMP assignments.
   /*! The \a smpAssignable compilation flag indicates whether the vector can be used in SMP
       (shared memory parallel) assignments (both on the left-hand and right-hand side of the
       assignment). */
   static constexpr bool smpAssignable = !IsSMPAssignable_v<Type>;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
            inline CompressedVector() noexcept;
   explicit inline CompressedVector( size_t size ) noexcept;
            inline CompressedVector( size_t size, size_t nonzeros );
            inline CompressedVector( initializer_list<Type> list );

   inline CompressedVector( const CompressedVector& sv );
   inline CompressedVector( CompressedVector&& sv ) noexcept;

   template< typename VT > inline CompressedVector( const DenseVector<VT,TF>&  dv );
   template< typename VT > inline CompressedVector( const SparseVector<VT,TF>& sv );
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   inline ~CompressedVector();
   //@}
   //**********************************************************************************************

   //**Data access functions***********************************************************************
   /*!\name Data access functions */
   //@{
   inline Reference      operator[]( size_t index ) noexcept;
   inline ConstReference operator[]( size_t index ) const noexcept;
   inline Reference      at( size_t index );
   inline ConstReference at( size_t index ) const;
   inline Iterator       begin () noexcept;
   inline ConstIterator  begin () const noexcept;
   inline ConstIterator  cbegin() const noexcept;
   inline Iterator       end   () noexcept;
   inline ConstIterator  end   () const noexcept;
   inline ConstIterator  cend  () const noexcept;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   inline CompressedVector& operator=( initializer_list<Type> list ) &;
   inline CompressedVector& operator=( const CompressedVector& rhs ) &;
   inline CompressedVector& operator=( CompressedVector&& rhs ) & noexcept;

   template< typename VT > inline CompressedVector& operator= ( const DenseVector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator= ( const SparseVector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator+=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator-=( const Vector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator*=( const DenseVector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator*=( const SparseVector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator/=( const DenseVector<VT,TF>& rhs ) &;
   template< typename VT > inline CompressedVector& operator%=( const Vector<VT,TF>& rhs ) &;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t size() const noexcept;
   inline size_t capacity() const noexcept;
   inline size_t nonZeros() const;
   inline void   reset();
   inline void   clear();
   inline void   resize( size_t n, bool preserve=true );
          void   reserve( size_t n );
   inline void   shrinkToFit();
   inline void   swap( CompressedVector& sv ) noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions*************************************************************************
   /*!\name Insertion functions */
   //@{
   inline Iterator set   ( size_t index, const Type& value );
   inline Iterator insert( size_t index, const Type& value );
   inline void     append( size_t index, const Type& value, bool check=false );
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
   template< typename Other > inline CompressedVector& scale( const Other& scalar );
   //@}
   //**********************************************************************************************

   //**Expression template evaluation functions****************************************************
   /*!\name Expression template evaluation functions */
   //@{
   template< typename Other > inline bool canAlias ( const Other* alias ) const noexcept;
   template< typename Other > inline bool isAliased( const Other* alias ) const noexcept;

   inline bool canSMPAssign() const noexcept;

   template< typename VT > inline void assign    ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void assign    ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void addAssign ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void addAssign ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void subAssign ( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void subAssign ( const SparseVector<VT,TF>& rhs );
   template< typename VT > inline void multAssign( const DenseVector <VT,TF>& rhs );
   template< typename VT > inline void divAssign ( const DenseVector <VT,TF>& rhs );
   //@}
   //**********************************************************************************************

 private:
   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline size_t       extendCapacity() const noexcept;
   inline Iterator     castDown( IteratorBase it ) const noexcept;
   inline IteratorBase castUp  ( Iterator     it ) const noexcept;
   //@}
   //**********************************************************************************************

   //**Insertion functions***************************************************************************
   /*!\name Insertion functions */
   //@{
   Iterator insert( Iterator pos, size_t index, const Type& value );
   //@}
   //**********************************************************************************************

   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   size_t size_;             //!< The current size/dimension of the compressed vector.
   size_t capacity_;         //!< The maximum capacity of the compressed vector.
   Iterator begin_;          //!< Pointer to the first non-zero element of the compressed vector.
   Iterator end_;            //!< Pointer one past the last non-zero element of the compressed vector.

   static const Type zero_;  //!< Neutral element for accesses to zero elements.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   BLAZE_CONSTRAINT_MUST_HAVE_SAME_SIZE( ElementBase, Element );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  DEFINITION AND INITIALIZATION OF THE STATIC MEMBER VARIABLES
//
//=================================================================================================

template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
const Type CompressedVector<Type,TF,Tag>::zero_{};




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The default constructor for CompressedVector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector() noexcept
   : size_    ( 0UL )      // The current size/dimension of the compressed vector
   , capacity_( 0UL )      // The maximum capacity of the compressed vector
   , begin_   ( nullptr )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( nullptr )  // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a compressed vector of size \a n.
//
// \param n The size of the vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector( size_t n ) noexcept
   : size_    ( n   )      // The current size/dimension of the compressed vector
   , capacity_( 0UL )      // The maximum capacity of the compressed vector
   , begin_   ( nullptr )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( nullptr )  // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a compressed vector of size \a n.
//
// \param n The size of the vector.
// \param nonzeros The number of expected non-zero elements.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector( size_t n, size_t nonzeros )
   : size_    ( n )                               // The current size/dimension of the compressed vector
   , capacity_( nonzeros )                        // The maximum capacity of the compressed vector
   , begin_   ( allocate<Element>( capacity_ ) )  // Pointer to the first non-zero element of the compressed vector
   , end_     ( begin_ )                          // Pointer to the last non-zero element of the compressed vector
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief List initialization of all vector elements.
//
// \param list The initializer list.
//
// This assignment operator provides the option to explicitly initialize the elements of the
// vector within a constructor call:

   \code
   blaze::CompressedVector<double> v1{ 4.2, 6.3, -1.2 };
   \endcode

// The vector is sized according to the size of the initializer list and all its elements are
// initialized by the non-zero elements of the given initializer list.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector( initializer_list<Type> list )
   : CompressedVector( list.size(), blaze::nonZeros( list ) )
{
   size_t i( 0UL );

   for( const Type& element : list ) {
      if( !isDefault<strict>( element ) )
         append( i, element );
      ++i;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The copy constructor for CompressedVector.
//
// \param sv Compressed vector to be copied.
//
// The copy constructor is explicitly defined due to the required dynamic memory management
// and in order to enable/facilitate NRV optimization.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector( const CompressedVector& sv )
   : CompressedVector( sv.size_, sv.nonZeros() )
{
   end_ = begin_ + capacity_;
   std::copy( sv.begin_, sv.end_, castUp( begin_ ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief The move constructor for CompressedVector.
//
// \param sv The compressed vector to be moved into this instance.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::CompressedVector( CompressedVector&& sv ) noexcept
   : size_    ( sv.size_ )      // The current size/dimension of the compressed vector
   , capacity_( sv.capacity_ )  // The maximum capacity of the compressed vector
   , begin_   ( sv.begin_ )     // Pointer to the first non-zero element of the compressed vector
   , end_     ( sv.end_ )       // Pointer to the last non-zero element of the compressed vector
{
   sv.size_     = 0UL;
   sv.capacity_ = 0UL;
   sv.begin_    = nullptr;
   sv.end_      = nullptr;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from dense vectors.
//
// \param dv Dense vector to be copied.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the foreign dense vector
inline CompressedVector<Type,TF,Tag>::CompressedVector( const DenseVector<VT,TF>& dv )
   : CompressedVector( (*dv).size() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   assign( *this, *dv );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Conversion constructor from different sparse vectors.
//
// \param sv Sparse vector to be copied.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the foreign sparse vector
inline CompressedVector<Type,TF,Tag>::CompressedVector( const SparseVector<VT,TF>& sv )
   : CompressedVector( (*sv).size(), (*sv).nonZeros() )
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   assign( *this, *sv );
}
//*************************************************************************************************




//=================================================================================================
//
//  DESTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The destructor for CompressedVector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>::~CompressedVector()
{
   deallocate( begin_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
//
// This function returns a reference to the accessed value at position \a index. In case the
// compressed vector does not yet store an element for index \a index, a new element is inserted
// into the compressed vector. An alternative for traversing the non-zero elements of the sparse
// vector are the begin() and end() functions.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Reference
   CompressedVector<Type,TF,Tag>::operator[]( size_t index ) noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   return Reference( *this, index );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subscript operator for the direct access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstReference
   CompressedVector<Type,TF,Tag>::operator[]( size_t index ) const noexcept
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const ConstIterator pos( lowerBound( index ) );

   if( pos == end_ || pos->index_ != index )
      return zero_;
   else
      return pos->value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid compressed vector access index.
//
// This function returns a reference to the accessed value at position \a index. In case the
// compressed vector does not yet store an element for index \a index, a new element is inserted
// into the compressed vector. In contrast to the subscript operator this function always
// performs a check of the given access index.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Reference
   CompressedVector<Type,TF,Tag>::at( size_t index )
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid compressed vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checked access to the compressed vector elements.
//
// \param index Access index. The index has to be in the range \f$[0..N-1]\f$.
// \return Reference to the accessed value.
// \exception std::out_of_range Invalid compressed vector access index.
//
// In contrast to the subscript operator this function always performs a check of the given
// access indices.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstReference
   CompressedVector<Type,TF,Tag>::at( size_t index ) const
{
   if( index >= size_ ) {
      BLAZE_THROW_OUT_OF_RANGE( "Invalid compressed vector access index" );
   }
   return (*this)[index];
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the compressed vector.
//
// \return Iterator to the first non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::begin() noexcept
{
   return Iterator( begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the compressed vector.
//
// \return Iterator to the first non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::begin() const noexcept
{
   return ConstIterator( begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first non-zero element of the compressed vector.
//
// \return Iterator to the first non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::cbegin() const noexcept
{
   return ConstIterator( begin_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the compressed vector.
//
// \return Iterator just past the last non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::end() noexcept
{
   return Iterator( end_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the compressed vector.
//
// \return Iterator just past the last non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::end() const noexcept
{
   return ConstIterator( end_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator just past the last non-zero element of the compressed vector.
//
// \return Iterator just past the last non-zero element of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::cend() const noexcept
{
   return ConstIterator( end_ );
}
//*************************************************************************************************




//=================================================================================================
//
//  ASSIGNMENT OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief List assignment to all vector elements.
//
// \param list The initializer list.
//
// This assignment operator offers the option to directly assign to all elements of the vector
// by means of an initializer list:

   \code
   blaze::CompressedVector<double> v;
   v = { 4.2, 6.3, -1.2 };
   \endcode

// The vector is resized according to the size of the initializer list and all its elements are
// assigned the non-zero values from the given initializer list.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator=( initializer_list<Type> list ) &
{
   using blaze::nonZeros;

   resize( list.size(), false );
   reserve( nonZeros( list ) );

   size_t i( 0UL );

   for( const Type& element : list ) {
      if( !isDefault<strict>( element ) )
         append( i, element );
      ++i;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Copy assignment operator for CompressedVector.
//
// \param rhs Compressed vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The compressed vector is resized according to the given compressed vector and initialized
// as a copy of this vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator=( const CompressedVector& rhs ) &
{
   using std::swap;

   if( &rhs == this ) return *this;

   const size_t nonzeros( rhs.nonZeros() );

   if( nonzeros > capacity_ ) {
      Iterator newBegin( allocate<Element>( nonzeros ) );
      end_ = castDown( std::copy( rhs.begin_, rhs.end_, castUp( newBegin ) ) );
      swap( begin_, newBegin );
      deallocate( newBegin );

      size_     = rhs.size_;
      capacity_ = nonzeros;
   }
   else {
      end_  = castDown( std::copy( rhs.begin_, rhs.end_, castUp( begin_ ) ) );
      size_ = rhs.size_;
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Move assignment operator for CompressedVector.
//
// \param rhs The compressed vector to be moved into this instance.
// \return Reference to the assigned compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator=( CompressedVector&& rhs ) & noexcept
{
   deallocate( begin_ );

   size_     = rhs.size_;
   capacity_ = rhs.capacity_;
   begin_    = rhs.begin_;
   end_      = rhs.end_;

   rhs.size_     = 0UL;
   rhs.capacity_ = 0UL;
   rhs.begin_    = nullptr;
   rhs.end_      = nullptr;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for dense vectors.
//
// \param rhs Dense vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The vector is resized according to the given dense vector and initialized as a copy of
// this vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator=( const DenseVector<VT,TF>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).canAlias( this ) ) {
      CompressedVector tmp( *rhs );
      swap( tmp );
   }
   else {
      size_ = (*rhs).size();
      end_  = begin_;
      assign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different sparse vectors.
//
// \param rhs Sparse vector to be copied.
// \return Reference to the assigned compressed vector.
//
// The vector is resized according to the given sparse vector and initialized as a copy of
// this vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side sparse vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator=( const SparseVector<VT,TF>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).canAlias( this ) || (*rhs).nonZeros() > capacity_ ) {
      CompressedVector tmp( *rhs );
      swap( tmp );
   }
   else {
      size_ = (*rhs).size();
      end_  = begin_;

      if( !IsZero_v<VT> ) {
         assign( *this, *rhs );
      }
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment operator for the addition of a vector (\f$ \vec{a}+=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be added to the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator+=( const Vector<VT,TF>& rhs ) &
{
   using blaze::addAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsZero_v<VT> ) {
      addAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment operator for the subtraction of a vector (\f$ \vec{a}-=\vec{b} \f$).
//
// \param rhs The right-hand side vector to be subtracted from the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator-=( const Vector<VT,TF>& rhs ) &
{
   using blaze::subAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsZero_v<VT> ) {
      subAssign( *this, *rhs );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a dense vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector to be multiplied with the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator*=( const DenseVector<VT,TF>& rhs ) &
{
   using blaze::multAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      CompressedVector tmp( *this * (*rhs) );
      swap( tmp );
   }
   else {
      CompositeType_t<VT> tmp( *rhs );
      multAssign( *this, tmp );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment operator for the multiplication of a sparse vector
//        (\f$ \vec{a}*=\vec{b} \f$).
//
// \param rhs The right-hand side sparse vector to be multiplied with the compressed vector.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator*=( const SparseVector<VT,TF>& rhs ) &
{
   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( !IsZero_v<VT> ) {
      CompressedVector tmp( *this * (*rhs) );
      swap( tmp );
   }
   else {
      reset();
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment operator for the division of a dense vector (\f$ \vec{a}/=\vec{b} \f$).
//
// \param rhs The right-hand side dense vector divisor.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Vector sizes do not match.
//
// In case the current sizes of the two vectors don't match, a \a std::invalid_argument exception
// is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator/=( const DenseVector<VT,TF>& rhs ) &
{
   using blaze::divAssign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   if( (*rhs).size() != size_ ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Vector sizes do not match" );
   }

   if( (*rhs).canAlias( this ) ) {
      CompressedVector tmp( *this / (*rhs) );
      swap( tmp );
   }
   else {
      CompositeType_t<VT> tmp( *rhs );
      divAssign( *this, tmp );
   }

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Cross product assignment operator for the multiplication of a vector
//        (\f$ \vec{a}\times=\vec{b} \f$).
//
// \param rhs The right-hand side vector for the cross product.
// \return Reference to the compressed vector.
// \exception std::invalid_argument Invalid vector size for cross product.
//
// In case the current size of any of the two vectors is not equal to 3, a \a std::invalid_argument
// exception is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side vector
inline CompressedVector<Type,TF,Tag>&
   CompressedVector<Type,TF,Tag>::operator%=( const Vector<VT,TF>& rhs ) &
{
   using blaze::assign;

   BLAZE_CONSTRAINT_MUST_BE_SAME_TAG( Tag, TagType_t<VT> );

   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( ResultType_t<VT>, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( ResultType_t<VT> );

   using CrossType = CrossTrait_t< This, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( CrossType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( CrossType, TF );
   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( CrossType );

   if( size_ != 3UL || (*rhs).size() != 3UL ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid vector size for cross product" );
   }

   if( !IsZero_v<VT> ) {
      const CrossType tmp( *this % (*rhs) );
      reset();
      assign( *this, tmp );
   }
   else {
      reset();
   }

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current size/dimension of the compressed vector.
//
// \return The size of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline size_t CompressedVector<Type,TF,Tag>::size() const noexcept
{
   return size_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the maximum capacity of the compressed vector.
//
// \return The capacity of the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline size_t CompressedVector<Type,TF,Tag>::capacity() const noexcept
{
   return capacity_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the number of non-zero elements in the compressed vector.
//
// \return The number of non-zero elements in the compressed vector.
//
// Note that the number of non-zero elements is always smaller than the current size of the
// compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline size_t CompressedVector<Type,TF,Tag>::nonZeros() const
{
   return end_ - begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reset to the default initial values.
//
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::reset()
{
   end_ = begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the compressed vector.
//
// \return void
//
// After the clear() function, the size of the compressed vector is 0.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::clear()
{
   size_ = 0UL;
   end_  = begin_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Changing the size of the compressed vector.
//
// \param n The new size of the compressed vector.
// \param preserve \a true if the old values of the vector should be preserved, \a false if not.
// \return void
//
// This function resizes the compressed vector using the given size to \a n. During this
// operation, new dynamic memory may be allocated in case the capacity of the compressed
// vector is too small. Note that this function may invalidate all existing views (subvectors,
// ...) on the vector if it is used to shrink the vector. Additionally, the resize operation
// potentially changes all vector elements. In order to preserve the old vector values, the
// \a preserve flag can be set to \a true.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::resize( size_t n, bool preserve )
{
   if( preserve ) {
      end_ = lowerBound( n );
   }
   else {
      end_ = begin_;
   }

   size_ = n;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the minimum capacity of the compressed vector.
//
// \param n The new minimum capacity of the compressed vector.
// \return void
//
// This function increases the capacity of the compressed vector to at least \a n elements. The
// current values of the vector elements are preserved.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
void CompressedVector<Type,TF,Tag>::reserve( size_t n )
{
   using std::swap;

   if( n > capacity_ ) {
      const size_t newCapacity( n );

      // Allocating a new data and index array
      Iterator newBegin  = allocate<Element>( newCapacity );

      // Replacing the old data and index array
      end_ = castDown( transfer( begin_, end_, castUp( newBegin ) ) );
      swap( newBegin, begin_ );
      capacity_ = newCapacity;
      deallocate( newBegin );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Requesting the removal of unused capacity.
//
// \return void
//
// This function minimizes the capacity of the vector by removing unused capacity. Please note
// that in case a reallocation occurs, all iterators (including end() iterators), all pointers
// and references to elements of this vector are invalidated.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::shrinkToFit()
{
   if( nonZeros() < capacity_ ) {
      CompressedVector( *this ).swap( *this );
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two compressed vectors.
//
// \param sv The compressed vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::swap( CompressedVector& sv ) noexcept
{
   using std::swap;

   swap( size_, sv.size_ );
   swap( capacity_, sv.capacity_ );
   swap( begin_, sv.begin_ );
   swap( end_, sv.end_ );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Calculating a new vector capacity.
//
// \return The new compressed vector capacity.
//
// This function calculates a new vector capacity based on the current capacity of the sparse
// vector. Note that the new capacity is restricted to the interval \f$[7..size]\f$.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline size_t CompressedVector<Type,TF,Tag>::extendCapacity() const noexcept
{
   using blaze::max;
   using blaze::min;

   size_t nonzeros( 2UL*capacity_+1UL );
   nonzeros = max( nonzeros, 7UL   );
   nonzeros = min( nonzeros, size_ );

   BLAZE_INTERNAL_ASSERT( nonzeros > capacity_, "Invalid capacity value" );

   return nonzeros;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs a down-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs a down-cast of the given iterator to base elements to an iterator to
// derived elements.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::castDown( IteratorBase it ) const noexcept
{
   return static_cast<Iterator>( it );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Performs an up-cast of the given iterator.
//
// \return The casted iterator.
//
// This function performs an up-cast of the given iterator to derived elements to an iterator
// to base elements.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::IteratorBase
   CompressedVector<Type,TF,Tag>::castUp( Iterator it ) const noexcept
{
   return static_cast<IteratorBase>( it );
}
//*************************************************************************************************




//=================================================================================================
//
//  INSERTION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Setting an element of the compressed vector.
//
// \param index The index of the element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be set.
// \return Reference to the set value.
//
// This function sets the value of an element of the compressed vector. In case the sparse vector
// already contains an element with index \a index its value is modified, else a new element with
// the given \a value is inserted.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::set( size_t index, const Type& value )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const Iterator pos( lowerBound( index ) );

   if( pos != end_ && pos->index_ == index ) {
      pos->value() = value;
      return pos;
   }
   else return insert( pos, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the compressed vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid compressed vector access index.
//
// This function inserts a new element into the compressed vector. However, duplicate elements
// are not allowed. In case the sparse vector already contains an element with index \a index,
// a \a std::invalid_argument exception is thrown.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::insert( size_t index, const Type& value )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const Iterator pos( lowerBound( index ) );

   if( pos != end_ && pos->index_ == index ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Bad access index" );
   }

   return insert( pos, index, value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inserting an element into the compressed vector.
//
// \param pos The position of the new element.
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be inserted.
// \return Reference to the inserted value.
// \exception std::invalid_argument Invalid compressed vector access index.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::insert( Iterator pos, size_t index, const Type& value )
{
   using std::swap;

   if( nonZeros() != capacity_ ) {
      std::move_backward( pos, end_, castUp( end_+1 ) );
      pos->value_ = value;
      pos->index_ = index;
      ++end_;

      return pos;
   }
   else {
      size_t newCapacity( extendCapacity() );

      Iterator newBegin = allocate<Element>( newCapacity );
      Iterator tmp = castDown( std::move( begin_, pos, castUp( newBegin ) ) );
      tmp->value_ = value;
      tmp->index_ = index;
      end_ = castDown( std::move( pos, end_, castUp( tmp+1 ) ) );

      swap( newBegin, begin_ );
      deallocate( newBegin );
      capacity_ = newCapacity;

      return tmp;
   }
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Appending an element to the compressed vector.
//
// \param index The index of the new element. The index has to be in the range \f$[0..N-1]\f$.
// \param value The value of the element to be appended.
// \param check \a true if the new value should be checked for default values, \a false if not.
// \return void
//
// This function provides a very efficient way to fill a compressed vector with elements. It
// appends a new element to the end of the compressed vector without any memory allocation.
// Therefore it is strictly necessary to keep the following preconditions in mind:
//
//  - the index of the new element must be strictly larger than the largest index of non-zero
//    elements in the compressed vector
//  - the current number of non-zero elements must be smaller than the capacity of the vector
//
// Ignoring these preconditions might result in undefined behavior! The optional \a check
// parameter specifies whether the new value should be tested for a default value. If the new
// value is a default value (for instance 0 in case of an integral element type) the value is
// not appended. Per default the values are not tested.
//
// \note Although append() does not allocate new memory, it still invalidates all iterators
// returned by the end() functions!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::append( size_t index, const Type& value, bool check )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );
   BLAZE_USER_ASSERT( nonZeros() < capacity(), "Not enough reserved capacity" );
   BLAZE_USER_ASSERT( begin_ == end_ || (end_-1UL)->index_ < index, "Index is not strictly increasing" );

   end_->value_ = value;

   if( !check || !isDefault<strict>( end_->value_ ) ) {
      end_->index_ = index;
      ++end_;
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  ERASE FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Erasing an element from the compressed vector.
//
// \param index The index of the element to be erased. The index has to be in the range \f$[0..N-1]\f$.
// \return void
//
// This function erases an element from the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void CompressedVector<Type,TF,Tag>::erase( size_t index )
{
   BLAZE_USER_ASSERT( index < size_, "Invalid compressed vector access index" );

   const Iterator pos( find( index ) );
   if( pos != end_ )
      end_ = castDown( std::move( pos+1, end_, castUp( pos ) ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing an element from the compressed vector.
//
// \param pos Iterator to the element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases an element from the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::erase( Iterator pos )
{
   BLAZE_USER_ASSERT( pos >= begin_ && pos <= end_, "Invalid compressed vector iterator" );

   if( pos != end_ )
      end_ = castDown( std::move( pos+1, end_, castUp( pos ) ) );
   return pos;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing a range of elements from the compressed vector.
//
// \param first Iterator to first element to be erased.
// \param last Iterator just past the last element to be erased.
// \return Iterator to the element after the erased element.
//
// This function erases a range of elements from the compressed vector.
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::erase( Iterator first, Iterator last )
{
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range" );
   BLAZE_USER_ASSERT( first >= begin_ && first <= end_, "Invalid compressed vector iterator" );
   BLAZE_USER_ASSERT( last  >= begin_ && last  <= end_, "Invalid compressed vector iterator" );

   if( first != last )
      end_ = castDown( std::move( last, end_, castUp( first ) ) );
   return first;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing specific elements from the compressed vector.
//
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from the compressed vector. The elements are selected
// by the given unary predicate \a predicate, which is expected to accept a single argument of
// the type of the elements and to be pure. The following example demonstrates how to remove all
// elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization

   a.erase( []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename Pred  // Type of the unary predicate
        , typename >     // Type restriction on the unary predicate
inline void CompressedVector<Type,TF,Tag>::erase( Pred predicate )
{
   end_ = castDown( std::remove_if( castUp( begin_ ), castUp( end_ ),
                                    [predicate=predicate]( const ElementBase& element ) {
                                       return predicate( element.value() );
                                    } ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Erasing specific elements from a range of the compressed vector.
//
// \param first Iterator to first element of the range.
// \param last Iterator just past the last element of the range.
// \param predicate The unary predicate for the element selection.
// \return void.
//
// This function erases specific elements from a range of elements of the compressed vector.
// The elements are selected by the given unary predicate \a predicate, which is expected to
// accept a single argument of the type of the elements and to be pure. The following example
// demonstrates how to remove all elements that are smaller than a certain threshold value:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization

   a.erase( a.begin(), a.end(), []( double value ){ return value < 1E-8; } );
   \endcode

// \note The predicate is required to be pure, i.e. to produce deterministic results for elements
// with the same value. The attempt to use an impure predicate leads to undefined behavior!
*/
template< typename Type    // Data type of the vector
        , bool TF          // Transpose flag
        , typename Tag >   // Type tag
template< typename Pred >  // Type of the unary predicate
inline void CompressedVector<Type,TF,Tag>::erase( Iterator first, Iterator last, Pred predicate )
{
   BLAZE_USER_ASSERT( first <= last, "Invalid iterator range" );
   BLAZE_USER_ASSERT( first >= begin_ && first <= end_, "Invalid compressed vector iterator" );
   BLAZE_USER_ASSERT( last  >= begin_ && last  <= end_, "Invalid compressed vector iterator" );

   const auto pos = std::remove_if( castUp( first ), castUp( last  ),
                                    [predicate=predicate]( const ElementBase& element ) {
                                       return predicate( element.value() );
                                    } );

   end_ = castDown( std::move( last, end_, pos ) );
}
//*************************************************************************************************




//=================================================================================================
//
//  LOOKUP FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the compressed vector (the end() iterator) is returned. Note
// that the returned compressed vector iterator is subject to invalidation due to inserting
// operations via the subscript operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::find( size_t index )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).find( index ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Searches for a specific vector element.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the element in case the index is found, end() iterator otherwise.
//
// This function can be used to check whether a specific element is contained in the sparse
// vector. It specifically searches for the element with index \a index. In case the element
// is found, the function returns an iterator to the element. Otherwise an iterator just past
// the last non-zero element of the compressed vector (the end() iterator) is returned. Note
// that the returned compressed vector iterator is subject to invalidation due to inserting
// operations via the subscript operator, the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::find( size_t index ) const
{
   const ConstIterator pos( lowerBound( index ) );
   if( pos != end_ && pos->index_ == index )
      return pos;
   else return end_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed vector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::lowerBound( size_t index )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).lowerBound( index ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index not less then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index not less then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index not less then the given
// index. In combination with the upperBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed vector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::lowerBound( size_t index ) const
{
   return std::lower_bound( begin_, end_, index,
                            []( const Element& element, size_t i )
                            {
                               return element.index() < i;
                            } );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed vector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::Iterator
   CompressedVector<Type,TF,Tag>::upperBound( size_t index )
{
   return const_cast<Iterator>( const_cast<const This&>( *this ).upperBound( index ) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns an iterator to the first index greater then the given index.
//
// \param index The index of the search element. The index has to be in the range \f$[0..N-1]\f$.
// \return Iterator to the first index greater then the given index, end() iterator otherwise.
//
// This function returns an iterator to the first element with an index greater then the given
// index. In combination with the lowerBound() function this function can be used to create a
// pair of iterators specifying a range of indices. Note that the returned compressed vector
// iterator is subject to invalidation due to inserting operations via the subscript operator,
// the set() function or the insert() function!
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline typename CompressedVector<Type,TF,Tag>::ConstIterator
   CompressedVector<Type,TF,Tag>::upperBound( size_t index ) const
{
   return std::upper_bound( begin_, end_, index,
                            []( size_t i, const Element& element )
                            {
                               return i < element.index();
                            } );
}
//*************************************************************************************************




//=================================================================================================
//
//  NUMERIC FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Scaling of the compressed vector by the scalar value \a scalar (\f$ \vec{a}=\vec{b}*s \f$).
//
// \param scalar The scalar value for the vector scaling.
// \return Reference to the compressed vector.
//
// This function scales the vector by applying the given scalar value \a scalar to each element
// of the vector. For built-in and \c complex data types it has the same effect as using the
// multiplication assignment operator:

   \code
   blaze::CompressedVector<int> a;
   // ... Resizing and initialization
   a *= 4;        // Scaling of the vector
   a.scale( 4 );  // Same effect as above
   \endcode
*/
template< typename Type     // Data type of the vector
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the scalar value
inline CompressedVector<Type,TF,Tag>& CompressedVector<Type,TF,Tag>::scale( const Other& scalar )
{
   for( auto element=begin_; element!=end_; ++element )
      element->value_ *= scalar;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  EXPRESSION TEMPLATE EVALUATION FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns whether the vector can alias with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address can alias with the vector. In contrast
// to the isAliased() function this function is allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool CompressedVector<Type,TF,Tag>::canAlias( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector is aliased with the given address \a alias.
//
// \param alias The alias to be checked.
// \return \a true in case the alias corresponds to this vector, \a false if not.
//
// This function returns whether the given address is aliased with the vector. In contrast
// to the canAlias() function this function is not allowed to use compile time expressions
// to optimize the evaluation.
*/
template< typename Type     // Data type of the vector
        , bool TF           // Transpose flag
        , typename Tag >    // Type tag
template< typename Other >  // Data type of the foreign expression
inline bool CompressedVector<Type,TF,Tag>::isAliased( const Other* alias ) const noexcept
{
   return static_cast<const void*>( this ) == static_cast<const void*>( alias );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the vector can be used in SMP assignments.
//
// \return \a true in case the vector can be used in SMP assignments, \a false if not.
//
// This function returns whether the vector can be used in SMP assignments. In contrast to the
// \a smpAssignable member enumeration, which is based solely on compile time information, this
// function additionally provides runtime information (as for instance the current size of the
// vector).
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline bool CompressedVector<Type,TF,Tag>::canSMPAssign() const noexcept
{
   return false;
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF,Tag>::assign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   size_t nonzeros( 0UL );

   for( size_t i=0UL; i<size_; ++i )
   {
      if( nonzeros == capacity_ )
         reserve( extendCapacity() );

      end_->value_ = (*rhs)[i];

      if( !isDefault<strict>( end_->value_ ) ) {
         end_->index_ = i;
         ++end_;
         ++nonzeros;
      }
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF,Tag>::assign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );
   BLAZE_INTERNAL_ASSERT( nonZeros() == 0UL, "Invalid non-zero elements detected" );

   // Using the following formulation instead of a std::copy function call of the form
   //
   //          end_ = std::copy( (*rhs).begin(), (*rhs).end(), begin_ );
   //
   // results in much less requirements on the ConstIterator type provided from the right-hand
   // sparse vector type
   for( auto element=(*rhs).begin(); element!=(*rhs).end(); ++element )
      append( element->index(), element->value() );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF,Tag>::addAssign( const DenseVector<VT,TF>& rhs )
{
   using AddType = AddTrait_t< This, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( AddType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( AddType, TF );
   BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( CompositeType_t<AddType> );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const AddType tmp( serial( *this + (*rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF,Tag>::addAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   CompressedVector tmp( serial( *this + (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF,Tag>::subAssign( const DenseVector<VT,TF>& rhs )
{
   using SubType = SubTrait_t< This, ResultType_t<VT> >;

   BLAZE_CONSTRAINT_MUST_BE_DENSE_VECTOR_TYPE( SubType );
   BLAZE_CONSTRAINT_MUST_BE_VECTOR_WITH_TRANSPOSE_FLAG( SubType, TF );
   BLAZE_CONSTRAINT_MUST_BE_REFERENCE_TYPE( CompositeType_t<SubType> );

   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   const SubType tmp( serial( *this - (*rhs) ) );
   reset();
   assign( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side sparse vector
inline void CompressedVector<Type,TF,Tag>::subAssign( const SparseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   CompressedVector tmp( serial( *this - (*rhs) ) );
   swap( tmp );
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF,Tag>::multAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   for( auto element=begin_; element!=end_; ++element ) {
      element->value_ *= (*rhs)[element->index_];
   }
}
//*************************************************************************************************


//*************************************************************************************************
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
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
template< typename VT >  // Type of the right-hand side dense vector
inline void CompressedVector<Type,TF,Tag>::divAssign( const DenseVector<VT,TF>& rhs )
{
   BLAZE_INTERNAL_ASSERT( size_ == (*rhs).size(), "Invalid vector sizes" );

   BLAZE_CONSTRAINT_MUST_NOT_REQUIRE_EVALUATION( VT );

   for( auto element=begin_; element!=end_; ++element ) {
      element->value_ /= (*rhs)[element->index_];
   }
}
//*************************************************************************************************




//=================================================================================================
//
//  COMPRESSEDVECTOR OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name CompressedVector operators */
//@{
template< typename Type, bool TF, typename Tag >
void reset( CompressedVector<Type,TF,Tag>& v );

template< typename Type, bool TF, typename Tag >
void clear( CompressedVector<Type,TF,Tag>& v );

template< RelaxationFlag RF, typename Type, bool TF, typename Tag >
bool isDefault( const CompressedVector<Type,TF,Tag>& v );

template< typename Type, bool TF, typename Tag >
bool isIntact( const CompressedVector<Type,TF,Tag>& v ) noexcept;

template< typename Type, bool TF, typename Tag >
void swap( CompressedVector<Type,TF,Tag>& a, CompressedVector<Type,TF,Tag>& b ) noexcept;
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the given compressed vector.
// \ingroup compressed_vector
//
// \param v The compressed vector to be resetted.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void reset( CompressedVector<Type,TF,Tag>& v )
{
   v.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the given compressed vector.
// \ingroup compressed_vector
//
// \param v The compressed vector to be cleared.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void clear( CompressedVector<Type,TF,Tag>& v )
{
   v.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the given compressed vector is in default state.
// \ingroup compressed_vector
//
// \param v The compressed vector to be tested for its default state.
// \return \a true in case the given vector's size is zero, \a false otherwise.
//
// This function checks whether the compressed vector is in default (constructed) state, i.e. if
// it's size is 0. In case it is in default state, the function returns \a true, else it will
// return \a false. The following example demonstrates the use of the \a isDefault() function:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isDefault( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isDefault<relaxed>( a ) ) { ... }
   \endcode
*/
template< RelaxationFlag RF  // Relaxation flag
        , typename Type      // Data type of the vector
        , bool TF            // Transpose flag
        , typename Tag >     // Type tag
inline bool isDefault( const CompressedVector<Type,TF,Tag>& v )
{
   return ( v.size() == 0UL );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the invariants of the given compressed vector are intact.
// \ingroup compressed_vector
//
// \param v The compressed vector to be tested.
// \return \a true in case the given vector's invariants are intact, \a false otherwise.
//
// This function checks whether the invariants of the compressed vector are intact, i.e. if
// its state is valid. In case the invariants are intact, the function returns \a true, else
// it will return \a false. The following example demonstrates the use of the \a isIntact()
// function:

   \code
   blaze::CompressedVector<double> a;
   // ... Resizing and initialization
   if( isIntact( a ) ) { ... }
   \endcode
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline bool isIntact( const CompressedVector<Type,TF,Tag>& v ) noexcept
{
   return ( v.nonZeros() <= v.capacity() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Swapping the contents of two compressed vectors.
// \ingroup compressed_vector
//
// \param a The first compressed vector to be swapped.
// \param b The second compressed vector to be swapped.
// \return void
*/
template< typename Type   // Data type of the vector
        , bool TF         // Transpose flag
        , typename Tag >  // Type tag
inline void swap( CompressedVector<Type,TF,Tag>& a, CompressedVector<Type,TF,Tag>& b ) noexcept
{
   a.swap( b );
}
//*************************************************************************************************




//=================================================================================================
//
//  ISRESIZABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF, typename Tag >
struct IsResizable< CompressedVector<T,TF,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ISSHRINKABLE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, bool TF, typename Tag >
struct IsShrinkable< CompressedVector<T,TF,Tag> >
   : public TrueType
{};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ADDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct AddTraitEval2< T1, T2
                    , EnableIf_t< IsSparseVector_v<T1> && IsSparseVector_v<T2> > >
{
   using Type = CompressedVector< AddTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , TransposeFlag_v<T1>
                                , AddTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct SubTraitEval2< T1, T2
                    , EnableIf_t< IsSparseVector_v<T1> && IsSparseVector_v<T2> > >
{
   using Type = CompressedVector< SubTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , TransposeFlag_v<T1>
                                , SubTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MULTTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsSparseVector_v<T1> && IsScalar_v<T2> > >
{
   using Type = CompressedVector< MultTrait_t< ElementType_t<T1>, T2 >
                                , TransposeFlag_v<T1>
                                , MultTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsScalar_v<T1> && IsSparseVector_v<T2> > >
{
   using Type = CompressedVector< MultTrait_t< T1, ElementType_t<T2> >
                                , TransposeFlag_v<T2>
                                , MultTrait_t< T1, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< ( ( IsRowVector_v<T1> && IsRowVector_v<T2> ) ||
                                     ( IsColumnVector_v<T1> && IsColumnVector_v<T2> ) ) &&
                                   ( IsSparseVector_v<T1> || IsSparseVector_v<T2> ) > >
{
   using Type = CompressedVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , TransposeFlag_v<T1>
                                , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsSparseMatrix_v<T1> &&
                                   IsSparseVector_v<T2> &&
                                   IsColumnVector_v<T2> > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = CompressedVector< AddTrait_t<MultType,MultType>
                                , false
                                , AddTrait_t<MultTag,MultTag> >;
};

template< typename T1, typename T2 >
struct MultTraitEval2< T1, T2
                     , EnableIf_t< IsSparseVector_v<T1> &&
                                   IsRowVector_v<T1> &&
                                   IsSparseMatrix_v<T2> > >
{
   using MultType = MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >;
   using MultTag  = MultTrait_t< TagType_t<T1>, TagType_t<T2> >;

   using Type = CompressedVector< AddTrait_t<MultType,MultType>
                                , true
                                , AddTrait_t<MultTag,MultTag> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  KRONTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1    // Type of the left-hand side operand
        , typename T2 >  // Type of the right-hand side operand
struct KronTraitEval2< T1, T2
                     , EnableIf_t< IsVector_v<T1> &&
                                   IsVector_v<T2> &&
                                   ( IsSparseVector_v<T1> || IsSparseVector_v<T2> ) > >
{
   using Type = CompressedVector< MultTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , ( IsDenseVector_v<T2> ? TransposeFlag_v<T1> : TransposeFlag_v<T2> )
                                , MultTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  DIVTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsSparseVector_v<T1> && IsScalar_v<T2> > >
{
   using Type = CompressedVector< DivTrait_t< ElementType_t<T1>, T2 >
                                , TransposeFlag_v<T1>
                                , DivTrait_t< TagType_t<T1>, T2 > >;
};

template< typename T1, typename T2 >
struct DivTraitEval2< T1, T2
                    , EnableIf_t< IsSparseVector_v<T1> && IsDenseVector_v<T2> > >
{
   using Type = CompressedVector< DivTrait_t< ElementType_t<T1>, ElementType_t<T2> >
                                , TransposeFlag_v<T1>
                                , DivTrait_t< TagType_t<T1>, TagType_t<T2> > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MAPTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, typename OP >
struct UnaryMapTraitEval2< T, OP
                         , EnableIf_t< IsSparseVector_v<T> > >
{
   using ElementType = decltype( std::declval<OP>()( std::declval< ElementType_t<T> >() ) );

   using Type = CompressedVector< EvaluateTrait_t<ElementType>
                                , TransposeFlag_v<T>
                                , MapTrait_t< TagType_t<T>, OP > >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  REPEATTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T, size_t R0 >
struct RepeatTraitEval2< T, R0, inf, inf
                       , EnableIf_t< IsSparseVector_v<T> > >
{
   using Type = CompressedVector< ElementType_t<T>
                                , TransposeFlag_v<T>
                                , TagType_t<T> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  HIGHTYPE SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool TF, typename Tag, typename T2 >
struct HighType< CompressedVector<T1,TF,Tag>, CompressedVector<T2,TF,Tag> >
{
   using Type = CompressedVector< typename HighType<T1,T2>::Type, TF, Tag >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  MATHTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename T1, bool TF, typename Tag, typename T2 >
struct LowType< CompressedVector<T1,TF,Tag>, CompressedVector<T2,TF,Tag> >
{
   using Type = CompressedVector< typename LowType<T1,T2>::Type, TF, Tag >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  SUBVECTORTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, size_t I, size_t N >
struct SubvectorTraitEval2< VT, I, N
                          , EnableIf_t< IsSparseVector_v<VT> > >
{
   using Type = CompressedVector< RemoveConst_t< ElementType_t<VT> >
                                , TransposeFlag_v<VT>
                                , TagType_t<VT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ELEMENTSTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename VT, size_t N >
struct ElementsTraitEval2< VT, N
                         , EnableIf_t< IsSparseVector_v<VT> > >
{
   using Type = CompressedVector< RemoveConst_t< ElementType_t<VT> >
                                , TransposeFlag_v<VT>
                                , TagType_t<VT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  ROWTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I >
struct RowTraitEval2< MT, I
                    , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedVector< RemoveConst_t< ElementType_t<MT> >
                                , true
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  COLUMNTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, size_t I >
struct ColumnTraitEval2< MT, I
                       , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedVector< RemoveConst_t< ElementType_t<MT> >
                                , false
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************




//=================================================================================================
//
//  BANDTRAIT SPECIALIZATIONS
//
//=================================================================================================

//*************************************************************************************************
/*! \cond BLAZE_INTERNAL */
template< typename MT, ptrdiff_t I >
struct BandTraitEval2< MT, I
                     , EnableIf_t< IsSparseMatrix_v<MT> > >
{
   using Type = CompressedVector< RemoveConst_t< ElementType_t<MT> >
                                , defaultTransposeFlag
                                , TagType_t<MT> >;
};
/*! \endcond */
//*************************************************************************************************

} // namespace blaze

#endif
