//=================================================================================================
/*!
//  \file blaze/math/sparse/ValueIndexPair.h
//  \brief Header file for the ValueIndexPair class
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

#ifndef _BLAZE_MATH_SPARSE_VALUEINDEXPAIR_H_
#define _BLAZE_MATH_SPARSE_VALUEINDEXPAIR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/sparse/SparseElement.h>
#include <blaze/math/typetraits/IsSparseElement.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/And.h>
#include <blaze/util/mpl/Not.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/RemoveReference.h>
#include <blaze/util/typetraits/IsRValueReference.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Index-value-pair for sparse vectors and matrices.
// \ingroup math
//
// The ValueIndexPair class represents a single index-value-pair of a sparse vector or sparse
// matrix.
*/
template< typename Type >  // Type of the value element
class ValueIndexPair
   : private SparseElement
{
 public:
   //**Type definitions****************************************************************************
   using ValueType      = Type;         //!< The value type of the value-index-pair.
   using IndexType      = size_t;       //!< The index type of the value-index-pair.
   using Reference      = Type&;        //!< Reference return type.
   using ConstReference = const Type&;  //!< Reference-to-const return type.
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   constexpr ValueIndexPair();
   constexpr ValueIndexPair( const Type& v, size_t i );

   ValueIndexPair( const ValueIndexPair& ) = default;
   ValueIndexPair( ValueIndexPair&& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ValueIndexPair() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
   ValueIndexPair& operator=( const ValueIndexPair& ) = default;
   ValueIndexPair& operator=( ValueIndexPair&& ) = default;

   template< typename Other >
   constexpr auto operator=( const Other& rhs )
      -> EnableIf_t< IsSparseElement_v<Other>, ValueIndexPair& >;

   template< typename Other >
   constexpr auto operator=( Other&& rhs )
      -> EnableIf_t< IsSparseElement_v< RemoveReference_t<Other> > &&
                     IsRValueReference_v<Other&&>, ValueIndexPair& >;

   template< typename Other >
   constexpr auto operator=( const Other& v )
      -> EnableIf_t< !IsSparseElement_v<Other>, ValueIndexPair& >;

   template< typename Other >
   constexpr auto operator=( Other&& v )
      -> EnableIf_t< !IsSparseElement_v< RemoveReference_t<Other> > &&
                     IsRValueReference_v<Other&&>, ValueIndexPair& >;

   template< typename Other > constexpr ValueIndexPair& operator+=( const Other& v );
   template< typename Other > constexpr ValueIndexPair& operator-=( const Other& v );
   template< typename Other > constexpr ValueIndexPair& operator*=( const Other& v );
   template< typename Other > constexpr ValueIndexPair& operator/=( const Other& v );
   //@}
   //**********************************************************************************************

   //**Acess functions*****************************************************************************
   /*!\name Access functions */
   //@{
   constexpr Reference      value() noexcept;
   constexpr ConstReference value() const noexcept;
   constexpr IndexType      index() const noexcept;
   //@}
   //**********************************************************************************************

 protected:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   Type   value_;  //!< Value of the value-index-pair.
   size_t index_;  //!< Index of the value-index-pair.
   //@}
   //**********************************************************************************************

 private:
   //**Friend declarations*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   template< typename Other > friend class ValueIndexPair;
   /*! \endcond */
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE  ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST         ( Type );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE      ( Type );
   /*! \endcond */
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Default constructor for value-index-pairs.
*/
template< typename Type >  // Type of the value element
constexpr ValueIndexPair<Type>::ValueIndexPair()
   : value_()  // Value of the value-index-pair
   , index_()  // Index of the value-index-pair
{}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Constructor for a direct initialization of value-index-pairs.
//
// \param v The value of the value-index-pair.
// \param i The index of the value-index-pair.
*/
template< typename Type >  // Type of the value element
constexpr ValueIndexPair<Type>::ValueIndexPair( const Type& v, size_t i )
   : value_( v )  // Value of the value-index-pair
   , index_( i )  // Index of the value-index-pair
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Assignment operator for different value-index-pair types.
//
// \param rhs Value-index-pair to be copied.
// \return Reference to the assigned value-index-pair.
//
// This assignment operator enables the assignment of other value-index-pair types. The given
// \a Other data type qualifies as value-index-pair type in case it provides a value() and an
// index() member function.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value-index-pair
constexpr auto ValueIndexPair<Type>::operator=( const Other& rhs )
   -> EnableIf_t< IsSparseElement_v<Other>, ValueIndexPair& >
{
   value_ = rhs.value();
   index_ = rhs.index();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment operator for different value-index-pair types.
//
// \param rhs Value-index-pair to be moved.
// \return Reference to the assigned value-index-pair.
//
// This assignment operator enables the assignment of other value-index-pair types. The given
// \a Other data type qualifies as value-index-pair type in case it provides a value() and an
// index() member function.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value-index-pair
constexpr auto ValueIndexPair<Type>::operator=( Other&& rhs )
   -> EnableIf_t< IsSparseElement_v< RemoveReference_t<Other> > &&
                  IsRValueReference_v<Other&&>, ValueIndexPair& >
{
   value_ = std::move( rhs.value() );
   index_ = rhs.index();
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the value of the value-index-pair.
//
// \param v The new value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr auto ValueIndexPair<Type>::operator=( const Other& v )
   -> EnableIf_t< !IsSparseElement_v<Other>, ValueIndexPair& >
{
   value_ = v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the value of the value-index-pair.
//
// \param v The new value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr auto ValueIndexPair<Type>::operator=( Other&& v )
   -> EnableIf_t< !IsSparseElement_v< RemoveReference_t<Other> > &&
                  IsRValueReference_v<Other&&>, ValueIndexPair& >
{
   value_ = std::move( v );
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the value of the value-index-pair.
//
// \param v The right-hand side value to be added to the value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr ValueIndexPair<Type>& ValueIndexPair<Type>::operator+=( const Other& v )
{
   value_ += v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the value of the value-index-pair.
//
// \param v The right-hand side value to be subtracted from the value-index-pair value.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr ValueIndexPair<Type>& ValueIndexPair<Type>::operator-=( const Other& v )
{
   value_ -= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the value of the value-index-pair.
//
// \param v The right-hand side value for the multiplication.
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr ValueIndexPair<Type>& ValueIndexPair<Type>::operator*=( const Other& v )
{
   value_ *= v;
   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the value of the value-index-pair.
//
// \param v The right-hand side value for the division
// \return Reference to the assigned value-index-pair.
*/
template< typename Type >   // Type of the value element
template< typename Other >  // Data type of the right-hand side value
constexpr ValueIndexPair<Type>& ValueIndexPair<Type>::operator/=( const Other& v )
{
   value_ /= v;
   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access to the current value of the value-index-pair.
//
// \return The current value of the value-index-pair.
*/
template< typename Type >  // Type of the value element
constexpr typename ValueIndexPair<Type>::Reference
   ValueIndexPair<Type>::value() noexcept
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current value of the value-index-pair.
//
// \return The current value of the value-index-pair.
*/
template< typename Type >  // Type of the value element
constexpr typename ValueIndexPair<Type>::ConstReference
   ValueIndexPair<Type>::value() const noexcept
{
   return value_;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Access to the current index of the value-index-pair.
//
// \return The current index of the value-index-pair.
*/
template< typename Type >  // Type of the value element
constexpr typename ValueIndexPair<Type>::IndexType
   ValueIndexPair<Type>::index() const noexcept
{
   return index_;
}
//*************************************************************************************************

} // namespace blaze

#endif
