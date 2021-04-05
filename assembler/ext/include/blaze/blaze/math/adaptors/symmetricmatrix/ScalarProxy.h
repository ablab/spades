//=================================================================================================
/*!
//  \file blaze/math/adaptors/symmetricmatrix/ScalarProxy.h
//  \brief Header file for the ScalarProxy class
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

#ifndef _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_SCALARPROXY_H_
#define _BLAZE_MATH_ADAPTORS_SYMMETRICMATRIX_SCALARPROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/math/Aliases.h>
#include <blaze/math/constraints/Computation.h>
#include <blaze/math/constraints/Hermitian.h>
#include <blaze/math/constraints/Lower.h>
#include <blaze/math/constraints/Matrix.h>
#include <blaze/math/constraints/Scalar.h>
#include <blaze/math/constraints/Symmetric.h>
#include <blaze/math/constraints/Transformation.h>
#include <blaze/math/constraints/Upper.h>
#include <blaze/math/constraints/View.h>
#include <blaze/math/proxy/Proxy.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Clear.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/IsDefault.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/IsReal.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Reset.h>
#include <blaze/util/constraints/Const.h>
#include <blaze/util/constraints/Pointer.h>
#include <blaze/util/constraints/Reference.h>
#include <blaze/util/constraints/Volatile.h>
#include <blaze/util/InvalidType.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/Types.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Access proxy for symmetric, square matrices with scalar element types.
// \ingroup symmetric_matrix
//
// The ScalarProxy provides controlled access to the elements of a non-const symmetric matrix
// with scalar element type (e.g. integral values, floating point values, and complex values).
// It guarantees that a modification of element \f$ a_{ij} \f$ of the accessed matrix is also
// applied to element \f$ a_{ji} \f$. The following example illustrates this by means of a
// \f$ 3 \times 3 \f$ dense symmetric matrix:

   \code
   // Creating a 3x3 symmetric dense matrix
   blaze::SymmetricMatrix< blaze::DynamicMatrix<int> > A( 3UL );

   A(0,2) = -2;  //        (  0 0 -2 )
   A(1,1) =  3;  // => A = (  0 3  5 )
   A(1,2) =  5;  //        ( -2 5  0 )
   \endcode
*/
template< typename MT >  // Type of the adapted matrix
class ScalarProxy
   : public Proxy< ScalarProxy<MT> >
{
 private:
   //**struct BuiltinType**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary struct to determine the value type of the represented complex element.
   */
   template< typename T >
   struct BuiltinType { using Type = INVALID_TYPE; };
   /*! \endcond */
   //**********************************************************************************************

   //**struct ComplexType**************************************************************************
   /*! \cond BLAZE_INTERNAL */
   /*!\brief Auxiliary struct to determine the value type of the represented complex element.
   */
   template< typename T >
   struct ComplexType { using Type = typename T::value_type; };
   /*! \endcond */
   //**********************************************************************************************

 public:
   //**Type definitions****************************************************************************
   using RepresentedType = ElementType_t<MT>;     //!< Type of the represented matrix element.
   using Reference       = Reference_t<MT>;       //!< Reference to the represented element.
   using ConstReference  = ConstReference_t<MT>;  //!< Reference-to-const to the represented element.
   using Pointer         = ScalarProxy*;         //!< Pointer to the represented element.
   using ConstPointer    = const ScalarProxy*;   //!< Pointer-to-const to the represented element.

   //! Value type of the represented complex element.
   using ValueType = typename If_t< IsComplex_v<RepresentedType>
                                  , ComplexType<RepresentedType>
                                  , BuiltinType<RepresentedType> >::Type;
   //**********************************************************************************************

   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline ScalarProxy( MT& matrix, size_t row, size_t column );

   ScalarProxy( const ScalarProxy& ) = default;
   //@}
   //**********************************************************************************************

   //**Destructor**********************************************************************************
   /*!\name Destructor */
   //@{
   ~ScalarProxy() = default;
   //@}
   //**********************************************************************************************

   //**Assignment operators************************************************************************
   /*!\name Assignment operators */
   //@{
                          inline ScalarProxy& operator= ( const ScalarProxy& sp );
   template< typename T > inline ScalarProxy& operator= ( const T& value );
   template< typename T > inline ScalarProxy& operator+=( const T& value );
   template< typename T > inline ScalarProxy& operator-=( const T& value );
   template< typename T > inline ScalarProxy& operator*=( const T& value );
   template< typename T > inline ScalarProxy& operator/=( const T& value );
   template< typename T > inline ScalarProxy& operator%=( const T& value );
   //@}
   //**********************************************************************************************

   //**Access operators****************************************************************************
   /*!\name Access operators */
   //@{
   inline Pointer      operator->();
   inline ConstPointer operator->() const;
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline void reset () const;
   inline void clear () const;
   inline void invert() const;

   inline ConstReference get() const noexcept;
   //@}
   //**********************************************************************************************

   //**Conversion operator*************************************************************************
   /*!\name Conversion operator */
   //@{
   inline operator ConstReference() const noexcept;
   //@}
   //**********************************************************************************************

   //**Complex data access functions***************************************************************
   /*!\name Complex data access functions */
   //@{
   inline ValueType real() const;
   inline void      real( ValueType value ) const;
   inline ValueType imag() const;
   inline void      imag( ValueType value ) const;
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   MT&    matrix_;  //!< Reference to the adapted matrix.
   size_t row_;     //!< Row index of the accessed matrix element.
   size_t column_;  //!< Column index of the accessed matrix element.
   //@}
   //**********************************************************************************************

   //**Compile time checks*************************************************************************
   /*! \cond BLAZE_INTERNAL */
   BLAZE_CONSTRAINT_MUST_BE_MATRIX_TYPE              ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_REFERENCE_TYPE       ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_POINTER_TYPE         ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_CONST                ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VOLATILE             ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_VIEW_TYPE            ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_COMPUTATION_TYPE     ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_TRANSFORMATION_TYPE  ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_SYMMETRIC_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_HERMITIAN_MATRIX_TYPE( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_LOWER_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_NOT_BE_UPPER_MATRIX_TYPE    ( MT );
   BLAZE_CONSTRAINT_MUST_BE_SCALAR_TYPE              ( RepresentedType );
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
/*!\brief Initialization constructor for a ScalarProxy.
//
// \param matrix Reference to the adapted matrix.
// \param row The row-index of the accessed matrix element.
// \param column The column-index of the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline ScalarProxy<MT>::ScalarProxy( MT& matrix, size_t row, size_t column )
   : matrix_( matrix )  // Reference to the adapted matrix
   , row_   ( row    )  // Row index of the accessed matrix element
   , column_( column )  // Column index of the accessed matrix element
{}
//*************************************************************************************************




//=================================================================================================
//
//  OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Copy assignment operator for ScalarProxy.
//
// \param sp Scalar proxy to be copied.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
inline ScalarProxy<MT>& ScalarProxy<MT>::operator=( const ScalarProxy& sp )
{
   matrix_(row_,column_) = sp.matrix_(sp.row_,sp.column_);
   matrix_(column_,row_) = sp.matrix_(sp.row_,sp.column_);

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Assignment to the accessed matrix element.
//
// \param value The new value of the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator=( const T& value )
{
   matrix_(row_,column_) = value;
   if( row_ != column_ )
      matrix_(column_,row_) = value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment to the accessed matrix element.
//
// \param value The right-hand side value to be added to the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator+=( const T& value )
{
   matrix_(row_,column_) += value;
   if( row_ != column_ )
      matrix_(column_,row_) += value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment to the accessed matrix element.
//
// \param value The right-hand side value to be subtracted from the matrix element.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator-=( const T& value )
{
   matrix_(row_,column_) -= value;
   if( row_ != column_ )
      matrix_(column_,row_) -= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment to the accessed matrix element.
//
// \param value The right-hand side value for the multiplication.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator*=( const T& value )
{
   matrix_(row_,column_) *= value;
   if( row_ != column_ )
      matrix_(column_,row_) *= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment to the accessed matrix element.
//
// \param value The right-hand side value for the division.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator/=( const T& value )
{
   matrix_(row_,column_) /= value;
   if( row_ != column_ )
      matrix_(column_,row_) /= value;

   return *this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Modulo assignment to the accessed matrix element.
//
// \param value The right-hand side value for the modulo operation.
// \return Reference to the assigned proxy.
*/
template< typename MT >  // Type of the adapted matrix
template< typename T >   // Type of the right-hand side value
inline ScalarProxy<MT>& ScalarProxy<MT>::operator%=( const T& value )
{
   matrix_(row_,column_) %= value;
   if( row_ != column_ )
      matrix_(column_,row_) %= value;

   return *this;
}
//*************************************************************************************************




//=================================================================================================
//
//  ACCESS OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Direct access to the represented matrix element.
//
// \return Pointer to the represented matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename ScalarProxy<MT>::Pointer ScalarProxy<MT>::operator->()
{
   return this;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Direct access to the represented matrix element.
//
// \return Pointer to the represented matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename ScalarProxy<MT>::ConstPointer ScalarProxy<MT>::operator->() const
{
   return this;
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Reset the represented element to its default initial value.
//
// \return void
//
// This function resets the element represented by the proxy to its default initial value.
*/
template< typename MT >  // Type of the adapted matrix
inline void ScalarProxy<MT>::reset() const
{
   using blaze::reset;

   reset( matrix_(row_,column_) );
   if( row_ != column_ )
      reset( matrix_(column_,row_) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
//
// \return void
//
// This function clears the element represented by the proxy to its default initial state.
*/
template< typename MT >  // Type of the adapted matrix
inline void ScalarProxy<MT>::clear() const
{
   using blaze::clear;

   clear( matrix_(row_,column_) );
   if( row_ != column_ )
      clear( matrix_(column_,row_) );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element
//
// \return void
*/
template< typename MT >  // Type of the adapted matrix
inline void ScalarProxy<MT>::invert() const
{
   using blaze::invert;

   invert( matrix_(row_,column_) );
   if( row_ != column_ )
      matrix_(column_,row_) = matrix_(row_,column_);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returning the value of the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline typename ScalarProxy<MT>::ConstReference ScalarProxy<MT>::get() const noexcept
{
   return const_cast<const MT&>( matrix_ )(row_,column_);
}
//*************************************************************************************************




//=================================================================================================
//
//  CONVERSION OPERATOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Conversion to the accessed matrix element.
//
// \return Direct/raw reference to the accessed matrix element.
*/
template< typename MT >  // Type of the adapted matrix
inline ScalarProxy<MT>::operator ConstReference() const noexcept
{
   return get();
}
//*************************************************************************************************




//=================================================================================================
//
//  COMPLEX DATA ACCESS FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the real part of the represented complex number.
//
// \return The current real part of the represented complex number.
//
// In case the proxy represents a complex number, this function returns the current value of its
// real part.
*/
template< typename MT >  // Type of the adapted matrix
inline typename ScalarProxy<MT>::ValueType ScalarProxy<MT>::real() const
{
   return matrix_(row_,column_).real();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the real part of the represented complex number.
//
// \param value The new value for the real part.
// \return void
//
// In case the proxy represents a complex number, this function sets a new value to its real part.
*/
template< typename MT >  // Type of the adapted matrix
inline void ScalarProxy<MT>::real( ValueType value ) const
{
   matrix_(row_,column_).real( value );
   if( row_ != column_ )
      matrix_(column_,row_).real( value );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns the imaginary part of the represented complex number.
//
// \return The current imaginary part of the represented complex number.
//
// In case the proxy represents a complex number, this function returns the current value of its
// imaginary part.
*/
template< typename MT >  // Type of the adapted matrix
inline typename ScalarProxy<MT>::ValueType ScalarProxy<MT>::imag() const
{
   return matrix_(row_,column_).imag();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Setting the imaginary part of the represented complex number.
//
// \param value The new value for the imaginary part.
// \return void
//
// In case the proxy represents a complex number, this function sets a new value to its imaginary
// part.
*/
template< typename MT >  // Type of the adapted matrix
inline void ScalarProxy<MT>::imag( ValueType value ) const
{
   matrix_(row_,column_).imag( value );
   if( row_ != column_ )
      matrix_(column_,row_).imag( value );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name ScalarProxy global functions */
//@{
template< typename MT >
void reset( const ScalarProxy<MT>& proxy );

template< typename MT >
void clear( const ScalarProxy<MT>& proxy );

template< typename MT >
void invert( const ScalarProxy<MT>& proxy );

template< RelaxationFlag RF, typename MT >
bool isDefault( const ScalarProxy<MT>& proxy );

template< RelaxationFlag RF, typename MT >
bool isReal( const ScalarProxy<MT>& proxy );

template< RelaxationFlag RF, typename MT >
bool isZero( const ScalarProxy<MT>& proxy );

template< RelaxationFlag RF, typename MT >
bool isOne( const ScalarProxy<MT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Resetting the represented element to the default initial values.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function resets the element represented by the scalar proxy to its default initial
// value.
*/
template< typename MT >
inline void reset( const ScalarProxy<MT>& proxy )
{
   proxy.reset();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Clearing the represented element.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return void
//
// This function clears the element represented by the scalar proxy to its default initial
// state.
*/
template< typename MT >
inline void clear( const ScalarProxy<MT>& proxy )
{
   proxy.clear();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element.
// \ingroup symmetric_matrix
//
// \param proxy The given proxy instance.
// \return void
*/
template< typename MT >
inline void invert( const ScalarProxy<MT>& proxy )
{
   proxy.invert();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is in default state.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy
// \return \a true in case the represented element is in default state, \a false otherwise.
//
// This function checks whether the element represented by the access proxy is in default state.
// In case it is in default state, the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename MT >
inline bool isDefault( const ScalarProxy<MT>& proxy )
{
   using blaze::isDefault;

   return isDefault<RF>( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the matrix element represents a real number.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return \a true in case the matrix element represents a real number, \a false otherwise.
//
// This function checks whether the element represented by the access proxy represents the a
// real number. In case the element is of built-in type, the function returns \a true. In case
// the element is of complex type, the function returns \a true if the imaginary part is equal
// to 0. Otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename MT >
inline bool isReal( const ScalarProxy<MT>& proxy )
{
   using blaze::isReal;

   return isReal<RF>( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 0.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return \a true in case the represented element is 0, \a false otherwise.
//
// This function checks whether the element represented by the access proxy represents the numeric
// value 0. In case it is 0, the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename MT >
inline bool isZero( const ScalarProxy<MT>& proxy )
{
   using blaze::isZero;

   return isZero<RF>( proxy.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 1.
// \ingroup symmetric_matrix
//
// \param proxy The given access proxy.
// \return \a true in case the represented element is 1, \a false otherwise.
//
// This function checks whether the element represented by the access proxy represents the numeric
// value 1. In case it is 1, the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename MT >
inline bool isOne( const ScalarProxy<MT>& proxy )
{
   using blaze::isOne;

   return isOne<RF>( proxy.get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
