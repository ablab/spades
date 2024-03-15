//=================================================================================================
/*!
//  \file blaze/math/Vector.h
//  \brief Header file for all basic Vector functionality
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

#ifndef _BLAZE_MATH_VECTOR_H_
#define _BLAZE_MATH_VECTOR_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iomanip>
#include <iosfwd>
#include <blaze/math/Aliases.h>
#include <blaze/math/expressions/Vector.h>
#include <blaze/math/TransposeFlag.h>
#include <blaze/math/views/Elements.h>


namespace blaze {

//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Vector functions */
//@{
template< typename VT, bool TF >
bool isUniform( const Vector<VT,TF>& v );

template< typename T1, typename T2 >
decltype(auto) inner( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
decltype(auto) inner( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs );

template< typename T1, typename T2 >
decltype(auto) inner( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
decltype(auto) inner( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
decltype(auto) dot( const Vector<T1,TF1>& lhs, const Vector<T2,TF2>& rhs );

template< typename T1, bool TF1, typename T2, bool TF2 >
decltype(auto) operator,( const Vector<T1,TF1>& lhs, const Vector<T2,TF2>& rhs );

template< typename T1, typename T2 >
decltype(auto) outer( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
decltype(auto) outer( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs );

template< typename T1, typename T2 >
decltype(auto) outer( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs );

template< typename T1, typename T2 >
decltype(auto) outer( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs );

template< typename VT >
decltype(auto) reverse( VT&& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Checks if the given vector is a uniform vector.
// \ingroup vector
//
// \param v The vector to be checked.
// \return \a true if the vector is a uniform vector, \a false if not.
//
// This function checks if the given dense or sparse vector is a uniform vector. The vector
// is considered to be uniform if all its elements are identical. The following code example
// demonstrates the use of the function:

   \code
   blaze::DynamicVector<int,blaze::columnVector> a, b;
   // ... Initialization
   if( isUniform( a ) ) { ... }
   \endcode

// Optionally, it is possible to switch between strict semantics (blaze::strict) and relaxed
// semantics (blaze::relaxed):

   \code
   if( isUniform<relaxed>( a ) ) { ... }
   \endcode

// It is also possible to check if a vector expression results in a uniform vector:

   \code
   if( isUniform( a + b ) ) { ... }
   \endcode

// However, note that this might require the complete evaluation of the expression, including
// the generation of a temporary vector.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline bool isUniform( const Vector<VT,TF>& v )
{
   return isUniform<relaxed>( *v );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) inner( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs )
{
   return trans(*lhs) * (*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) inner( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs )
{
   return trans(*lhs) * trans(*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) inner( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs )
{
   return (*lhs) * (*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) inner( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs )
{
   return (*lhs) * trans(*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1  // Type of the left-hand side vector
        , bool TF1     // Transpose flag of the left-hand side vector
        , typename T2  // Type of the right-hand side vector
        , bool TF2 >   // Transpose flag of the right-hand side vector
inline decltype(auto) dot( const Vector<T1,TF1>& lhs, const Vector<T2,TF2>& rhs )
{
   return inner( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the scalar product (dot/inner product) of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the scalar product.
// \param rhs The right-hand side vector for the scalar product.
// \return The scalar product.
*/
template< typename T1  // Type of the left-hand side vector
        , bool TF1     // Transpose flag of the left-hand side vector
        , typename T2  // Type of the right-hand side vector
        , bool TF2 >   // Transpose flag of the right-hand side vector
inline decltype(auto) operator,( const Vector<T1,TF1>& lhs, const Vector<T2,TF2>& rhs )
{
   return inner( *lhs, *rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the outer product.
// \param rhs The right-hand side vector for the outer product.
// \return The outer product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) outer( const Vector<T1,false>& lhs, const Vector<T2,false>& rhs )
{
   return (*lhs) * trans(*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the outer product.
// \param rhs The right-hand side vector for the outer product.
// \return The outer product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) outer( const Vector<T1,false>& lhs, const Vector<T2,true>& rhs )
{
   return (*lhs) * (*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the outer product.
// \param rhs The right-hand side vector for the outer product.
// \return The outer product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) outer( const Vector<T1,true>& lhs, const Vector<T2,false>& rhs )
{
   return trans(*lhs) * trans(*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication operator for the outer product of two vectors
//        (\f$ s=(\vec{a},\vec{b}) \f$).
// \ingroup vector
//
// \param lhs The left-hand side vector for the outer product.
// \param rhs The right-hand side vector for the outer product.
// \return The outer product.
*/
template< typename T1    // Type of the left-hand side vector
        , typename T2 >  // Type of the right-hand side vector
inline decltype(auto) outer( const Vector<T1,true>& lhs, const Vector<T2,true>& rhs )
{
   return trans(*lhs) * (*rhs);
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Reverse the elements of a vector.
// \ingroup vector
//
// \param v The vector to be reversed.
// \return The reversed vector.
//
// This function reverses the elements of a dense or sparse vector. The following examples
// demonstrates this by means of a dense vector:

   \code
   blaze::DynamicVector<int> a{ 1, 2, 3, 4, 5 };
   blaze::DynamicVector<int> b;

   b = reverse( a );  // Results in ( 5 4 3 2 1 )
   \endcode
*/
template< typename VT >  // Type of the vector
inline decltype(auto) reverse( VT&& v )
{
   return elements( *v, [max=(*v).size()-1UL]( size_t i ){ return max - i; }, (*v).size() );
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Vector operators */
//@{
template< typename VT, bool TF >
std::ostream& operator<<( std::ostream& os, const Vector<VT,TF>& v );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for dense and sparse vectors.
// \ingroup vector
//
// \param os Reference to the output stream.
// \param v Reference to a constant vector object.
// \return Reference to the output stream.
*/
template< typename VT  // Type of the vector
        , bool TF >    // Transpose flag
inline std::ostream& operator<<( std::ostream& os, const Vector<VT,TF>& v )
{
   CompositeType_t<VT> tmp( *v );

   if( tmp.size() == 0UL ) {
      os << "( )\n";
   }
   else if( TF == rowVector ) {
      os << "(";
      for( size_t i=0UL; i<tmp.size(); ++i )
         os << " " << tmp[i];
      os << " )\n";
   }
   else {
      for( size_t i=0UL; i<tmp.size(); ++i )
         os << "( " << std::setw( 11UL ) << tmp[i] << " )\n";
   }

   return os;
}
//*************************************************************************************************

} // namespace blaze

#endif
