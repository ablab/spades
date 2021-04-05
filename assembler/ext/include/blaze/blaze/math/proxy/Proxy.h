//=================================================================================================
/*!
//  \file blaze/math/proxy/Proxy.h
//  \brief Header file for the Proxy class
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

#ifndef _BLAZE_MATH_PROXY_PROXY_H_
#define _BLAZE_MATH_PROXY_PROXY_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <iosfwd>
#include <blaze/math/Aliases.h>
#include <blaze/math/Exception.h>
#include <blaze/math/InversionFlag.h>
#include <blaze/math/proxy/ComplexProxy.h>
#include <blaze/math/proxy/DefaultProxy.h>
#include <blaze/math/proxy/DenseMatrixProxy.h>
#include <blaze/math/proxy/DenseVectorProxy.h>
#include <blaze/math/proxy/SparseMatrixProxy.h>
#include <blaze/math/proxy/SparseVectorProxy.h>
#include <blaze/math/RelaxationFlag.h>
#include <blaze/math/shims/Abs.h>
#include <blaze/math/shims/Acos.h>
#include <blaze/math/shims/Acosh.h>
#include <blaze/math/shims/Asin.h>
#include <blaze/math/shims/Asinh.h>
#include <blaze/math/shims/Atan.h>
#include <blaze/math/shims/Atan2.h>
#include <blaze/math/shims/Atanh.h>
#include <blaze/math/shims/Cbrt.h>
#include <blaze/math/shims/Ceil.h>
#include <blaze/math/shims/Conjugate.h>
#include <blaze/math/shims/Cos.h>
#include <blaze/math/shims/Cosh.h>
#include <blaze/math/shims/Erf.h>
#include <blaze/math/shims/Erfc.h>
#include <blaze/math/shims/Evaluate.h>
#include <blaze/math/shims/Exp.h>
#include <blaze/math/shims/Exp2.h>
#include <blaze/math/shims/Exp10.h>
#include <blaze/math/shims/Floor.h>
#include <blaze/math/shims/Imaginary.h>
#include <blaze/math/shims/InvCbrt.h>
#include <blaze/math/shims/Invert.h>
#include <blaze/math/shims/InvSqrt.h>
#include <blaze/math/shims/IsFinite.h>
#include <blaze/math/shims/IsInf.h>
#include <blaze/math/shims/IsNaN.h>
#include <blaze/math/shims/IsOne.h>
#include <blaze/math/shims/IsReal.h>
#include <blaze/math/shims/IsZero.h>
#include <blaze/math/shims/Log.h>
#include <blaze/math/shims/Log2.h>
#include <blaze/math/shims/Log10.h>
#include <blaze/math/shims/Pow.h>
#include <blaze/math/shims/Pow2.h>
#include <blaze/math/shims/Pow3.h>
#include <blaze/math/shims/Pow4.h>
#include <blaze/math/shims/Real.h>
#include <blaze/math/shims/Round.h>
#include <blaze/math/shims/Sign.h>
#include <blaze/math/shims/Sin.h>
#include <blaze/math/shims/Sinh.h>
#include <blaze/math/shims/Sqrt.h>
#include <blaze/math/shims/Tan.h>
#include <blaze/math/shims/Tanh.h>
#include <blaze/math/shims/Trunc.h>
#include <blaze/math/typetraits/IsDenseMatrix.h>
#include <blaze/math/typetraits/IsDenseVector.h>
#include <blaze/math/typetraits/IsMatrix.h>
#include <blaze/math/typetraits/IsProxy.h>
#include <blaze/math/typetraits/IsVector.h>
#include <blaze/util/algorithms/Max.h>
#include <blaze/util/algorithms/Min.h>
#include <blaze/util/EnableIf.h>
#include <blaze/util/mpl/If.h>
#include <blaze/util/typetraits/IsComplex.h>


namespace blaze {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Proxy base class.
// \ingroup math
//
// The Proxy class is a base class for all proxy classes within the \b Blaze library that may
// represent non-numeric data types (vectors, matrices, ...). It augments the interface of the
// deriving proxy class depending on the data type represented by the proxy. In addition, it
// provides an abstraction from the actual type of the proxy, but enables a type-safe conversion
// back to this type via the 'Curiously Recurring Template Pattern' (CRTP).
//
// In order to use the Proxy class it is necessary to publicly derive from it and to provide
// an accessible member function called \a get(), which grants access to the represented element
// via non-const reference. The following example demonstrates these requirements by means of
// the VectorAccessProxy class:

   \code
   template< typename VT >
   class VectorAccessProxy : public Proxy< VectorAccessProxy<VT>, typename VT::ElementType >
   {
      // ...
      using RepresentedType = typename VT::ElementType;
      inline RepresentedType& get() const;
      // ...
   };
   \endcode

// The first template parameter specifies the type of the deriving proxy class (CRTP), the second
// template parameter specifies the type of the element represented by the proxy. Within the
// context of the VectorAccessProxy this is the type of the elements of the vector to be accessed.
// Depending on this type the proxy selects the additional interface to provide to the deriving
// class.
*/
template< typename PT          // Type of the proxy
        , typename RT = int >  // Type of the represented element
class Proxy
   : public If_t< IsVector_v<RT>
                , If_t< IsDenseVector_v<RT>
                      , DenseVectorProxy<PT,RT>
                      , SparseVectorProxy<PT,RT> >
                , If_t< IsMatrix_v<RT>
                      , If_t< IsDenseMatrix_v<RT>
                            , DenseMatrixProxy<PT,RT>
                            , SparseMatrixProxy<PT,RT> >
                      , If_t< IsComplex_v<RT>
                            , ComplexProxy<PT,RT>
                            , DefaultProxy<PT,RT> > > >
{
 protected:
   //**Special member functions********************************************************************
   /*!\name Special member functions */
   //@{
   Proxy() = default;
   Proxy( const Proxy& ) = default;
   Proxy( Proxy&& ) = default;
   ~Proxy() = default;
   Proxy& operator=( const Proxy& ) = default;
   Proxy& operator=( Proxy&& ) = default;
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL OPERATORS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Proxy operators */
//@{
template< typename T, typename PT, typename RT >
decltype(auto) operator+=( T& lhs, const Proxy<PT,RT>& rhs );

template< typename T, typename PT, typename RT >
decltype(auto) operator-=( T& lhs, const Proxy<PT,RT>& rhs );

template< typename T, typename PT, typename RT >
decltype(auto) operator/=( T& lhs, const Proxy<PT,RT>& rhs );

template< typename T, typename PT, typename RT >
decltype(auto) operator*=( T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) operator+( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator+( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator+( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) operator-( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator-( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator-( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) operator*( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator*( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator*( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) operator/( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator/( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) operator/( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator==( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator==( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator==( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator!=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator!=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator!=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator<( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator<( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator<( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator>( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator>( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator>( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator<=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator<=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator<=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
bool operator>=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool operator>=( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool operator>=( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT, typename RT >
std::ostream& operator<<( std::ostream& os, const Proxy<PT,RT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition assignment of an object of different type with a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the addition assignment.
*/
template< typename T, typename PT, typename RT >
inline decltype(auto) operator+=( T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs += (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction assignment of an object of different type with a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the subtraction assignment.
*/
template< typename T, typename PT, typename RT >
inline decltype(auto) operator-=( T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs -= (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication assignment of an object of different type with a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the multiplication assignment.
*/
template< typename T, typename PT, typename RT >
inline decltype(auto) operator*=( T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs *= (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division assignment of an object of different type with a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the division assignment.
*/
template< typename T, typename PT, typename RT >
inline decltype(auto) operator/=( T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs /= (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the addition.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) operator+( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (*lhs).get() + (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the addition.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) operator+( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (*lhs).get() + rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Addition between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the addition.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) operator+( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs + (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the subtraction.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) operator-( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (*lhs).get() - (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the subtraction.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) operator-( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (*lhs).get() - rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Subtraction between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the subtraction.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) operator-( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs - (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the multiplication.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) operator*( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (*lhs).get() * (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the multiplication.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) operator*( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (*lhs).get() * rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Multiplication between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the multiplication.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) operator*( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs * (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The result of the division.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) operator/( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return (*lhs).get() / (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The result of the division.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) operator/( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return (*lhs).get() / rhs;
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Division between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The result of the division.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) operator/( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return lhs / (*rhs).get();
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator==( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() == (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator==( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() == rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator==( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs == (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if both referenced values are not equal, \a false if they are.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator!=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() != (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inequality comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are not equal, \a false if they are.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator!=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() != rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Inquality comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the other object and the referenced value are not equal, \a false if they are.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator!=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs != (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() < (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller, \a false if not.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator<( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() < rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is smaller, \a false if not.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator<( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs < rhs.get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() > (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater, \a false if not.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator>( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() > rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is greater, \a false if not.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator>( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs > (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator<=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() <= (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is smaller or equal, \a false if not.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator<=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() <= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Less-or-equal-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is smaller or equal, \a false if not.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator<=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs <= (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool operator>=( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return ( (*lhs).get() >= (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the left-hand side referenced value is greater or equal, \a false if not.
*/
template< typename PT, typename RT, typename T, typename >
inline bool operator>=( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return ( (*lhs).get() >= rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Greater-or-equal-than comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the left-hand side other object is greater or equal, \a false if not.
*/
template< typename T, typename PT, typename RT, typename >
inline bool operator>=( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return ( lhs >= (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Global output operator for the Proxy class template.
// \ingroup math
//
// \param os Reference to the output stream.
// \param proxy Reference to a constant proxy object.
// \return Reference to the output stream.
*/
template< typename PT, typename RT >
inline std::ostream& operator<<( std::ostream& os, const Proxy<PT,RT>& proxy )
{
   return os << (*proxy).get();
}
//*************************************************************************************************




//=================================================================================================
//
//  GLOBAL FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\name Proxy global functions */
//@{
template< typename PT, typename RT >
decltype(auto) abs( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) sign( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) floor( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) ceil( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) trunc( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) round( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) conj( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) trans( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) ctrans( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) real( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) imag( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) sqrt( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) invsqrt( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) cbrt( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) invcbrt( const Proxy<PT,RT>& proxy );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) hypot( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) hypot( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) hypot( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT, typename RT, typename ET >
decltype(auto) pow( const Proxy<PT,RT>& proxy, const ET& exp );

template< typename PT, typename RT >
decltype(auto) pow2( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) pow3( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) pow4( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) exp( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) exp2( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) exp10( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) log( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) log2( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) log10( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) sin( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) asin( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) sinh( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) asinh( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) cos( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) acos( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) cosh( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) acosh( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) tan( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) atan( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) tanh( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) atanh( const Proxy<PT,RT>& proxy );

template< typename PT1, typename RT1, typename PT2, typename RT2 >
decltype(auto) atan2( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) atan2( const Proxy<PT,RT>& lhs, const T& rhs );

template< typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
decltype(auto) atan2( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT, typename RT >
decltype(auto) erf( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
decltype(auto) erfc( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
void transpose( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
void ctranspose( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
void invert( const Proxy<PT,RT>& proxy );

template< InversionFlag IF, typename PT, typename RT >
void invert( const Proxy<PT,RT>& proxy );

template< RelaxationFlag RF, typename PT, typename RT >
bool isReal( const Proxy<PT,RT>& proxy );

template< RelaxationFlag RF, typename PT, typename RT >
bool isZero( const Proxy<PT,RT>& proxy );

template< RelaxationFlag RF, typename PT, typename RT >
bool isOne( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
bool isnan( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
bool isinf( const Proxy<PT,RT>& proxy );

template< typename PT, typename RT >
bool isfinite( const Proxy<PT,RT>& proxy );

template< RelaxationFlag RF, typename PT1, typename RT1, typename PT2, typename RT2 >
bool equal( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs );

template< RelaxationFlag RF, typename PT, typename RT, typename T, typename = DisableIf_t< IsProxy_v<T> > >
bool equal( const Proxy<PT,RT>& lhs, const T& rhs );

template< RelaxationFlag RF, typename T, typename PT, typename RT, typename = DisableIf_t< IsProxy_v<T> > >
bool equal( const T& lhs, const Proxy<PT,RT>& rhs );

template< typename PT, typename RT >
decltype(auto) evaluate( const Proxy<PT,RT>& proxy );
//@}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the absolute value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The absolute value of the represented element.
//
// This function computes the absolute value of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the absolute values of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) abs( const Proxy<PT,RT>& proxy )
{
   using blaze::abs;

   return abs( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Evaluating the sign of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return 1 if the value is greater than zero, 0 if it is zero, and -1 if it is less than zero.
//
// This function evaluates the sign of the element represented by the proxy. In case the proxy
// represents a vector- or matrix-like data structure the function returns an expression
// representing the operation.
*/
template< typename PT, typename RT >
inline decltype(auto) sign( const Proxy<PT,RT>& proxy )
{
   using blaze::sign;

   return sign( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the largest integral value that is not greater than the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The largest integral value that is not greater than the represented element.
//
// This function computes the largest integral value that is not greater than the element
// represented by the proxy. In case the proxy represents a vector- or matrix-like data
// structure the function returns an expression representing the operation.
*/
template< typename PT, typename RT >
inline decltype(auto) floor( const Proxy<PT,RT>& proxy )
{
   using blaze::floor;

   return floor( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the smallest integral value that is not less than the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The smallest integral value that is not less than the represented element.
//
// This function computes the smallest integral value that is not less than the element
// represented by the proxy. In case the proxy represents a vector- or matrix-like data
// structure the function returns an expression representing the operation.
*/
template< typename PT, typename RT >
inline decltype(auto) ceil( const Proxy<PT,RT>& proxy )
{
   using blaze::ceil;

   return ceil( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the nearest integral value that is not greater than the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The nearest integral value that is not greater than the represented element.
//
// This function computes the nearest integral value that is not greater than the element
// represented by the proxy. In case the proxy represents a vector- or matrix-like data
// structure the function returns an expression representing the operation.
*/
template< typename PT, typename RT >
inline decltype(auto) trunc( const Proxy<PT,RT>& proxy )
{
   using blaze::trunc;

   return trunc( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the nearest integral value to the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The nearest integral value to the represented element.
//
// This function computes the nearest integral value to the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the operation.
*/
template< typename PT, typename RT >
inline decltype(auto) round( const Proxy<PT,RT>& proxy )
{
   using blaze::round;

   return round( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the complex conjugate of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The complex conjugate of the represented element.
//
// This function computes the complex conjugate of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns an
// expression representing the complex conjugate of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) conj( const Proxy<PT,RT>& proxy )
{
   using blaze::conj;

   return conj( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the transpose of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The transpose of the represented element.
//
// This function returns an expression representing the transpose of the element represented by
// the proxy.
*/
template< typename PT, typename RT >
inline decltype(auto) trans( const Proxy<PT,RT>& proxy )
{
   using blaze::trans;

   return trans( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the conjugate transpose of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The conjugate transpose of the represented element.
//
// This function returns an expression representing the conjugate transpose of the element
// represented by the proxy.
*/
template< typename PT, typename RT >
inline decltype(auto) ctrans( const Proxy<PT,RT>& proxy )
{
   using blaze::ctrans;

   return ctrans( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the real part of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The real part of the represented element.
//
// This function returns the real part of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the real part of each each element of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) real( const Proxy<PT,RT>& proxy )
{
   using blaze::real;

   return real( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the imaginary part of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The imaginary part of the represented element.
//
// This function returns the imaginary part of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the real part of each each element of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) imag( const Proxy<PT,RT>& proxy )
{
   using blaze::imag;

   return imag( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the square root of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The square root of the represented element.
//
// This function computes the square root of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the square roots of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) sqrt( const Proxy<PT,RT>& proxy )
{
   using blaze::sqrt;

   return sqrt( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse square root of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse square root of the represented element.
//
// This function computes the inverse square root of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the inverse square roots of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) invsqrt( const Proxy<PT,RT>& proxy )
{
   using blaze::invsqrt;

   return invsqrt( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the cubic root of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The cubic root of the represented element.
//
// This function computes the cubic root of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the cubic roots of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) cbrt( const Proxy<PT,RT>& proxy )
{
   using blaze::cbrt;

   return cbrt( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse cubic root of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse cubic root of the represented element.
//
// This function computes the inverse cubic root of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the inverse cubic roots of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) invcbrt( const Proxy<PT,RT>& proxy )
{
   using blaze::invcbrt;

   return invcbrt( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hypotenous of the two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The hypotenous of the given objects.
//
// This function computes the hypotenous of the elements represented by the two proxies \a lhs
// and \a rhs. In case the objects represent vector- or matrix-like data structures the function
// returns an expression representing the hypotenous of the elements of the vectors/matrices.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) hypot( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   using blaze::hypot;

   return hypot( (*lhs).get(), (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hypotenous of a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The hypotenous of the given objects.
//
// This function computes the hypotenous of the element represented by the proxy \a lhs and
// the object \a rhs. In case the objects represent vector- or matrix-like data structures
// the function returns an expression representing the hypotenous of the elements of the
// vectors/matrices.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) hypot( const Proxy<PT,RT>& lhs, const T& rhs )
{
   using blaze::hypot;

   return hypot( (*lhs).get(), rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the hypotenous of an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The hypotenous of the given objects.
//
// This function computes the hypotenous of the element represented by the object \a lhs and
// the proxy \a rhs. In case the objects represent vector- or matrix-like data structures
// the function returns an expression representing the hypotenous of the elements of the
// vectors/matrices.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) hypot( const T& lhs, const Proxy<PT,RT>& rhs )
{
   using blaze::hypot;

   return hypot( lhs, (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the exponential value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \param exp The exponent.
// \return The exponential value of the represented element.
//
// This function computes the exponential value of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the exponential value of the elements of the vector/matrix.
*/
template< typename PT, typename RT, typename ET >
inline decltype(auto) pow( const Proxy<PT,RT>& proxy, const ET& exp )
{
   using blaze::pow;

   return pow( (*proxy).get(), exp );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the square value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The squared value of the represented element.
//
// This function squares the element represented by the proxy. In case the proxy represents a
// vector- or matrix-like data structure the function returns an expression representing the
// squared value of the elements of the vector/matrix.
*/
template< typename PT, typename RT, typename ET >
inline decltype(auto) pow2( const Proxy<PT,RT>& proxy )
{
   using blaze::pow2;

   return pow2( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the cube value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The cubed value of the represented element.
//
// This function cubes the element represented by the proxy. In case the proxy represents a
// vector- or matrix-like data structure the function returns an expression representing the
// cubed value of the elements of the vector/matrix.
*/
template< typename PT, typename RT, typename ET >
inline decltype(auto) pow3( const Proxy<PT,RT>& proxy )
{
   using blaze::pow3;

   return pow3( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the quadruple value of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The quadrupled value of the represented element.
//
// This function quadruples the element represented by the proxy. In case the proxy represents
// a vector- or matrix-like data structure the function returns an expression representing the
// quadrupled value of the elements of the vector/matrix.
*/
template< typename PT, typename RT, typename ET >
inline decltype(auto) pow4( const Proxy<PT,RT>& proxy )
{
   using blaze::pow4;

   return pow4( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the base-e exponential of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The base-e exponential of the represented element.
//
// This function computes the base-e exponential of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the base-e exponentials of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) exp( const Proxy<PT,RT>& proxy )
{
   using blaze::exp;

   return exp( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the base-2 exponential of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The base-2 exponential of the represented element.
//
// This function computes the base-2 exponential of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the base-2 exponentials of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) exp2( const Proxy<PT,RT>& proxy )
{
   using blaze::exp2;

   return exp2( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the base-10 exponential of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The base-10 exponential of the represented element.
//
// This function computes the base-10 exponential of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the base-10 exponentials of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) exp10( const Proxy<PT,RT>& proxy )
{
   using blaze::exp10;

   return exp10( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the natural logarithm of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The natural logarithm of the represented element.
//
// This function computes the natural logarithm of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the natural logarithm of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) log( const Proxy<PT,RT>& proxy )
{
   using blaze::log;

   return log( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the binary logarithm of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The binary logarithm of the represented element.
//
// This function computes the binary logarithm of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the binary logarithm of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) log2( const Proxy<PT,RT>& proxy )
{
   using blaze::log2;

   return log2( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the common logarithm of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The common logarithm of the represented element.
//
// This function computes the common logarithm of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the common logarithm of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) log10( const Proxy<PT,RT>& proxy )
{
   using blaze::log10;

   return log10( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the natural logarithm of x+1 of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The natural logarithm of x+1 of the represented element.
//
// This function computes the natural logarithm of x+1 of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the natural logarithm of x+1 of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) log1p( const Proxy<PT,RT>& proxy )
{
   using blaze::log1p;

   return log1p( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the natural logarithm of the absolute value of the gamma function of the
//        represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The natural logarithm of the absolute value of the gamma function of the represented element.
//
// This function computes the natural logarithm of the absolute value of the gamma function of
// the element represented by the proxy. In case the proxy represents a vector- or matrix-like
// data structure the function returns an expression representing the natural logarithm of the
// elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) lgamma( const Proxy<PT,RT>& proxy )
{
   using blaze::lgamma;

   return lgamma( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the sine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The sine of the represented element.
//
// This function computes the sine of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the sines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) sin( const Proxy<PT,RT>& proxy )
{
   using blaze::sin;

   return sin( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse sine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse sine of the represented element.
//
// This function computes the inverse sine of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the inverse sines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) asin( const Proxy<PT,RT>& proxy )
{
   using blaze::asin;

   return asin( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the hyperbolic sine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The hyperbolic sine of the represented element.
//
// This function computes the hyperbolic sine of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the hyperbolic sines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) sinh( const Proxy<PT,RT>& proxy )
{
   using blaze::sinh;

   return sinh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse hyperbolic sine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse hyperbolic sine of the represented element.
//
// This function computes the inverse hyperbolic sine of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the inverse hyperbolic sines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) asinh( const Proxy<PT,RT>& proxy )
{
   using blaze::asinh;

   return asinh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the cosine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The cosine of the represented element.
//
// This function computes the cosine of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the cosines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) cos( const Proxy<PT,RT>& proxy )
{
   using blaze::cos;

   return cos( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse cosine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse cosine of the represented element.
//
// This function computes the inverse cosine of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the inverse cosines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) acos( const Proxy<PT,RT>& proxy )
{
   using blaze::acos;

   return acos( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the hyperbolic cosine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The hyperbolic cosine of the represented element.
//
// This function computes the hyperbolic cosine of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the hyperbolic cosines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) cosh( const Proxy<PT,RT>& proxy )
{
   using blaze::cosh;

   return cosh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse hyperbolic cosine of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse hyperbolic cosine of the represented element.
//
// This function computes the inverse hyperbolic cosine of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns an
// expression representing the inverse hyperbolic cosines of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) acosh( const Proxy<PT,RT>& proxy )
{
   using blaze::acosh;

   return acosh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the tangent of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The tangent of the represented element.
//
// This function computes the tangent of the element represented by the proxy. In case the
// proxy represents a vector- or matrix-like data structure the function returns an expression
// representing the tangents of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) tan( const Proxy<PT,RT>& proxy )
{
   using blaze::tan;

   return tan( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse tangent of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse tangent of the represented element.
//
// This function computes the inverse tangent of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the inverse tangents of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) atan( const Proxy<PT,RT>& proxy )
{
   using blaze::atan;

   return atan( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the hyperbolic tangent of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The hyperbolic tangent of the represented element.
//
// This function computes the hyperbolic tangent of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the hyperbolic tangents of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) tanh( const Proxy<PT,RT>& proxy )
{
   using blaze::tanh;

   return tanh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the inverse hyperbolic tangent of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The inverse hyperbolic tangent of the represented element.
//
// This function computes the inverse hyperbolic tangent of the element represented by the proxy.
// In case the proxy represents a vector- or matrix-like data structure the function returns an
// expression representing the inverse hyperbolic tangents of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) atanh( const Proxy<PT,RT>& proxy )
{
   using blaze::atanh;

   return atanh( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the multi-valued inverse tangent of two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return The multi-valued inverse tangent of the given objects.
//
// This function computes the multi-valued inverse tangent of the elements represented by the two
// proxies \a lhs and \a rhs. In case the objects represent vector- or matrix-like data structures
// the function returns an expression representing the multi-valued inverse tangent of the elements
// of the vectors/matrices.
*/
template< typename PT1, typename RT1, typename PT2, typename RT2 >
inline decltype(auto) atan2( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return atan2( (*lhs).get(), (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the multi-valued inverse tangent of a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return The multi-valued inverse tangent of the given objects.
//
// This function computes the multi-valued inverse tangent of the element represented by the
// proxy \a lhs and the object \a rhs. In case the objects represent vector- or matrix-like data
// structures the function returns an expression representing the multi-valued inverse tangent of
// the elements of the vectors/matrices.
*/
template< typename PT, typename RT, typename T, typename >
inline decltype(auto) atan2( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return atan2( (*lhs).get(), rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computes the multi-valued inverse tangent of an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return The multi-valued inverse tangent of the given objects.
//
// This function computes the multi-valued inverse tangent of the element represented by the
// object \a lhs and the proxy \a rhs. In case the objects represent vector- or matrix-like data
// structures the function returns an expression representing the multi-valued inverse tangent of
// the elements of the vectors/matrices.
*/
template< typename T, typename PT, typename RT, typename >
inline decltype(auto) atan2( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return atan2( lhs, (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the error function of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The error function of the represented element.
//
// This function computes the error function of the element represented by the proxy. In
// case the proxy represents a vector- or matrix-like data structure the function returns
// an expression representing the error functions of the elements of the vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) erf( const Proxy<PT,RT>& proxy )
{
   using blaze::erf;

   return erf( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Computing the complementary error function of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return The complementary error function of the represented element.
//
// This function computes the complementary error function of the element represented by the
// proxy. In case the proxy represents a vector- or matrix-like data structure the function
// returns an expression representing the complementary error functions of the elements of the
// vector/matrix.
*/
template< typename PT, typename RT >
inline decltype(auto) erfc( const Proxy<PT,RT>& proxy )
{
   using blaze::erfc;

   return erfc( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place transpose of the represented matrix element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::logic_error Matrix cannot be transposed.
//
// This function transposes the represented matrix in-place. The transpose operation fails if ...
//
//  - ... the represented matrix has a fixed size and is non-square;
//  - ... the represented matrix is a triangular matrix.
//
// In all failure cases a \a std::logic_error exception is thrown. Additionally, in case the
// represented matrix cannot be modified, a \a std::invalid_argument exception is thrown.
*/
template< typename PT, typename RT >
inline void transpose( const Proxy<PT,RT>& proxy )
{
   if( (*proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   transpose( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place conjugate transpose of the represented matrix element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::logic_error Matrix cannot be transposed.
//
// This function transposes the represented matrix in-place. The transpose operation fails if ...
//
//  - ... the represented matrix has a fixed size and is non-square;
//  - ... the represented matrix is a triangular matrix.
//
// In all failure cases a \a std::logic_error exception is thrown. Additionally, in case the
// represented matrix cannot be modified, a \a std::invalid_argument exception is thrown.
*/
template< typename PT, typename RT >
inline void ctranspose( const Proxy<PT,RT>& proxy )
{
   if( (*proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   ctranspose( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the represented scalar or dense matrix element. The inversion fails if
// the represented element is a dense matrix, which ...
//
//  - ... is not a square matrix;
//  - ... is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \a std::invalid_argument exception is thrown. Additionally, in case the
// represented scalar or matrix cannot be modified, a \a std::invalid_argument exception is thrown.
//
// \note In case the represented element is a dense matrix, this function does not provide any
// exception safety guarantee, i.e. in case an exception is thrown the matrix may already have
// been modified.
//
// \note In case the represented element is a dense matrix, this function can only be used if the
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error will
// be created.
*/
template< typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy )
{
   if( (*proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   invert( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief In-place inversion of the represented element.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return void
// \exception std::invalid_argument Invalid access to restricted element.
// \exception std::invalid_argument Inversion of singular matrix failed.
// \exception std::invalid_argument Invalid non-square matrix provided.
//
// This function inverts the represented dense matrix element by means of the specified matrix
// inversion algorithm \c IF:

   \code
   invert<byLU>( A );    // Inversion of a general matrix
   invert<byLDLT>( A );  // Inversion of a symmetric indefinite matrix
   invert<byLDLH>( A );  // Inversion of a Hermitian indefinite matrix
   invert<byLLH>( A );   // Inversion of a Hermitian positive definite matrix
   \endcode

// The inversion fails if the represented dense matrix element ...
//
//  - ... is not a square matrix;
//  - ... is singular and not invertible.
//
// In all failure cases either a compilation error is created if the failure can be predicted at
// compile time or a \a std::invalid_argument exception is thrown. Additionally, in case the
// represented scalar or matrix cannot be modified, a \a std::invalid_argument exception is thrown.
//
// \note In case the represented element is a dense matrix, this function does not provide any
// exception safety guarantee, i.e. in case an exception is thrown the matrix may already have
// been modified.
//
// \note In case the represented element is a dense matrix, this function can only be used if the
// fitting LAPACK library is available and linked to the executable. Otherwise a linker error will
// be created.
*/
template< InversionFlag IF, typename PT, typename RT >
inline void invert( const Proxy<PT,RT>& proxy )
{
   if( (*proxy).isRestricted() ) {
      BLAZE_THROW_INVALID_ARGUMENT( "Invalid access to restricted element" );
   }

   invert<IF>( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the element represents a real number.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the element represents a real number, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the a real
// number. In case the element is of built-in type, the function returns \a true. In case
// the element is of complex type, the function returns \a true if the imaginary part is
// equal to 0. Otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename PT, typename RT >
inline bool isReal( const Proxy<PT,RT>& proxy )
{
   using blaze::isReal;

   return isReal<RF>( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 0.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is 0, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the numeric
// value 0. In case it is 0, the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename PT, typename RT >
inline bool isZero( const Proxy<PT,RT>& proxy )
{
   using blaze::isZero;

   return isZero<RF>( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is 1.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is 1, \a false otherwise.
//
// This function checks whether the element represented by the proxy represents the numeric
// value 1. In case it is 1, the function returns \a true, otherwise it returns \a false.
*/
template< RelaxationFlag RF, typename PT, typename RT >
inline bool isOne( const Proxy<PT,RT>& proxy )
{
   using blaze::isOne;

   return isOne<RF>( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is not-a-number.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is not-a-number, \a false otherwise.
//
// This function checks whether the element represented by the proxy is not-a-number (NaN).
// In case it is not-a-number, the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isnan( const Proxy<PT,RT>& proxy )
{
   using blaze::isnan;

   return isnan( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is infinite.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is infinite, \a false otherwise.
//
// This function checks whether the element represented by the proxy is infinite (inf).
// In case it is infinite, the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isinf( const Proxy<PT,RT>& proxy )
{
   using blaze::isinf;

   return isinf( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is finite.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is finite, \a false otherwise.
//
// This function checks whether the element represented by the proxy is finite (i.e.
// normal, subnormal or zero elements, but not infinite or NaN). In case it is finite,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline bool isfinite( const Proxy<PT,RT>& proxy )
{
   using blaze::isfinite;

   return isfinite( (*proxy).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between two Proxy objects.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side Proxy object.
// \return \a true if both referenced values are equal, \a false if they are not.
*/
template< RelaxationFlag RF, typename PT1, typename RT1, typename PT2, typename RT2 >
inline bool equal( const Proxy<PT1,RT1>& lhs, const Proxy<PT2,RT2>& rhs )
{
   return equal<RF>( (*lhs).get(), (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between a Proxy object and an object of different type.
// \ingroup math
//
// \param lhs The left-hand side Proxy object.
// \param rhs The right-hand side object of other type.
// \return \a true if the referenced value and the other object are equal, \a false if they are not.
*/
template< RelaxationFlag RF, typename PT, typename RT, typename T, typename >
inline bool equal( const Proxy<PT,RT>& lhs, const T& rhs )
{
   return equal<RF>( (*lhs).get(), rhs );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Equality comparison between an object of different type and a Proxy object.
// \ingroup math
//
// \param lhs The left-hand side object of other type.
// \param rhs The right-hand side Proxy object.
// \return \a true if the other object and the referenced value are equal, \a false if they are not.
*/
template< RelaxationFlag RF, typename T, typename PT, typename RT, typename >
inline bool equal( const T& lhs, const Proxy<PT,RT>& rhs )
{
   return equal<RF>( lhs, (*rhs).get() );
}
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Returns whether the represented element is finite.
// \ingroup math
//
// \param proxy The given proxy instance.
// \return \a true in case the represented element is finite, \a false otherwise.
//
// This function checks whether the element represented by the proxy is finite (i.e.
// normal, subnormal or zero elements, but not infinite or NaN). In case it is finite,
// the function returns \a true, otherwise it returns \a false.
*/
template< typename PT, typename RT >
inline decltype(auto) evaluate( const Proxy<PT,RT>& proxy )
{
   using blaze::evaluate;

   return evaluate( (*proxy).get() );
}
//*************************************************************************************************

} // namespace blaze

#endif
