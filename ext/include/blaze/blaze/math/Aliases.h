//=================================================================================================
/*!
//  \file blaze/math/Aliases.h
//  \brief Header file for auxiliary alias declarations
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

#ifndef _BLAZE_MATH_ALIASES_H_
#define _BLAZE_MATH_ALIASES_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/Types.h>


namespace blaze {

//=================================================================================================
//
//  ALIAS DECLARATION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Alias declaration for nested \c AllocatorType type definitions.
// \ingroup aliases
//
// The AllocatorType_t alias declaration provides a convenient shortcut to access the nested
// \a AllocatorType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::AllocatorType;
   using Type2 = AllocatorType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using AllocatorType_t = typename T::AllocatorType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c BaseType type definitions.
// \ingroup aliases
//
// The BaseType_t alias declaration provides a convenient shortcut to access the nested
// \a BaseType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::BaseType;
   using Type2 = BaseType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using BaseType_t = typename T::BaseType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c CompositeType type definitions.
// \ingroup aliases
//
// The CompositeType_t alias declaration provides a convenient shortcut to access the nested
// \a CompositeType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::CompositeType;
   using Type2 = CompositeType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using CompositeType_t = typename T::CompositeType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ConstIterator type definitions.
// \ingroup aliases
//
// The ConstIterator_t alias declaration provides a convenient shortcut to access the nested
// \a ConstIterator type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ConstIterator;
   using Type2 = ConstIterator_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ConstIterator_t = typename T::ConstIterator;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ConstPointer type definitions.
// \ingroup aliases
//
// The ConstPointer_t alias declaration provides a convenient shortcut to access the nested
// \a ConstPointer type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ConstPointer;
   using Type2 = ConstPointer_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ConstPointer_t = typename T::ConstPointer;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ConstReference type definitions.
// \ingroup aliases
//
// The ConstReference_t alias declaration provides a convenient shortcut to access the nested
// \a ConstReference type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ConstReference;
   using Type2 = ConstReference_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ConstReference_t = typename T::ConstReference;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ElementType type definitions.
// \ingroup aliases
//
// The ElementType_t alias declaration provides a convenient shortcut to access the nested
// \a ElementType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ElementType;
   using Type2 = ElementType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ElementType_t = typename T::ElementType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Iterator type definitions.
// \ingroup aliases
//
// The Iterator_t alias declaration provides a convenient shortcut to access the nested
// \a Iterator type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::Iterator;
   using Type2 = Iterator_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using Iterator_t = typename T::Iterator;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c LeftOperand type definitions.
// \ingroup aliases
//
// The LeftOperand_t alias declaration provides a convenient shortcut to access the nested
// \a LeftOperand type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::LeftOperand;
   using Type2 = LeftOperand_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using LeftOperand_t = typename T::LeftOperand;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c MatrixType type definitions.
// \ingroup aliases
//
// The MatrixType_t alias declaration provides a convenient shortcut to access the nested
// \a MatrixType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::MatrixType;
   using Type2 = MatrixType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using MatrixType_t = typename T::MatrixType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Operand type definitions.
// \ingroup aliases
//
// The Operand_t alias declaration provides a convenient shortcut to access the nested \a Operand
// type definition of the given type \a T. The following code example shows both ways to access
// the nested type definition:

   \code
   using Type1 = typename T::Operand;
   using Type2 = Operand_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using Operand_t = typename T::Operand;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Operation type definitions.
// \ingroup aliases
//
// The Operation_t alias declaration provides a convenient shortcut to access the nested
// \a Operation type definition of the given type \a T. The following code example shows both
// ways to access the nested type definition:

   \code
   using Type1 = typename T::Operation;
   using Type2 = Operation_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using Operation_t = typename T::Operation;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c OppositeType type definitions.
// \ingroup aliases
//
// The OppositeType_t alias declaration provides a convenient shortcut to access the nested
// \a OppositeType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::OppositeType;
   using Type2 = OppositeType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using OppositeType_t = typename T::OppositeType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Pointer type definitions.
// \ingroup aliases
//
// The Pointer_t alias declaration provides a convenient shortcut to access the nested
// \a Pointer type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::Pointer;
   using Type2 = Pointer_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using Pointer_t = typename T::Pointer;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Rebind class templates.
// \ingroup aliases
//
// The Rebind_t alias declaration provides a convenient shortcut to access the nested \a Rebind
// class template of the given type \a T1. The following code example shows both ways to access
// the nested class template:

   \code
   using Type1 = typename T1::template Rebind<T2>::Other;
   using Type2 = Rebind_t<T1,T2>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T1, typename T2 >
using Rebind_t = typename T1::template Rebind<T2>::Other;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c rebind class templates.
// \ingroup aliases
//
// The rebind_t alias declaration provides a convenient shortcut to access the nested \a rebind
// class template of the given type \a T1. The following code example shows both ways to access
// the nested class template:

   \code
   using Type1 = typename T1::template rebind<T2>::Other;
   using Type2 = rebind_t<T1,T2>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T1, typename T2 >
using rebind_t = typename T1::template rebind<T2>::other;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Reference type definitions.
// \ingroup aliases
//
// The Reference_t alias declaration provides a convenient shortcut to access the nested
// \a Reference type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::Reference;
   using Type2 = Reference_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using Reference_t = typename T::Reference;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c RepresentedType type definitions.
// \ingroup aliases
//
// The RepresentedType_t alias declaration provides a convenient shortcut to access the nested
// \a RepresentedType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::RepresentedType;
   using Type2 = RepresentedType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using RepresentedType_t = typename T::RepresentedType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c Resize class templates.
// \ingroup aliases
//
// The Resize_t alias declaration provides a convenient shortcut to access the nested \a Resize
// class template of the given type \a T1. The following code example shows both ways to access
// the nested class template:

   \code
   using Type1 = typename T::template Resize<N>::Other;
   using Type2 = Resize_t<T,N>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T, size_t... Ns >
using Resize_t = typename T::template Resize<Ns...>::Other;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ResultType type definitions.
// \ingroup aliases
//
// The ResultType_t alias declaration provides a convenient shortcut to access the nested
// \a ResultType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ResultType;
   using Type2 = ResultType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ResultType_t = typename T::ResultType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ReturnType type definitions.
// \ingroup aliases
//
// The ReturnType_t alias declaration provides a convenient shortcut to access the nested
// \a ReturnType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ReturnType;
   using Type2 = ReturnType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ReturnType_t = typename T::ReturnType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c RightOperand type definitions.
// \ingroup aliases
//
// The RightOperand_t alias declaration provides a convenient shortcut to access the nested
// \a RightOperand type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::RightOperand;
   using Type2 = RightOperand_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using RightOperand_t = typename T::RightOperand;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c SIMDType type definitions.
// \ingroup aliases
//
// The SIMDType_t alias declaration provides a convenient shortcut to access the nested
// \a SIMDType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::SIMDType;
   using Type2 = SIMDType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using SIMDType_t = typename T::SIMDType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c TagType type definitions.
// \ingroup aliases
//
// The TagType_t alias declaration provides a convenient shortcut to access the nested
// \a TagType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::TagType;
   using Type2 = TagType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using TagType_t = typename ResultType_t<T>::TagType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c TransposeType type definitions.
// \ingroup aliases
//
// The TransposeType_t alias declaration provides a convenient shortcut to access the nested
// \a TransposeType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::TransposeType;
   using Type2 = TransposeType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using TransposeType_t = typename T::TransposeType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ValueType type definitions.
// \ingroup aliases
//
// The ValueType_t alias declaration provides a convenient shortcut to access the nested
// \a ValueType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ValueType;
   using Type2 = ValueType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ValueType_t = typename T::ValueType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c VectorType type definitions.
// \ingroup aliases
//
// The VectorType_t alias declaration provides a convenient shortcut to access the nested
// \a VectorType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::VectorType;
   using Type2 = VectorType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using VectorType_t = typename T::VectorType;
//*************************************************************************************************


//*************************************************************************************************
/*!\brief Alias declaration for nested \c ViewedType type definitions.
// \ingroup aliases
//
// The ViewedType_t alias declaration provides a convenient shortcut to access the nested
// \a ViewedType type definition of the given type \a T. The following code example shows
// both ways to access the nested type definition:

   \code
   using Type1 = typename T::ViewedType;
   using Type2 = ViewedType_t<T>;

   BLAZE_CONSTRAINT_MUST_BE_STRICTLY_SAME_TYPE( Type1, Type2 );
   \endcode
*/
template< typename T >
using ViewedType_t = typename T::ViewedType;
//*************************************************************************************************

} // namespace blaze

#endif
