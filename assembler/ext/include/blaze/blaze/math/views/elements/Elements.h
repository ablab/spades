//=================================================================================================
/*!
//  \file blaze/math/views/elements/Elements.h
//  \brief Elements documentation
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

#ifndef _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTS_H_
#define _BLAZE_MATH_VIEWS_ELEMENTS_ELEMENTS_H_


//=================================================================================================
//
//  DOXYGEN DOCUMENTATION
//
//=================================================================================================

//*************************************************************************************************
/*!\defgroup elements Elements
// \ingroup views
//
// Element selections provide views on arbitrary compositions of elements of dense and sparse
// vectors. These views act as a reference to the selected elements and represent them as another
// dense or sparse vector. This reference is valid and can be used in every way any other dense
// or sparse vector can be used as long as the vector containing the elements is not resized or
// entirely destroyed. The element selection also acts as an alias to the vector elements in the
// specified range: Changes made to the elements (e.g. modifying values, inserting or erasing
// elements) are immediately visible in the vector and changes made via the vector are immediately
// visible in the elements.
//
//
// \n \section elements_setup Setup of Element Selections
//
// An element selection can be created very conveniently via the \c elements() function. It can
// be included via the header files

   \code
   #include <blaze/Blaze.h>
   // or
   #include <blaze/Math.h>
   // or
   #include <blaze/math/Elements.h>
   \endcode

// and forward declared via the header file

   \code
   #include <blaze/Forward.h>
   \endcode

// The indices of the elements to be selected can be specified either at compile time or at runtime
// (by means of an initializer list, array or vector):

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   // ... Resizing and initialization

   // Selecting the elements 4, 6, 8, and 10 (compile time arguments)
   auto e1 = elements<4UL,6UL,8UL,10UL>( x );

   // Selecting the elements 3, 2, and 1 (runtime arguments via an initializer list)
   const std::initializer_list<size_t> list{ 3UL, 2UL, 1UL };
   auto e2 = elements( x, { 3UL, 2UL, 1UL } );
   auto e3 = elements( x, list );

   // Selecting the elements 1, 2, 3, 3, 2, and 1 (runtime arguments via a std::array)
   const std::array<size_t> array{ 1UL, 2UL, 3UL, 3UL, 2UL, 1UL };
   auto e4 = elements( x, array );
   auto e5 = elements( x, array.data(), array.size() );

   // Selecting the element 4 fives times (runtime arguments via a std::vector)
   const std::vector<size_t> vector{ 4UL, 4UL, 4UL, 4UL, 4UL };
   auto e6 = elements( x, vector );
   auto e7 = elements( x, vector.data(), vector.size() );
   \endcode

// Note that it is possible to alias the elements of the underlying vector in any order. Also note
// that it is possible to use the same index multiple times.
//
// Alternatively it is possible to pass a callable such as a lambda or functor that produces the
// indices:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

   // Selecting all even elements of the vector, i.e. selecting (1,3,5,7,9)
   auto e1 = elements( x, []( size_t i ){ return i*2UL; }, 5UL );

   // Selecting all odd elements of the vector, i.e. selecting (2,4,6,8)
   auto e2 = elements( x, []( size_t i ){ return i*2UL+1UL; }, 4UL );

   // Reversing the elements of the vector, i.e. selecting (9,8,7,6,5,4,3,2,1)
   auto e3 = elements( v, [max=v.size()-1UL]( size_t i ){ return max-i; }, 9UL );
   \endcode

// The \c elements() function returns an expression representing the view on the selected elements.
// The type of this expression depends on the given arguments, primarily the type of the vector
// and the compile time arguments. If the type is required, it can be determined via \c decltype
// specifier:

   \code
   using VectorType = blaze::DynamicVector<int>;
   using ElementsType = decltype( blaze::elements<4UL,12UL>( std::declval<VectorType>() ) );
   \endcode

// The resulting view can be treated as any other dense or sparse vector, i.e. it can be assigned
// to, it can be copied from, and it can be used in arithmetic operations. An element selection
// created from a row vector can be used as any other row vector, an element selection created
// from a column vector can be used as any other column vector. The view can also be used on both
// sides of an assignment: It can either be used as an alias to grant write access to specific
// elements of a vector primitive on the left-hand side of an assignment or to grant read-access
// to specific elements of a vector primitive or expression on the right-hand side of an assignment.
// The following example demonstrates this in detail:

   \code
   blaze::DynamicVector<double,blaze::rowVector> x;
   blaze::CompressedVector<double,blaze::rowVector> y;
   blaze::DynamicMatrix<double,blaze::rowMajor> A;
   // ... Resizing and initialization

   // Selecting the elements 1, 3, 5, and 7
   auto e = elements( x, { 1UL, 3UL, 5UL, 7UL } );

   // Setting the elements 1, 3, 5, and 7 of x to the 2nd row of matrix A
   e = row( A, 2UL );

   // Setting the elements 2, 4, 6, and 8 of x to y
   elements( x, { 2UL, 4UL, 6UL, 8UL } ) = y;

   // Setting the 3rd row of A to the elements 5, 4, 3, and 2 of x
   row( A, 3UL ) = elements( x, { 5UL, 4UL, 3UL, 2UL } );

   // Rotating the result of the addition between y and the 1st row of A
   x = elements( y + row( A, 1UL ), { 2UL, 3UL, 0UL, 1UL } )
   \endcode

// Please note that using an element selection, which refers to an index multiple times, on the
// left-hand side of an assignment leads to undefined behavior:

   \code
   blaze::DynamicVector<int,blaze::rowVector> a{ 1, 2, 3 };
   blaze::DynamicVector<int,blaze::rowVector> b{ 1, 2, 3, 4 };

   auto e = elements( a, { 1, 1, 1, 1 } );  // Selecting the element 1 four times
   e = b;  // Undefined behavior
   \endcode

// In this example both vectors have the same size, which results in a correct vector assignment,
// but the final value of the element at index 1 is unspecified.
//
//
// \n \section elements_element_access Element access
//
// The elements of an element selection can be directly accessed via the subscript operator:

   \code
   blaze::DynamicVector<double,blaze::rowVector> v;
   // ... Resizing and initialization

   // Selecting the elements 2, 4, 6, and 8
   auto e = elements( v, { 2UL, 4UL, 6UL, 8UL } );

   // Setting the 1st element of the element selection, which corresponds to
   // the element at index 4 in vector v
   e[1] = 2.0;
   \endcode

// The numbering of the selected elements is

                             \f[\left(\begin{array}{*{5}{c}}
                             0 & 1 & 2 & \cdots & N-1 \\
                             \end{array}\right),\f]

// where N is the number of selected elements. Alternatively, the elements of an element selection
// can be traversed via iterators. Just as with vectors, in case of non-const element selections,
// \c begin() and \c end() return an iterator, which allows to manipulate the elements, in case of
// constant element selections an iterator to immutable elements is returned:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating an element selection including specific elements of dense vector v
   auto e = elements( v, { 0UL, 3UL, 6UL, 9UL, 12UL } );

   // Traversing the elements via iterators to non-const elements
   for( auto it=e.begin(); it!=e.end(); ++it ) {
      *it = ...;  // OK: Write access to the dense vector value.
      ... = *it;  // OK: Read access to the dense vector value.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=e.cbegin(); it!=e.cend(); ++it ) {
      *it = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = *it;  // OK: Read access to the dense vector value.
   }
   \endcode

   \code
   blaze::CompressedVector<int,blaze::rowVector> v( 256UL );
   // ... Resizing and initialization

   // Creating an element selection including specific elements of sparse vector v
   auto e = elements( v, { 0UL, 3UL, 6UL, 9UL, 12UL } );

   // Traversing the elements via iterators to non-const elements
   for( auto it=e.begin(); it!=e.end(); ++it ) {
      it->value() = ...;  // OK: Write access to the value of the non-zero element.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }

   // Traversing the elements via iterators to const elements
   for( auto it=e.cbegin(); it!=e.cend(); ++it ) {
      it->value() = ...;  // Compilation error: Assignment to the value via iterator-to-const is invalid.
      ... = it->value();  // OK: Read access to the value of the non-zero element.
      it->index() = ...;  // Compilation error: The index of a non-zero element cannot be changed.
      ... = it->index();  // OK: Read access to the index of the sparse element.
   }
   \endcode

// \n \section elements_element_insertion Element Insertion
//
// Inserting/accessing elements in a sparse element selection can be done by several alternative
// functions. The following example demonstrates all options:

   \code
   blaze::CompressedVector<double,blaze::rowVector> v( 256UL );  // Non-initialized vector of size 256

   std::vector<size_t> indices;
   // ... Selecting indices of the sparse vector

   auto e = elements( v, indices );

   // The subscript operator provides access to the selected elements of the sparse vector,
   // including the zero elements. In case the subscript operator is used to access an element
   // that is currently not stored in the sparse vector, the element is inserted.
   e[42] = 2.0;

   // The second operation for inserting elements via the element selection is the set() function.
   // In case the element is not contained in the vector it is inserted into the vector, if it is
   // already contained in the vector its value is modified.
   e.set( 45UL, -1.2 );

   // An alternative for inserting elements into the vector is the insert() function. However, it
   // inserts the element only in case the element is not already contained in the vector.
   e.insert( 50UL, 3.7 );

   // Just as in case of vectors, elements can also be inserted via the append() function. In case
   // of element selections, append() also requires that the appended element's index is strictly
   // larger than the currently largest non-zero index of the selection and that the selections's
   // capacity is large enough to hold the new element. Note however that due to the nature of an
   // element selection, which is an alias to arbitrary elements of a sparse vector, the append()
   // function does not work as efficiently for an element selection as it does for a vector.
   e.reserve( 10UL );
   e.append( 51UL, -2.1 );
   \endcode

// \n \section elements_common_operations Common Operations
//
// An element selection can be used like any other dense or sparse vector. For instance, the
// number of selected elements can be obtained via the \c size() function, the current capacity
// via the \c capacity() function, and the number of non-zero elements via the \c nonZeros()
// function. However, since element selections are references to a specific range of a vector,
// several operations are not possible, such as resizing and swapping. The following example
// shows this by means of an element selection on a dense vector:

   \code
   blaze::DynamicVector<int,blaze::rowVector> v( 42UL );
   // ... Resizing and initialization

   // Selecting the elements 5 and 10
   auto e = elements( v, { 5UL, 10UL } );

   e.size();          // Returns the number of elements in the element selection
   e.capacity();      // Returns the capacity of the element selection
   e.nonZeros();      // Returns the number of non-zero elements contained in the element selection

   e.resize( 84UL );  // Compilation error: Cannot resize an element selection

   auto e2 = elements( v, { 15UL, 10UL } );
   swap( e, e2 );   // Compilation error: Swap operation not allowed
   \endcode

// \n \section elements_arithmetic_operations Arithmetic Operations
//
// Both dense and sparse element selections can be used in all arithmetic operations that any other
// dense or sparse vector can be used in. The following example gives an impression of the use of
// dense element selections within arithmetic operations. All operations (addition, subtraction,
// multiplication, scaling, ...) can be performed on all possible combinations of dense and sparse
// element selections with fitting element types:

   \code
   blaze::DynamicVector<double,blaze::rowVector> d1, d2, d3;
   blaze::CompressedVector<double,blaze::rowVector> s1, s2;
   blaze::DynamicMatrix<double,blaze::rowMajor> A;

   // ... Resizing and initialization

   std::initializer_list<size_t> indices1{ 0UL, 3UL, 6UL,  9UL, 12UL, 15UL, 18UL, 21UL };
   std::initializer_list<size_t> indices2{ 1UL, 4UL, 7UL, 10UL, 13UL, 16UL, 19UL, 22UL };
   std::initializer_list<size_t> indices3{ 2UL, 5UL, 8UL, 11UL, 14UL, 17UL, 20UL, 23UL };

   auto e( elements( d1, indices1 ) );  // Selecting the every third element of d1 in the range [0..21]

   e = d2;                         // Dense vector assignment to the selected elements
   elements( d1, indices2 ) = s1;  // Sparse vector assignment to the selected elements

   d3 = e + d2;                         // Dense vector/dense vector addition
   s2 = s1 + elements( d1, indices2 );  // Sparse vector/dense vector addition
   d2 = e * elements( d1, indices3 );   // Component-wise vector multiplication

   elements( d1, indices2 ) *= 2.0;      // In-place scaling of the second selection of elements
   d2 = elements( d1, indices3 ) * 2.0;  // Scaling of the elements in the third selection of elements
   d2 = 2.0 * elements( d1, indices3 );  // Scaling of the elements in the third selection of elements

   elements( d1, indices1 ) += d2;  // Addition assignment
   elements( d1, indices2 ) -= s2;  // Subtraction assignment
   elements( d1, indices3 ) *= e;   // Multiplication assignment

   double scalar = elements( d1, indices2 ) * trans( s1 );  // Scalar/dot/inner product between two vectors

   A = trans( s1 ) * elements( d1, { 3UL, 6UL } );  // Outer product between two vectors
   \endcode
*/
//*************************************************************************************************

#endif
