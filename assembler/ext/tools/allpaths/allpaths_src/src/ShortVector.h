///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#ifndef SHORTVECTOR
#define SHORTVECTOR

#include "system/Assert.h"
#include "system/Types.h"

// ================================================================================
//
// A shortvector holds a list of up to 255 things of any type T.  
//
// Construction of shortvector (by example):
//
// *** shortvector(5) gives you a length 5 vector filled with junk
//
// *** shortvector(1, t) gives you a vector with one object in it, namely t,
//     assuming that t is of type T.
//
// Accessing the elements of a vector v (by example):
//
// *** v(5) refers to the fifth (really sixth) element of v,
//     either for READ or WRITE access
//
// Modifying vector (EXPENSIVE, in all cases):
//
// *** Append(t) puts t on the end
// *** Prepend(t) inserts t at the beginning
// *** Setsize(n) changes the length to n, destroying the vector contents
//
// ================================================================================

template<class T> class shortvector {

   public:

     T* x;                     // the elements
     unsigned short length;    // number of elements

     shortvector( ) :
       x(0),
       length(0)
     {  }

     shortvector(unsigned short n)  // Construct shortvector with n default elements.
     {    x = new T[n];
          length = n;    
          AssertLe( length, 255 );    }

     shortvector(int n, const T& t) // n should be 1
     {    x = new T[1];
          x[0] = t;
          length = 1;    }

     shortvector(int n, const T& t1, const T& t2) // n should be 2
     {    x = new T[2];
          x[0] = t1;
          x[1] = t2;
          length = 2;    }

     shortvector(const shortvector& v) 
     {    if ( ! v.x ) x = 0;
          else
	  {    length = v.length;
	       x = new T[length];
	       for ( unsigned short i = 0; i < length; i++ )
	            x[i] = v(i);    }    }

     ~shortvector( ) { if ( x ) delete [ ] x; }

     void Setsize(unsigned short n) // resize, destroying contents
     {    if ( x && length == n ) return;
          if ( x ) delete [ ] x;
          x = new T[n];
          length = n;    }
     void resize(unsigned short n)
     {    T* x_new = new T[n];
          for ( unsigned short i = 0; i < min( n, length ); i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          length = n;    }
     void Append( const T& t )
     {    AssertLt( length, 255 );
          T* x_new = new T[ length + 1 ];
          for ( unsigned short i = 0; i < length; i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          ++length;
          (*this)(length-1) = t;    }
     void Prepend( const T& t )
     {    AssertLt( length, 255 );
          ++length;
          T* x_new = new T[length];
          for ( unsigned short i = 1; i < length; i++ )
               x_new[i] = x[i-1];
          delete [ ] x;
          x = x_new;    
          (*this)(0) = t;    }
     shortvector& operator=(const shortvector& v)
     {    if ( ! v.x )
          {    if ( x ) delete [ ] x;
               x = 0;    }
          else
          {    if ( ! x )
               {    length = v.length;
                    x = new T[length];    }
               else if ( length != v.length )
               {    delete [ ] x;
                    length = v.length;
                    x = new T[length];    }
               for ( unsigned short i = 0; i < v.length; i++ )
                    (*this).x[i] = v(i);    }
           return *this;    }
     T& operator( )(unsigned short i) const 
     {    AssertLt( i, length );
          return x[i];    }

};

template<class T> class avector {

public:

  T* x;                   ///< the elements
  unsigned int length;    ///< number of elements

  avector( ) :
    x(0),
    length(0)
  {  }
  
  /// Construct avector with n default elements.
  explicit avector(unsigned int n)  
  {    x = new T[n];
  length = n;    }

  avector(int n, const T & t): x(new T[n]), length(n) {    
    uninitialized_fill_n(x, n, t);
  }

  avector(const avector& v)
  {    if ( ! v.x )
    {    x = 0;
    length = 0;   }
  else
    {    length = v.length;
    x = new T[length];
    for ( unsigned int i = 0; i < length; i++ )
      x[i] = v(i);    }    }

  template<class ForwardIterator>
  avector(ForwardIterator begin, ForwardIterator end) {
    length = std::distance(begin, end);
    x = new T[length];
    for (unsigned int i=0; i != length; ++i) {
      x[i] = *(begin+i);
    }
  }

  ~avector( ) { if ( x ) delete [ ] x; }

     void Setsize(unsigned int n) // resize, destroying contents
     {    if ( x && length == n ) return;
          if ( x ) delete [ ] x;
          x = new T[n];
          length = n;    }
     void resize(unsigned int n)
     {    T* x_new = new T[n];
          for ( unsigned int i = 0; i < min( n, length ); i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          length = n;    }
     // Below version of Append creates a zero-length vector that
     // flags a warning in gcc 4.2 (and rightly so).
     /* void Append( const T& t )
     {    if ( x == 0 )
          {    x = new T[0];
               length = 0;    }
          T* x_new = new T[ length + 1 ];
          for ( unsigned int i = 0; i < length; i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          ++length;
          (*this)(length-1) = t;    } */
     void Append( const T& t )
     {
        if ( x == 0 ) length = 0;
        T* x_new = new T[ length + 1 ];
        for ( unsigned int i = 0; i < length; i++ )
             x_new[i] = x[i];  // Will never hit when length = 0.
        if ( x != 0 ) delete [ ] x;
	x = x_new;
	++length;
        (*this)(length-1) = t;
     }
     void Prepend( const T& t )
     {    ++length;
          T* x_new = new T[length];
          for ( unsigned int i = 1; i < length; i++ )
               x_new[i] = x[i-1];
          delete [ ] x;
          x = x_new;    
          (*this)(0) = t;    }
     avector& operator=(const avector& v)
     {    if ( !v.x )
          {    if ( x ) delete [ ] x;
               x = 0;    }
          else
          {    if ( !x )
               {    length = v.length;
                    x = new T[length];    }
               else if ( length != v.length )
               {    delete [ ] x;
                    length = v.length;
                    x = new T[length];    }
               for ( unsigned int i = 0; i < v.length; i++ )
                    (*this).x[i] = v(i);    }
           return *this;    }
     T& operator( )(unsigned int i) const 
     {    AssertLt( i, length );
          return x[i];    }
     friend Bool operator<( const avector& v1, const avector& v2 )
     {    if ( v1.length < v2.length ) return True;
          if ( v1.length > v2.length ) return False;
          for ( unsigned int j = 0; j < v1.length; j++ )
          {    if ( v1(j) < v2(j) ) return True;
               if ( v1(j) > v2(j) ) return False;    }
          return False;    }
     friend Bool operator>( const avector& v1, const avector& v2 )
     {     return v2 < v1;    }

};

template<class T> class bvector {

   public:

     T* x;                   // the elements
     longlong length;        // number of elements

     bvector( ) :
       x(0),
       length(0)
     {  }

     bvector(longlong n)  // Construct avector with n default elements.
     {    x = new T[n];
          length = n;    }

     bvector(int n, const T& t) // n should be 1
     {    x = new T[1];
          x[0] = t;
          length = 1;    }

     bvector(int n, const T& t1, const T& t2) // n should be 2
     {    x = new T[2];
          x[0] = t1;
          x[1] = t2;
          length = 2;    }

     bvector(const bvector& v)
     {    if ( ! v.x )
          {    x = 0;
               length = 0;   }
          else
          {    length = v.length;
               x = new T[length];
               for ( longlong i = 0; i < length; i++ )
                    x[i] = v(i);    }    }

     ~bvector( ) { if ( x ) delete [ ] x; }

     void Setsize(longlong n) // resize, destroying contents
     {    if ( x && length == n ) return;
          if ( x ) delete [ ] x;
          x = new T[n];
          length = n;    }
     void resize(longlong n)
     {    T* x_new = new T[n];
          for ( longlong i = 0; i < min( n, length ); i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          length = n;    }
     void Append( const T& t )
     {    T* x_new = new T[ length + 1 ];
          for ( longlong i = 0; i < length; i++ )
               x_new[i] = x[i];
          delete [ ] x;
          x = x_new;    
          ++length;
          (*this)(length-1) = t;    }
     void Prepend( const T& t )
     {    ++length;
          T* x_new = new T[length];
          for ( longlong i = 1; i < length; i++ )
               x_new[i] = x[i-1];
          delete [ ] x;
          x = x_new;    
          (*this)(0) = t;    }
     bvector& operator=(const bvector& v)
     {    if ( !v.x )
          {    if ( x ) delete [ ] x;
               x = 0;    }
          else
          {    if ( !x )
               {    length = v.length;
                    x = new T[length];    }
               else if ( length != v.length )
               {    delete [ ] x;
                    length = v.length;
                    x = new T[length];    }
               for ( longlong i = 0; i < v.length; i++ )
                    (*this).x[i] = v(i);    }
           return *this;    }
     T& operator( )(longlong i) const 
     {    AssertLt( i, length );
          return x[i];    }

};

#endif
