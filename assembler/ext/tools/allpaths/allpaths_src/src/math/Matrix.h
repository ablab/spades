///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines class "matrix" and some basic linear algebra functions.

#ifndef MATRIX_H
#define MATRIX_H

#include "CoreTools.h"
#include "math/Functions.h"
#include "math/Permutation.h"

template<class T> class matrix {

     public:

     int Nrows( ) const { return nrows_; }
     int Ncols( ) const { return ncols_; }
     Bool empty() const { return ncols_==0 || nrows_==0;}
     Bool Square( ) const { return Nrows( ) == Ncols( ); }

     vec<T>& operator[ ]( int i ) { return entries_[i]; }
     const vec<T>& operator[ ]( int i ) const { return entries_[i]; }
     T& operator( )( int i, int j ) { return entries_[i][j]; }
     const T& operator( )( int i, int j ) const { return entries_[i][j]; }

     matrix( ) : nrows_(0), ncols_(0) { }

     matrix( int r, int c );
     matrix( int r, int c, const T& t );

     void Resize( int r, int c );

     double Cofactor(int i, int j);

     void Invert(void);

     void Transpose( );

     void Identity(void) {
		if (Square()) {
			for (int rc = 0; rc < Nrows(); rc++) {
				entries_[rc][rc] = 1.0;
			}
		}
     }

     // lu_decompose: Compute the LU decomposition of an n by n invertible matrix 
     // A using row-oriented partial pivoting.  We find a unit lower triangular 
     // matrix L, an upper triangular matrix U, and a permutation matrix P such that 
     // A = PLU.  Upon entry *this is to be A.  Upon exit it is replaced by 
     // L - I + U.  If in the process of computation it is found that A is not 
     // invertible, False is returned; otherwise True is returned. 

     Bool lu_decompose(Permutation& P);

     // solve_from_lu: Start with an invertible matrix A and a vector
     // b.  Solve Ax = b for x.  Upon entry *this and P are to be the
     // LU-decomposition of A coming from lu_decompose.  b.size() and
     // x.size() must equal Nrows().  The code is based on the code
     // given in Forsythe and Moler, pp. 59--60.

     void solve_from_lu( const Permutation& P, const vec<T>& b, vec<T>& x )
     { x.resize(Nrows()); solve_from_lu(P, &b[0], &x[0]); }
    
     // For gcc 4.2 compliance, moved code body to below from MatrixTemplate.h.
     // This version should be usable no matter how the data are
     // stored as long as the elements are consecutive in memory.  The
     // vectors b and x must, of course, contain Nrows() elements each.

     void solve_from_lu( const Permutation& P, const T* b, T* x )
     {    matrix<T>& lu = *this;
          int i, j, n = Nrows( );
	  T dot;
	  for ( i = 0; i < n; i++ )
	  {    dot = 0;
	       T* v = &lu[i][0];
	       for ( j = 0; j < i; j++ )
		  dot += v[j] * x[j];
	       x[i] = b[P[i]-1] - dot;    }
	  for ( i = n-1; i >= 0; i-- )
	  {    dot = 0;
	       T* v = &lu[i][0];
	       for ( j = i+1; j < n; j++ )
	          dot += v[j] * x[j];
	       x[i] = (x[i] - dot) / v[i];    }    }

     void Print( ostream& o ) const;

     void PrettyPrint( ostream& o ) const;

     // Support for a bonehead implementation of inverse.

     T Det();
     void Inverse(matrix<T>& inverse);
     
     private:

     vec< vec<T> > entries_;
     int nrows_, ncols_;

};

// mul( A, x, Ax ): set the last argument to the product of the first two arguments.

template<class T> void mul( const matrix<T>& A, const vec<T>& x, vec<T>& Ax );

// mul( A, B, C ): set the last argument to the product of the first two arguments.

template<class T> void mul( const matrix<T>& A, const matrix<T>& B, matrix<T>& C );

// mulsub( A, x, b, Ax_minus_b ): compute Ax - b, returned as last argument.

template<class T> void mulsub( const matrix<T>& A, const vec<T>& x, 
     const vec<T>& b, vec<T>& Ax_minus_b );

// Solve: Find the solution x to the matrix equation Ax = b.  While
// useful for solving a single equation, this function is inefficient
// if you use the same A many times because it decomposes A every
// time.  Use the matrix class instead and compute the lu decomposition once. 
template<class T> void Solve( const matrix<T>& A, const vec<T>& b, vec<T>& x );

// SolveLeastSquares: For the matrix equation Ax = b, find the x that comes
// closest to solving it.

template<class T> void SolveLeastSquares( 
     const matrix<T>& A, const vec<T>& b, vec<T>& x );

// BestQuadratic: given vectors x and y of the same size, find the quadratic
// function f(t) = at^2 + bt + c which comes closest to satisfying the equation
// yi = f(xi) for each i.

template<class T> void BestQuadratic( const vec<T>& x, const vec<T>& y,
     T& a, T& b, T& c );

#endif
