///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This file defines class "matrix" and some basic linear algebra functions.

#include "CoreTools.h"
#include "math/Functions.h"
#include "math/Permutation.h"

template<class T> matrix<T>::matrix( int r, int c )
{    nrows_ = r;
     ncols_ = c;
     entries_.resize(r);
     for ( int i = 0; i < r; i++ )
          entries_[i].resize(c);    }

template<class T> matrix<T>::matrix( int r, int c, const T& t )
{    nrows_ = r;
     ncols_ = c;
     entries_.resize(r);
     for ( int i = 0; i < r; i++ )
     {    entries_[i].resize(c);
          for ( int j = 0; j < c; j++ )
               entries_[i][j] = t;    }    }

template<class T> void matrix<T>::Resize( int r, int c )
{    if ( r > nrows_ )
     {    entries_.resize(r);
          for ( int i = nrows_; i < r; i++ )
               entries_[i].resize(c);    }
     if ( c != ncols_ )
     {    for ( int i = 0; i < nrows_; i++ )
               entries_[i].resize(c);    }
     nrows_ = r;
     ncols_ = c;    }

template<class T> double matrix<T>::Cofactor(int i, int j) {
	matrix<T> subm(nrows_ - 1, nrows_ - 1);

	int rpos = 0;
	for (int subi = 0; subi < nrows_; subi++) {
		if (i == subi) { continue; }
		int cpos = 0;
		for (int subj = 0; subj < ncols_; subj++) {
			if (j == subj) { continue; }
			subm(rpos, cpos) = entries_[subi][subj];
			cpos++;
		}
		rpos++;
	}	

	return powf(-1, i+j)*subm.Det();
}

template<class T> void matrix<T>::Invert(void) {
	Assert(Square());

	double det = this->Det();

	matrix<double> cofactor(nrows_, nrows_);

	for (int i = 0; i < nrows_; i++) {
		for (int j = 0; j < nrows_; j++) {
			cofactor(i, j) = (1.0/det)*this->Cofactor(i, j);
		}
	}

	for (int i = 0; i < nrows_; i++) {
		for (int j = 0; j < nrows_; j++) {
			entries_[j][i] = static_cast<T>(cofactor(i, j));
		}
	}
}

template<class T> void matrix<T>::Transpose( )
{    vec< vec<T> > new_entries(ncols_);
     for ( int i = 0; i < ncols_; i++ )
     {    new_entries[i].resize(nrows_);
          for ( int j = 0; j < nrows_; j++ )
               new_entries[i][j] = entries_[j][i];    }
     entries_ = new_entries;
     swap( nrows_, ncols_ );    }

template<class T> Bool matrix<T>::lu_decompose(Permutation& P)
{    
     // We start with P = I.  As we proceed through the function, P changes, and
     // for expository purposes, we let M denote the matrix product
     // P^{-1} x (*this).

     matrix<T>& me = *this;
     if ( !Square( ) ) FatalErr( "lu_decompose called with non-square matrix" );
     int n = Nrows( );
     static vector<int> pivotindices;
     pivotindices.resize(n-1);
     P = Permutation(n);
     int i, j, k, m, pivotindex = -1;
     T x, biggest, pivot, mult;

     // The current column is k.

     for ( k = 0; k < n-1; k++ )
     {    
          // Find the entry of M in the current column (on or below the diagonal) 
          // which has the largest absolute value.  Let "pivotindex" be its row 
          // index.

          biggest = 0;
          for ( i = k; i < n; i++ )
          {    x = Abs( me( P[i]-1, k ) );
               if ( x > biggest ) 
               {    biggest = x;
                    pivotindex = i;    }    }
          if ( biggest == T(0) ) return false;

          // Replace P by P \circ (k <--> pivotindex), thereby in effect 
          // swapping rows k and pivotindex of M.

          if ( pivotindex != k ) swap( P[k], P[pivotindex] );
          pivotindices[k] = pivotindex;

          // In effect execute the following code:
          // for ( i = k + 1; i < n; i++ )
          // {    M(i,k) /= M(k,k);
          //      for ( j = k + 1; j < n; j++ )
          //           M(i,j) -= M(i,k) * M(k,j);    }

          pivot = me( P[k]-1, k );
          m = P[k]-1;
          T* v = &me[m][0];
          for ( i = k+1; i < n; i++ )
          {    T* w = &me[ P[i]-1 ][0];
               mult = w[k] / pivot;
               w[k] = mult;
               if ( mult != T(0) )
               {    for ( j = k + 1; j < n; j++ )
                         w[j] -= mult * v[j];    }    }    }
     
     // At this point, *this is P(L - I + U).  We change it into
     // L - I + U.

     for ( i = 0; i <= n-2; i++ )
          if ( i != pivotindices[i] ) swap( me[i], me[pivotindices[i]] );
     
     return me( n-1, n-1 ) != T(0);    }

template<class T> T matrix<T>::Det()
{
    matrix<T> LU = *this;
    Permutation P_vec;
    LU.lu_decompose(P_vec);

    T det = 1; // I wish I could zero this, but I can't because it's a template. :/
    for (int i = 0; i < this->Nrows(); i++)
    {
        det *= LU(i,i);
    }
    return det;
}

template<class T> void matrix<T>::Inverse(matrix<T>& inverse)
{
    inverse = *this; // Just to get the same dimensions.

    // Get L and U (and P).
    matrix<T> LU = *this; 
    Permutation P_vec;
    LU.lu_decompose(P_vec);
    matrix<T> L = LU;
    matrix<T> U = LU;
    for (int i = 0; i < LU.Nrows(); i++)
    {
        for (int j = 0; j < LU.Ncols(); j++)
        {
            if      (j >= i) { L(i,j) = 0; } // Upper.
            else if (j <  i) { U(i,j) = 0; } // Lower.
            if      (i == j) { L(i,j) = 1; } // Got to restore the 1's in the lower.
        }
    }  
    matrix<T> P(LU.Nrows(), LU.Ncols(), 0);
    for (int i = 0; i < P.Nrows(); i++) 
    { 
        for (int j = 0; j < P.Ncols(); j++)
        {
            if (P_vec[i] == j+1) { P(i,j) = 1; }
        }
    }

#if 0
    L.PrettyPrint(cout); cout << endl;
    U.PrettyPrint(cout); cout << endl;
    P.PrettyPrint(cout); cout << endl;
#endif

    // Now solve for each column of A^-1.
    for (int i = 0; i < this->Ncols(); i++)
    {
        vec<T> bk(this->Nrows());
        for (int j = 0; j < this->Nrows(); j++)
        {
            if (i == j) { bk[j] = 1; }
            else        { bk[j] = 0; }
        }

        vec<T> Pbk;
        mul(P, bk, Pbk);

        vec<T> y;
        Solve(L, Pbk, y);

        vec<T> xk;
        Solve(U, y, xk);

        for(int j = 0; j < this->Nrows(); j++)
        {
            inverse(i,j) = xk[j];
        }
    }

	inverse.Transpose();
}

template<class T> void Solve( const matrix<T>& A, const vec<T>& b, vec<T>& x )
{    Permutation P;
     matrix<T> M = A;
     M.lu_decompose(P);
     M.solve_from_lu( P, b, x );    }

template<class T> void mul( const matrix<T>& A, const vec<T>& x, vec<T>& Ax )
{    int r = A.Nrows( ), c = A.Ncols( );
     ForceAssert( c == (int) x.size( ) );
     Ax.resize_and_set( r, 0 );
     for ( int i = 0; i < r; i++ )
          for ( int j = 0; j < c; j++ )
               Ax[i] += A[i][j] * x[j];    }

template<class T> void mul( const matrix<T>& A, const matrix<T>& B, matrix<T>& C )
{    ForceAssertEq( A.Ncols( ), B.Nrows( ) );
     C.Resize( A.Nrows( ), B.Ncols( ) );
     for ( int i = 0; i < A.Nrows( ); i++ )
     {    for ( int j = 0; j < B.Ncols( ); j++ )
          {    T x = 0;
               for ( int k = 0; k < A.Ncols( ); k++ )
                    x += A(i,k) * B(k,j);
               C(i,j) = x;    }    }    }

template<class T> void mulsub( const matrix<T>& A, const vec<T>& x, 
     const vec<T>& b, vec<T>& Ax_minus_b )
{    mul( A, x, Ax_minus_b );
     ForceAssert( b.size( ) == Ax_minus_b.size( ) );
     for ( int i = 0; i < (int) b.size( ); i++ )
          Ax_minus_b[i] -= b[i];    }

template<class T> void SolveLeastSquares( 
     const matrix<T>& A, const vec<T>& b, vec<T>& x )
{    matrix<T> At = A, AtA;
     At.Transpose( );
     mul( At, A, AtA );
     vec<T> Atb;
     mul( At, b, Atb );
     Solve( AtA, Atb, x );    }

template<class T> void BestQuadratic( const vec<T>& x, const vec<T>& y,
     T& a, T& b, T& c )
{    ForceAssertEq( x.size( ), y.size( ) );
     ForceAssertGt( x.size( ), 0u );
     matrix<T> M( x.size( ), 3 );
     for ( size_t i = 0; i < x.size(); i++ )
     {    M(i,0) = x[i] * x[i],
          M(i,1) = x[i];
          M(i,2) = 1;    }
     vec<T> abc;
     SolveLeastSquares( M, y, abc );
     a = abc[0], b = abc[1], c = abc[2];    }

template<class T> void matrix<T>::PrettyPrint( ostream& o ) const
{    vec< vec<String> > rows;
     for ( int i = 0; i < Nrows( ); i++ )
     {    vec<String> row;
          for ( int j = 0; j < Ncols( ); j++ )
               row.push_back( ToString( (double)(*this)(i,j), 10 ) ); //JRM
          rows.push_back(row);    }
     PrintTabular( o, rows, 2 );    }

#define MATRIX_DEF(T)                                                           \
     template matrix<T>::matrix( int r, int c );                                \
     template matrix<T>::matrix( int r, int c, const T& t );                    \
     template void matrix<T>::Resize( int r, int c );                           \
     template double matrix<T>::Cofactor(int i, int j);                         \
     template void matrix<T>::Invert(void);                                     \
     template void matrix<T>::Transpose( );                                     \
     template T matrix<T>::Det( );                                              \
     template void matrix<T>::Inverse(matrix<T>& inverse);                                     \
     template Bool matrix<T>::lu_decompose(Permutation& P);                     \
     template void matrix<T>::solve_from_lu(                                    \
          const Permutation& P, const vec<T>& b, vec<T>& x );                   \
     template void matrix<T>::PrettyPrint( ostream& ) const;                    \
     template void Solve( const matrix<T>& A, const vec<T>& b, vec<T>& x );     \
     template void mul( const matrix<T>& A, const vec<T>& x, vec<T>& Ax );      \
     template void mul( const matrix<T>& A, const matrix<T>& B, matrix<T>& C ); \
     template void mulsub( const matrix<T>& A, const vec<T>& x,                 \
          const vec<T>& b, vec<T>& Ax_minus_b );                                \
     template void SolveLeastSquares(                                           \
          const matrix<T>& A, const vec<T>& b, vec<T>& x );                     \
     template void BestQuadratic( const vec<T>& x, const vec<T>& y,             \
          T& a, T& b, T& c );
