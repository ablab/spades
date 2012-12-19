
#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/numeric/ublas/triangular.hpp>

namespace ublas = boost::numeric::ublas;


/// Decompose the symmetric positive definite matrix A into product L L^T.
/**
 * @param MATRIX type of input matrix 
 * @param TRIA type of lower triangular output matrix
 * @param A square symmetric positive definite input matrix (only the lower triangle is accessed)
 * @param L lower triangular output matrix 
 * @return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
 */
template < class MATRIX, class TRIA >
  size_t cholesky_decompose(const MATRIX& A, TRIA& L)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  assert( A.size1() == A.size2() );
  assert( A.size1() == L.size1() );
  assert( A.size2() == L.size2() );

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
        
    double qL_kk = A(k,k) - inner_prod( project( row(L, k), range(0, k) ),
                                        project( row(L, k), range(0, k) ) );
    
    if (qL_kk <= 0)
      return 1 + k;
    else
    {
      double L_kk = sqrt( qL_kk );
      L(k,k) = L_kk;
      
      matrix_column<TRIA> cLk(L, k);
      project( cLk, range(k+1, n) )
        = ( project( column(A, k), range(k+1, n) )
            - prod( project(L, range(k+1, n), range(0, k)), 
                    project(row(L, k), range(0, k) ) ) ) / L_kk;
    }
  }
  return 0;
}


/// Decompose the symmetric positive definite matrix A into product L L^T.
/**
 * @param MATRIX type of matrix A
 * @param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
 * @param A output: the lower triangle of A is replaced by the cholesky factor
 * @return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
 */
template < class MATRIX >
  size_t cholesky_decompose(MATRIX& A)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  const MATRIX& A_c(A);

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
        
    double qL_kk = A_c(k,k) - inner_prod( project( row(A_c, k), range(0, k) ),
                                          project( row(A_c, k), range(0, k) ) );
    
    if (qL_kk <= 0) {
      return 1 + k;
    } else {
      double L_kk = sqrt( qL_kk );
      
      matrix_column<MATRIX> cLk(A, k);
      project( cLk, range(k+1, n) )
        = ( project( column(A_c, k), range(k+1, n) )
            - prod( project(A_c, range(k+1, n), range(0, k)), 
                    project(row(A_c, k), range(0, k) ) ) ) / L_kk;
      A(k,k) = L_kk;
    }
  }
  return 0;      
}

/// Decompose the symmetric positive definit matrix A into product L L^T.
/**
 * @param MATRIX type of matrix A
 * @param A input: square symmetric positive definite matrix (only the lower triangle is accessed)
 * @param A output: the lower triangle of A is replaced by the cholesky factor
 * @return nonzero if decompositon fails (the value ist 1 + the numer of the failing row)
 */
template < class MATRIX >
  size_t incomplete_cholesky_decompose(MATRIX& A)
{
  using namespace ublas;

  typedef typename MATRIX::value_type T;
  
  // read access to a const matrix is faster
  const MATRIX& A_c(A);

  const size_t n = A.size1();
  
  for (size_t k=0 ; k < n; k++) {
    
    double qL_kk = A_c(k,k) - inner_prod( project( row( A_c, k ), range(0, k) ),
                                          project( row( A_c, k ), range(0, k) ) );
    
    if (qL_kk <= 0) {
      return 1 + k;
    } else {
      double L_kk = sqrt( qL_kk );

      // aktualisieren
      for (size_t i = k+1; i < A.size1(); ++i) {
        T* Aik = A.find_element(i, k);

        if (Aik != 0) {
          *Aik = ( *Aik - inner_prod( project( row( A_c, k ), range(0, k) ),
                                      project( row( A_c, i ), range(0, k) ) ) ) / L_kk;
        }
      }
        
      A(k,k) = L_kk;
    }
  }
        
  return 0;
}

/// Solve system L L^T x = b inplace.
/**
 * @param L a triangular matrix
 * @param x input: right hand side b; output: solution x
 */
template < class TRIA, class VEC > void cholesky_solve(const TRIA& L, VEC& x/*, ublas::lower*/)
{
  using namespace ublas;
//   ::inplace_solve(L, x, lower_tag(), typename TRIA::orientation_category () );
  inplace_solve(L, x, lower_tag() );
  inplace_solve(trans(L), x, upper_tag());
}

#endif
