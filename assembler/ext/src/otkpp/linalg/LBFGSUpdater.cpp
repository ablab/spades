
#include "cholesky.hpp"
#include "LBFGSUpdater.h"

#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

LBFGSUpdater::LBFGSUpdater(int n, int m) : 
  n_(n), m_(m),
  I_(identity_matrix< double >(n)),
  invD_(zero_matrix< double >(m, m)),
  L_(zero_matrix< double >(m)),
  S_(zero_matrix< double >(n, m)),
  sqrtD_(zero_matrix< double >(m, m)),
  sqrtInvD_(zero_matrix< double >(m, m)),
  SS_(matrix< double >(m, m)),
  Y_(zero_matrix< double >(n, m)),
  JJ_(matrix< double >(m, m)),
  JJ_chol_(matrix< double >(m, m)),
  M_(matrix< double >(2*m, 2*m)),
  P_(permutation_matrix< double >(m)),
  p_(vector< double >(2*m)),
  q_(vector< double >(2*m)),
  tri_(matrix< double >(2*m, 2*m)) { }

bool LBFGSUpdater::computeProduct(const vector< double > &s,
                                  const vector< double > &y,
                                  const vector< double > &v,
                                  vector< double > &Bv)
{
  double sigma = inner_prod(y, s) / inner_prod(s, s);
  JJ_ = sigma * prod(trans(S_), S_) + 
        prod(matrix< double >(prod(L_, invD_)), trans(L_));
  
  project(p_, range(0, m_))    = prod(trans(Y_), v);
  project(p_, range(m_, 2*m_)) = sigma * prod(trans(S_), v);
  
  int cr = cholesky_decompose(JJ_, JJ_chol_);
  
  if(cr != 0)
    return false;
  
  matrix_range< matrix< double > > M1(M_, range(0, m_), range(0, m_));
  matrix_range< matrix< double > > M2(M_, range(0, m_), range(m_, 2*m_));
  matrix_range< matrix< double > > M3(M_, range(m_, 2*m_), range(0, m_));
  matrix_range< matrix< double > > M4(M_, range(m_, 2*m_), range(m_, 2*m_));
  
  M1 = sqrtD_;
  M2 = zero_matrix< double >(m_, m_);
  M3 = -prod(L_, sqrtInvD_);
  M4 = JJ_chol_;
  
  try {
    inplace_solve(M_, p_, lower_tag());
  }
  catch(singular) {
    return false;
  }
  
  M1 = -sqrtD_;
  M2 = prod(sqrtInvD_, trans(L_));
  M3 = zero_matrix< double >(m_, m_);
  M4 = trans(JJ_chol_);
  
  try {
    inplace_solve(M_, p_, upper_tag());
  }
  catch(singular) {
    return false;
  }
  
  q_  = prod(Y_, project(p_, range(0, m_)));
  q_ += prod(sigma*S_, project(p_, range(m_, 2*m_)));
  
  Bv = sigma*v - q_;
  
  return true;
}

void LBFGSUpdater::updateVectors(const vector< double >&s,
                                 const vector< double >&y)
{
  updateInvD_(s, y);
  updateL_(s, y);
  updateSY_(s, y);
  updateSS_();
}

void LBFGSUpdater::updateInvD_(const vector< double >&s,
                               const vector< double >&y)
{
  for(int i = 0; i < m_ - 1; i++)
  {
    invD_(i, i)     = invD_(i+1, i+1);
    sqrtD_(i, i)    = sqrtD_(i+1, i+1);
    sqrtInvD_(i, i) = sqrtInvD_(i+1, i+1);
  }
  double sy = inner_prod(s, y);
  invD_(m_-1, m_-1)     = 1.0 / sy;
  sqrtD_(m_-1, m_-1)    = sqrt(sy);
  sqrtInvD_(m_-1, m_-1) = 1.0 / sqrt(sy);
}

void LBFGSUpdater::updateL_(const vector< double > &s,
                            const vector< double > &y)
{
  int i, j;
  
  /*L_ = prod(trans(S_), Y_);
  for(int i = 0; i < m_; i++)
  {
    for(int j = 0; j < m_; j++)
    {
      if(i <= j)
        L_(i, j) = 0.0;
    }
  }*/
  
  for(i = 0; i < m_ - 1; i++)
    for(j = 0; j < i - 1; j++)
      L_(i, j) = L_(i+1, j+1);
  
  for(i = 0; i < m_ - 1; i++)
    L_(m_-1, i) = inner_prod(s, column(Y_, i));
  L_(m_-1, m_-1) = inner_prod(s, y);
}

void LBFGSUpdater::updateSS_()
{
  SS_ = prod(trans(S_), S_);
}

void LBFGSUpdater::updateSY_(const vector< double > &s,
                             const vector< double > &y)
{
  int i, j;
  
  for(i = 0; i < n_; i++)
  {
    for(j = 0; j < m_ - 1; j++)
    {
      S_(i, j) = S_(i, j+1);
      Y_(i, j) = Y_(i, j+1);
    }
    S_(i, m_-1) = s[i];
    Y_(i, m_-1) = y[i];
  }
}
