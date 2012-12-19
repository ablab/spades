
#include "InvLBFGSUpdater.h"

using namespace boost::numeric::ublas;

InvLBFGSUpdater::InvLBFGSUpdater(int m, int n) : m_(m), n_(n)
{
  alpha_.resize(m);
  I_ = identity_matrix< double >(n);
  r_.resize(n);
  rho_.resize(m);
}

void InvLBFGSUpdater::clearHistory()
{
  P_.clear();
  Q_.clear();
}

void InvLBFGSUpdater::removeOldestPair()
{
  if(P_.size() > m_)
  {
    P_.pop_front();
    Q_.pop_front();
  }
}

void InvLBFGSUpdater::storePair(const vector< double > &p,
                                const vector< double > &q)
{
  P_.push_back(p);
  Q_.push_back(q);
}

vector< double > &InvLBFGSUpdater::update(const vector< double > &g,
                                          vector< double > &d)
{
  int i;
  int m;
  double beta;
  double gamma;
  std::list< vector< double > >::const_iterator p_i;
  std::list< vector< double > >::const_iterator q_i;
  
  m = std::min((int)P_.size(), m_);
  
  r_ = g;
  p_i = P_.end();
  q_i = Q_.end();
  for(i = 0; i < m; i++)
  {
    p_i--;
    q_i--;
    
    if(i == 0)
      gamma = inner_prod(*p_i, *q_i) / inner_prod(*q_i, *q_i);
    
    rho_[i] = 1.0 / inner_prod(*p_i, *q_i);
    alpha_[i] = rho_[i] * inner_prod(*p_i, r_);
    r_ -= alpha_[i] * *q_i;
  }
  
  d = -gamma * prod(I_, r_);
  p_i = P_.begin();
  q_i = Q_.begin();
  for(i = m - 1; i >= 0; i--)
  {
    beta = rho_[i] * inner_prod(*q_i, d);
    d -= *p_i * (alpha_[i] + beta);
    
    p_i++;
    q_i++;
  }
  
  return d;
}
