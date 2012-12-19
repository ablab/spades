
#include "BFGSUpdater.h"

#include <boost/numeric/ublas/blas.hpp>

using namespace boost::numeric::ublas;

BFGSUpdater::BFGSUpdater(BFGSUpdater::Type type) : type_(type) { }

matrix< double > &BFGSUpdater::update(const vector< double > &p,
                                      const vector< double > &q,
                                      matrix< double > &H)
{
  if(type_ == DIRECT)
  {
    vector< double > Hp = prod(H, p);
    H += outer_prod(q, q) / inner_prod(q, p) - 
         outer_prod(Hp, Hp) / inner_prod(p, Hp);
  }
  else if(type_ == INVERSE)
  {
    double pq = inner_prod(p, q);
    Hq_ = prod(H, q);
    double qHq = inner_prod(q, Hq_);
    double m = 1.0 + qHq / pq;
    
    //blas_2::sr(H, m / pq, p);
    //blas_2::sr2(H, -1.0 / pq, p, Hq_);
    H += m / pq * outer_prod(p, p) - (outer_prod(p, Hq_) + outer_prod(Hq_, p)) / pq;
  }
  
  return H;
}
