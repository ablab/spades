
#include "SR1Updater.h"

using namespace boost::numeric::ublas;

matrix< double > &SR1Updater::update(const vector< double > &p,
                                     const vector< double > &q,
                                     matrix< double > &H)
{
  vector< double > t = q - prod(H, p);
  double pt = inner_prod(p, t);
  
  if(fabs(pt) >= 1e-8 * norm_2(p)*norm_2(t))
    H += outer_prod(t, t) / pt;
  
  return H;
}
