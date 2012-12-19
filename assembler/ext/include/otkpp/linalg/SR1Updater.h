
#ifndef SR1UPDATER_H

#include <otkpp/linalg/QuasiNewtonUpdater.h>

#include <boost/numeric/ublas/matrix.hpp>

/// Implements the SR1 formula for updating Hessian approximations.
class SR1Updater : public QuasiNewtonUpdater
{
 public:
  boost::numeric::ublas::matrix< double > &update(const boost::numeric::ublas::vector< double > &p,
                           const boost::numeric::ublas::vector< double > &q,
                           boost::numeric::ublas::matrix< double > &H);
};

#define SR1UPDATER_H

#endif
