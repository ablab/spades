
#ifndef BFGSUPDATER_H

#include <otkpp/linalg/QuasiNewtonUpdater.h>

#include <boost/numeric/ublas/matrix.hpp>

/// Implements the standard BFGS update formula.
class BFGSUpdater : public QuasiNewtonUpdater
{
 public:
  /// BFGS update type.
  enum Type
  {
    DIRECT,  /// update Hessian approximation
    INVERSE  /// update Inverse Hessian approximation
  };
  
  /// Constructs a new BFGS updater of the given type.
  BFGSUpdater(Type type = INVERSE);
  
  boost::numeric::ublas::matrix< double > &update(const boost::numeric::ublas::vector< double > &p,
                           const boost::numeric::ublas::vector< double > &q,
                           boost::numeric::ublas::matrix< double > &H);
 private:
  boost::numeric::ublas::vector< double > Hq_;
  Type type_;
};

#define BFGSUPDATER_H

#endif
