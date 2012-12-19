
#ifndef QUASINEWTONUPDATER_H

#include <boost/numeric/ublas/matrix.hpp>

/// Defines a matrix updater for quasi-Newton -based algorithms.
class QuasiNewtonUpdater
{
 public:
  virtual ~QuasiNewtonUpdater() { }
  
  /// Updates the Hessian or inverse Hessian approximation.
  /**
   * The Hessian (or inverse Hessian) approximation H is updated 
   * by using the difference vectors \f$\mathbf{p}_{k+1}=\mathbf{x}_{k+1}-
   * \mathbf{x}_{k}\f$ and \f$\mathbf{q}_{k+1}=\mathbf{g}_{k+1}-\mathbf{g}_{k}\f$.
   */
  virtual boost::numeric::ublas::matrix< double > &update(const boost::numeric::ublas::vector< double > &p,
                                   const boost::numeric::ublas::vector< double > &q,
                                   boost::numeric::ublas::matrix< double > &H) = 0;
};

#define QUASINEWTONUPDATER_H

#endif
