
#ifndef LBFGSUPDATER_H

#include <otkpp/linalg/QuasiNewtonUpdater.h>

#include <boost/numeric/ublas/lu.hpp>

// Implements the L-BFGS formula for Hessian approximations.
class LBFGSUpdater
{
 public:
  /// Constructs a new L-BFGS updater.
  /**
   * Constructs a new L-BFGS updater for an n-dimensional 
   * Hessian approximations with iteration history 
   * of length m.
   */
  LBFGSUpdater(int n, int m);
  
  bool computeProduct(const boost::numeric::ublas::vector< double > &s,
                      const boost::numeric::ublas::vector< double > &y,
                      const boost::numeric::ublas::vector< double > &v,
                      boost::numeric::ublas::vector< double > &Bv);
  
  /// Removes the oldest iterates from iteration history.
  void removeOldest();
  
  /// Processes a new pair of correction vectors s, y.
  /** 
   * @param s represents \f$\mathbf{s}_{k+1}=\mathbf{x}_{k+1}-\mathbf{x}_{k}\f$.
   * @param y represents \f$\mathbf{y}_{k+1}=\mathbf{g}_{k+1}-\mathbf{g}_{k}\f$.
   */
  void updateVectors(const boost::numeric::ublas::vector< double >&s,
                     const boost::numeric::ublas::vector< double >&y);
 private:
  int n_, m_;
  const boost::numeric::ublas::matrix< double > I_;
  boost::numeric::ublas::matrix< double > invD_;
  boost::numeric::ublas::matrix< double > L_;
  boost::numeric::ublas::matrix< double > S_;
  boost::numeric::ublas::matrix< double > sqrtD_;
  boost::numeric::ublas::matrix< double > sqrtInvD_;
  boost::numeric::ublas::matrix< double > SS_;
  boost::numeric::ublas::matrix< double > Y_;
  
  boost::numeric::ublas::matrix< double > JJ_;
  boost::numeric::ublas::matrix< double > JJ_chol_;
  boost::numeric::ublas::matrix< double > M_;
  boost::numeric::ublas::permutation_matrix< double > P_;
  boost::numeric::ublas::vector< double > p_, q_;
  boost::numeric::ublas::matrix< double > tri_;
  
  void updateInvD_(const boost::numeric::ublas::vector< double >&s,
                   const boost::numeric::ublas::vector< double >&y);
  void updateL_(const boost::numeric::ublas::vector< double > &s,
                const boost::numeric::ublas::vector< double > &y);
  void updateSS_();
  void updateSY_(const boost::numeric::ublas::vector< double > &s,
                 const boost::numeric::ublas::vector< double > &y);
};

#define LBFGSUPDATER_H

#endif
