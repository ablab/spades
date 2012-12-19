
#ifndef INVLBFGSUPDATER_H

#include <boost/numeric/ublas/matrix.hpp>
#include <list>

/// Implements the two-loop recursion inverse L-BFGS formula.
/**
 * This class implements the two-loop recursion inverse L-BFGS 
 * formula. It never explicitly stores the inverse Hessian 
 * approximation. Instead it stores and updates its product 
 * with negative gradient.
 */
class InvLBFGSUpdater
{
 public:
  InvLBFGSUpdater() { }
  
  /// Constructs a new inverse L-BFGS updater.
  /**
   * Constructs a new L-BFGS updater for n-dimensional 
   * inverse Hessian approximations with iteration history 
   * of length m.
   */
  InvLBFGSUpdater(int m, int n);
  
  /// Clears iteration history.
  void clearHistory();
  
  /// Removes the oldest pair from the iteration history.
  void removeOldestPair();
  
  /// Stores a new correction pair.
  void storePair(const boost::numeric::ublas::vector< double > &p,
                 const boost::numeric::ublas::vector< double > &q);
  
  /// Takes one update step.
  /**
   * Updates the product of the inverse Hessian approximation 
   * and the negative gradient and returns the updated vector.
   * @param g gradient vector.
   * @param d the vector into which the updated 
   *          product of the inverse Hessian approximation 
   *          and the negative gradient is stored.
   * @return the computed d-vector
   */
  boost::numeric::ublas::vector< double > &update(const boost::numeric::ublas::vector< double > &g,
                           boost::numeric::ublas::vector< double > &d);
 private:
  boost::numeric::ublas::vector< double > alpha_;
  boost::numeric::ublas::matrix< double > I_;
  int m_, n_;
  std::list< boost::numeric::ublas::vector< double > > P_;
  std::list< boost::numeric::ublas::vector< double > > Q_;
  boost::numeric::ublas::vector< double > r_;
  boost::numeric::ublas::vector< double > rho_;
};

#define INVLBFGSUPDATER_H

#endif
