
#ifndef XDISTTOMINTEST_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/solvers/native/NativeSolver.h>
#include <otkpp/stopcrit/StoppingCriterion.h>

/// Implements stopping criterion for \f$\|\mathbf{x}_{k}-\mathbf{x}^*\|\f$.
class XDistToMinTest : public Cloneable< XDistToMinTest, StoppingCriterion >
{
 public:
  /// Constructs a new stopping criterion.
  /**
   * @param xMin the known minimizer
   * @param eps the threshold value
   * @param relative if set to true, the relative distance 
   *                 \f$\frac{\|\mathbf{x}_{k}-\mathbf{x}^*\|}{\|\mathbf{x}^*\|}\f$ 
   *                 is tested instead.
   */
  XDistToMinTest(const boost::numeric::ublas::vector< double > &xMin, double eps, bool relative);
  
  double getTestValue(const NativeSolver &s) const;
  const boost::numeric::ublas::vector< double > &getXMin() const;
  bool test(const NativeSolver &s) const;
 private:
  double eps_;
  bool relative_;
  const boost::numeric::ublas::vector< double > xMin_;
};

#define XDISTTOMINTEST_H

#endif
