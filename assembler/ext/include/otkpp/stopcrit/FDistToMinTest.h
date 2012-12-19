
#ifndef FDISTTOMINTEST_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/solvers/native/NativeSolver.h>
#include <otkpp/stopcrit/StoppingCriterion.h>

/// Implements stopping criterion for \f$f(\mathbf{x}_{k})-f(\mathbf{x}^*)\f$.
class FDistToMinTest : public Cloneable< FDistToMinTest, StoppingCriterion >
{
 public:
  /// Constructs a new stopping criterion.
  /**
   * @param fMin the known minimum function value
   * @param eps the threshold value
   * @param relative if set to true, the relative distance 
   *                 \f$\frac{f(\mathbf{x}_{k})-f(\mathbf{x}^*)}{f(\mathbf{x}^*)}\f$ 
   *                 is tested instead.
   */
  FDistToMinTest(double fMin, double eps, bool relative);
  
  double getEps() const;
  double getFMin() const;
  double getTestValue(const NativeSolver &s) const;
  bool isRelative() const;
  bool test(const NativeSolver &s) const;
 private:
  double eps_;
  double fMin_;
  bool relative_;
};

#define FDISTTOMINTEST_H

#endif
