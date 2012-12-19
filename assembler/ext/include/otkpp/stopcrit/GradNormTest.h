
#ifndef GRADNORMTEST

#include <otkpp/lib/Cloneable.h>
#include <otkpp/stopcrit/StoppingCriterion.h>

/// Implements stopping criterion for gradient norm \f$\|\nabla f(\mathbf{x}_k)\|\f$.
class GradNormTest : public Cloneable< GradNormTest, StoppingCriterion >
{
 public:
  /// Constructs gradient norm test with the given threshold value.
  GradNormTest(double eps);
  
  double getTestValue(const NativeSolver &s) const;
  bool test(const NativeSolver &s) const;
 private:
  double eps_;
};

#define GRADNORMTEST

#endif
