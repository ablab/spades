
#ifndef MAXNUMITERTEST_H

#include <otkpp/lib/Cloneable.h>
#include <otkpp/stopcrit/StoppingCriterion.h>

class MaxNumIterTest : public Cloneable< MaxNumIterTest, StoppingCriterion >
{
 public:
  MaxNumIterTest(int n);
  
  double getTestValue(const NativeSolver &s) const;
  bool test(const NativeSolver &s) const;
 private:
  double n_;
};

#define MAXNUMITERTEST_H

#endif
