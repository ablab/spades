
#include "CompoundStoppingCriterion.h"
#include "StoppingCriterion.h"

CompoundStoppingCriterion StoppingCriterion::operator+(const StoppingCriterion &sc)
{
  CompoundStoppingCriterion result;
  result = result + *this;
  return result + sc;
}
