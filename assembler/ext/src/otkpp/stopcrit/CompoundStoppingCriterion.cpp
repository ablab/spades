
#include "CompoundStoppingCriterion.h"

#include <stdexcept>

CompoundStoppingCriterion::CompoundStoppingCriterion(const CompoundStoppingCriterion &sc)
{
  std::list< const StoppingCriterion * >::const_iterator it;
  for(it = sc.stopCrit_.begin(); it != sc.stopCrit_.end(); it++)
    stopCrit_.push_back((*it)->clone());
}

CompoundStoppingCriterion::~CompoundStoppingCriterion()
{
  std::list< const StoppingCriterion * >::iterator it;
  for(it = stopCrit_.begin(); it != stopCrit_.end(); it++)
    delete (*it);
}

CompoundStoppingCriterion &CompoundStoppingCriterion::operator=(const CompoundStoppingCriterion &sc)
{
  std::list< const StoppingCriterion * >::iterator it1;
  std::list< const StoppingCriterion * >::const_iterator it2;
  
  for(it1 = stopCrit_.begin(); it1 != stopCrit_.end(); it1++)
    delete (*it1);
  stopCrit_.clear();
  for(it2 = sc.stopCrit_.begin(); it2 != sc.stopCrit_.end(); it2++)
    stopCrit_.push_back((*it2)->clone());
  
  return *this;
}

CompoundStoppingCriterion CompoundStoppingCriterion::operator+(const StoppingCriterion &sc)
{
  CompoundStoppingCriterion result(*this);
  result.stopCrit_.push_back(sc.clone());
  
  return result;
}

double CompoundStoppingCriterion::getTestValue(const NativeSolver &s) const
{
  return 0.0; //throw std::runtime_error("operation not implemented");
}

bool CompoundStoppingCriterion::test(const NativeSolver &s) const
{
  std::list< const StoppingCriterion * >::const_iterator it;
  for(it = stopCrit_.begin(); it != stopCrit_.end(); it++)
  {
    if((*it)->test(s))
      return true;
  }
  
  return false;
}
