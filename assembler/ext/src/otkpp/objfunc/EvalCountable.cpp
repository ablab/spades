
#include "EvalCountable.h"

EvalCountable::EvalCountable() : evalCounting_(false) { }

void EvalCountable::disableEvalCounting()
{
  evalCounting_ = false;
}

void EvalCountable::enableEvalCounting()
{
  evalCounting_ = true;
}

int EvalCountable::getEvalCounter() const
{
  return evalCounter_;
}

void EvalCountable::resetEvalCounter()
{
  evalCounter_ = 0;
}
