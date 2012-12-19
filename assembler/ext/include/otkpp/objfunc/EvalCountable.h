
#ifndef EVALCOUNTABLE_H

/// The base class of objects with evaluation counter.
class EvalCountable
{
 public:
  EvalCountable();
  
  /// Disables evaluation counting.
  void disableEvalCounting();
  
  /// Enables evaluation counting.
  void enableEvalCounting();
  
  /// Returns the evaluation counter.
  int getEvalCounter() const;
  
  /// Resets the evaluation counter.
  void resetEvalCounter();
 protected:
  mutable unsigned int evalCounter_;
  bool evalCounting_;
};

#define EVALCOUNTABLE_H

#endif
