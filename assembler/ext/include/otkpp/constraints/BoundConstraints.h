
#ifndef BOUNDCONSTRAINTS_H

#include <otkpp/constraints/Constraints.h>

#include <vector>

/// Defines bound constraints.
/**
 * This class defines constraints of the form 
 * \f$l_i\leq x_i\leq u_i\f$, \f$i=1,\dots,n\f$.
 */
struct BoundConstraints : public Cloneable< BoundConstraints, Constraints >
{
  /// Defines the type of a bound constraint.
  enum BoundType
  {
    NONE,  ///< no constraints
    LOWER, ///< lower bound
    UPPER, ///< upper bound
    BOTH   ///< both bounds
  };
  
  /// constraint types
  std::vector< BoundType > types;
  
  /// vector of lower bounds
  std::vector< double > L;
  
  /// vector of upper bounds
  std::vector< double > U;
  
  /// Default constructor, constructs empty constraints.
  BoundConstraints() { }
  
  /// Constructs bound constraints for n-dimensional problem.
  BoundConstraints(int n);
  
  bool isFeasible(const boost::numeric::ublas::vector< double > &x) const;
};

#define BOUNDCONSTRAINTS_H

#endif
