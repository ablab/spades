
#ifndef LINEMINIMIZER_H

#include <otkpp/objfunc/Function.h>
#include <otkpp/lib/Cloneable.h>

#include <boost/numeric/ublas/vector.hpp>

class LineMinimizer;

/// A base class for line minimization algorithms.
class LineMinimizer
{
 public:
  /// Defines the parameters for a line minimizer.
  struct Setup
  {
    virtual ~Setup() { }
  
    virtual Setup *clone() const = 0;
    virtual bool isCompatibleWith(const LineMinimizer &s) const = 0;
  };

  /// Defines the default parameters for a line minimizer.
  struct DefaultSetup : public Cloneable< DefaultSetup, Setup >
  {
    bool isCompatibleWith(const LineMinimizer &s) const { return true; }
  };
  
  virtual ~LineMinimizer() { }
  
  /// Computes initial line search step length.
  /**
   * This method computes the initial line search step length 
   * with Fletcher's quadratic interpolation formula.
   */
  static double fletcherInitStep(double dphi0, double deltaF);
  
  /// Minimizes the objective function along the given search line.
  /**
   * This method minimizes the objective function along the given 
   * search line.
   * @param[in] x the starting point of the search line
   * @param[in] d search direction
   * @param[in] alpha0 initial step length
   * @param[in] fx objective function value at x
   * @param[in] gx gradient value at x
   * @param[out] alpha the computed step length \f$\alpha\f$
   * @param[out] x_plus the minimum point \f$\mathbf{x}_+=\mathbf{x}+\alpha\mathbf{d}\f$
   * @param[out] f_plus the minimum function value \f$f_+=f(\mathbf{x}_+)\f$
   * @param[out] g_plus the minimum gradient \f$\mathbf{g}_+=\nabla f(\mathbf{x}_+)\f$
   */
  virtual int minimize(const boost::numeric::ublas::vector< double > &x,
                       const boost::numeric::ublas::vector< double > &d,
                       double alpha0,
                       double fx,
                       const boost::numeric::ublas::vector< double > &gx,
                       double &alpha,
                       boost::numeric::ublas::vector< double > &x_plus,
                       double &f_plus,
                       boost::numeric::ublas::vector< double > &g_plus) = 0;
  /// Initializes this line minimizer.
  /**
   * Initializes this line minimizer with the given objective 
   * function and parameters.
   * @param f the objective function
   * @param s setup parameters
   */
  void setup(const Function &f, const Setup &s);
 protected:
  const Function *f_;
  
  virtual void doSetup_(const Setup &s) = 0;
};

#define LINEMINIMIZER_H

#endif
