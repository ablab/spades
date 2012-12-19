
#ifndef TRUSTREGIONSOLVER_H

#include <otkpp/objfunc/Function.h>

/// A base class for trust region algorithms.
class TrustRegionSolver
{
 public:
  struct Setup { };
  
  /// Computes a step vector.
  /**
   * This method computes the step that minimizes the objective 
   * function inside the current trust region radius. The trust 
   * region radius may also be updated, if necessary.
   * @param[in] x the center point of the trust region.
   * @param[in] fx the function value \f$f(\mathbf{x})\f$.
   * @param[in] g the gradient \f$\nabla f(\mathbf{x})\f$.
   * @param[in] H the current Hessian approximation.
   * @param[out] p the computed step vector.
   * @param[out] nonzeroStep true, if the step computation was successful, false if not.
   * @param[out] xPlus the new iterate \f$\mathbf{x}_+=\mathbf{x}+\mathbf{p}\f$.
   * @param[out] fxPlus the new function value \f$f_+=f(\mathbf{x}+\mathbf{p})\f$.
   * @param[in] Hg optional precomputed product of H and g (gradient).
   */
  virtual boost::numeric::ublas::vector< double > &computeStep(const boost::numeric::ublas::vector< double > &x,
                                        double fx,
                                        const boost::numeric::ublas::vector< double > &g,
                                        const boost::numeric::ublas::matrix< double > &H,
                                        boost::numeric::ublas::vector< double > &p,
                                        bool &nonzeroStep,
                                        boost::numeric::ublas::vector< double > &xPlus,
                                        double &fxPlus,
                                        const boost::numeric::ublas::vector< double > *Hg = NULL) = 0;
  
  double getDelta() const;
  
  void setDelta(double delta);
  
  /// Initializes this trust region solver with the given objective function.
  void setup(const Function &f);
 protected:
  double delta_;
  double deltaMax_;
  double eps_;
  double eta_;
  const Function *f_;
  
  double computeReduction_(double fx,
                           double fxPlus,
                           const boost::numeric::ublas::vector< double > &x,
                           const boost::numeric::ublas::vector< double > &g,
                           const boost::numeric::ublas::matrix< double > &H,
                           const boost::numeric::ublas::vector< double > &p);
  
  virtual void doSetup_() = 0;
  
  bool updateRadius_(double fx,
                     double fxPlus,
                     const boost::numeric::ublas::vector< double > &x,
                     const boost::numeric::ublas::vector< double > &g,
                     const boost::numeric::ublas::matrix< double > &H,
                     const boost::numeric::ublas::vector< double > &p);
};

#define TRUSTREGIONSOLVER_H

#endif
