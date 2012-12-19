
#ifndef FLETCHER_H

#include <otkpp/solvers/native/linmin/LineMinimizer.h>

/// Implements Fletcher's line search algorithm.
class Fletcher : public LineMinimizer
{
 public:
  struct Setup : public Cloneable< Setup, LineMinimizer::Setup >
  {
    double eta;
    double mu;
    double chi;
    double tau;
  
    Setup(double eta = 0.1, double mu = 0.01, double chi = 9.0, double tau = 0.05);
    bool isCompatibleWith(const LineMinimizer &s) const;
  };

  int minimize(const boost::numeric::ublas::vector< double > &x,
               const boost::numeric::ublas::vector< double > &d,
               double alpha0,
               double fx,
               const boost::numeric::ublas::vector< double > &gx,
               double &alpha,
               boost::numeric::ublas::vector< double > &x_plus,
               double &f_plus,
               boost::numeric::ublas::vector< double > &g_plus);
 private:
  double alphal_, alphat_, alphau_;
  double phil_, phit_;
  double dphil_, dphit_;
  boost::numeric::ublas::vector< double > xt_;
  boost::numeric::ublas::vector< double > gt_;
  Setup setup_;
  
  void doSetup_(const LineMinimizer::Setup &s);
  double extrap_safequard_(double alpha);
  double interp_safequard_(double alpha);
};

#define FLETCHER_H

#endif
