
#ifndef MORETHUENTE_H

#include <otkpp/solvers/native/linmin/LineMinimizer.h>

/// Implements the More and Thuente line search algorithm.
class MoreThuente : public LineMinimizer
{
 public:
   struct Setup : public Cloneable< Setup, LineMinimizer::Setup >
  {
    double eta;
    double mu;
    double gamma;
    double chi;
    
    Setup(double eta = 0.1, double mu = 0.001, double gamma = 0.0, double chi = 1e-2);
    bool isCompatibleWith(const LineMinimizer &s) const;
  };

  MoreThuente();
  
  int minimize(const boost::numeric::ublas::vector< double > &x,
               const boost::numeric::ublas::vector< double > &d,
               double alpha0,
               double fx,
               const boost::numeric::ublas::vector< double > &gx,
               double &alpha,
               boost::numeric::ublas::vector< double > &x_plus,
               double &f_plus,
               boost::numeric::ublas::vector< double > &g_plus);
  void setGamma(double gamma);
 private:
  double alphal_, alphat_, alphau_;
  double psil_, psit_, psiu_;
  double dpsil_, dpsit_, dpsiu_;
  double phil_, phit_, phiu_;
  double dphil_, dphit_, dphiu_;
  boost::numeric::ublas::vector< double > xt_;
  boost::numeric::ublas::vector< double > gt_;
  boost::numeric::ublas::vector< double > gu_;
  int alphau_evaluated_;
  double gamma_;
  Setup setup_;
  
  double trialstep(double alphat, double alphal, double alphau,
                   double fl, double dfl,
                   double ft, double dft,
                   double fu, double dfu,
                   const Function &f,
                   const boost::numeric::ublas::vector< double > &x,
                   const boost::numeric::ublas::vector< double > &d);
  void bracket(double &alphal, double alphat, double &alphau,
               double fl, double ft, double dft,
               double &phil, double phit, double &phiu,
               double &dphil, double dphit, double &dphiu);
  double psi(double phi, double alpha, double phi0, double dphi0, double mu);
  void doSetup_(const LineMinimizer::Setup &s);
};

#define MORETHUENTE_H

#endif
