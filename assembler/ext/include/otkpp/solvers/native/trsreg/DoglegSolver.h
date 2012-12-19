
#ifndef DOGLEGSOLVER_H

#include <otkpp/solvers/native/trsreg/TrustRegionSolver.h>

/// Implements the dogleg trust region algorithm.
class DoglegSolver : public TrustRegionSolver
{
 public:
   DoglegSolver();
   
  boost::numeric::ublas::vector< double > &computeStep(const boost::numeric::ublas::vector< double > &x,
                                double fx,
                                const boost::numeric::ublas::vector< double > &g,
                                const boost::numeric::ublas::matrix< double > &H,
                                boost::numeric::ublas::vector< double > &p,
                                bool &nonzeroStep,
                                boost::numeric::ublas::vector< double > &xPlus,
                                double &fxPlus,
                                const boost::numeric::ublas::vector< double > *Hg = NULL);
 private:
  double computeTau_(const boost::numeric::ublas::vector< double > &pB,
                     const boost::numeric::ublas::vector< double > &pU);
  
  void doSetup_();
};

#define DOGLEGSOLVER_H

#endif
