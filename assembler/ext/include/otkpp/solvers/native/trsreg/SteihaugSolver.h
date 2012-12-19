
#ifndef STEIHAUGSOLVER_H

#include <otkpp/solvers/native/trsreg/TrustRegionSolver.h>

/// Implements Steihaug's trust region algorithm.
class SteihaugSolver : public TrustRegionSolver
{
 public:
  SteihaugSolver();
  
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
  double computeTau_(const boost::numeric::ublas::vector< double > &d,
                     const boost::numeric::ublas::vector< double > &p);
  
  void doSetup_();
};

#define STEIHAUGSOLVER_H

#endif

