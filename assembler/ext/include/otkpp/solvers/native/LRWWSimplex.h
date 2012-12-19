
#ifndef LRWWSIMPLEX_H

#include <otkpp/solvers/native/NativeSolver.h>

class LRWWSimplex : public NativeSolver
{
 public:
  class Simplex
  {
   public:
    enum OpType
    {
      EXPANSION,
      INSIDE_CONTRACTION,
      OUTSIDE_CONTRACTION,
      REFLECTION
    };
    
    Simplex() { }
    Simplex(int n);
    void computeCentroid();
    void computeExpansion(boost::numeric::ublas::vector< double > &xe);
    void computeInsideContraction(boost::numeric::ublas::vector< double > &xic);
    void computeOutsideContraction(boost::numeric::ublas::vector< double > &xoc);
    void computeReflection(boost::numeric::ublas::vector< double > &xr);
    const boost::numeric::ublas::vector< double > &getCentroidPoint() const;
    double getfx(int i) const;
    const boost::numeric::ublas::matrix_column< const boost::numeric::ublas::matrix< double > > getx(int i) const;
    void improveHighestVertex(const boost::numeric::ublas::vector< double > &x, double fx, OpType opType);
    void setVertex(int i, const boost::numeric::ublas::vector< double > &x, double fx);
    void shrink(Function &f);
    void sortVertices();
   private:
    int n_;                // the number of vertices
    boost::numeric::ublas::matrix< double > X_;   // (n+1)*n matrix of simplex vertices, one vertex per each matrix column
    boost::numeric::ublas::vector< double > fx_;  // (n+1)-dimensional vector of function values at the vertices
    double vol_;           // simplex volume
    boost::numeric::ublas::vector< double > xc_;  // the previously computed centroid point
  };
  
  struct State : public Cloneable< State, NativeSolver::State >
  {
    Simplex S;
  };
  
  unsigned int getM() const;
  std::string getName() const;
  const State &getState() const { return state_; }
  const boost::numeric::ublas::matrix< double > getXArray() const;
  bool hasBuiltInStoppingCriterion() const { return false; }
  bool usesGradient() const;
  bool usesHessian() const;
 private:
  State state_;
  boost::numeric::ublas::vector< double > xe_;
  boost::numeric::ublas::vector< double > xic_;
  boost::numeric::ublas::vector< double > xoc_;
  boost::numeric::ublas::vector< double > xr_;
  
  NativeSolver::IterationStatus iterate_();
  void doSetup_(const Function &objFunc,
                const boost::numeric::ublas::vector< double > &x0,
                const Solver::Setup &solverSetup,
                const Constraints &C);
};

#define LRWWSIMPLEX_H

#endif
