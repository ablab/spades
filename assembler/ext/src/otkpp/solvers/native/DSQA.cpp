
#include "DSQA.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

std::string DSQA::getName() const
{
  return "DSQA";
}

bool DSQA::hasBuiltInStoppingCriterion() const
{
  return false;
}

bool DSQA::usesGradient() const
{
  return false;
}

bool DSQA::usesHessian() const
{
  return false;
}

/*void DSQA::allocateState_()
{
  NativeSolver::state_ = &state_;*/
  /*NativeSolver::state_ = new DSQA::State();
  state_ = dynamic_cast< DSQA::State * >(NativeSolver::state_.get_pointer());*/
//}

double DSQA::computeReduction_(const vector< double > &x,
                               const vector< double > &xPlus,
                               double f,
                               double fPlus,
                               const vector< double > &p) const
{
  double ared, pred;
  
  ared = f - fPlus;
  pred = -inner_prod(state_.model.getG(), p) - 
      0.5*inner_prod(p, prod(state_.model.getH(), p));
  
  return ared / pred;
}

double DSQA::computeTau_(const vector< double > &d,
                         const vector< double > &p)
{
  double pd = inner_prod(p, d);
  double d2 = inner_prod(d, d);
  double p2 = inner_prod(p, p);
  
  return (-pd + sqrt(pd*pd - d2*(p2 - state_.delta*state_.delta))) / d2;
}

vector< double > &DSQA::computeTrsRegStep_(const vector< double > &g,
                                           const matrix< double > &H,
                                           double delta,
                                           vector< double > &p)
{
  vector< double > r;
  vector< double > rPlus;
  vector< double > d;
  double dHd;
  double alpha, beta;
  double tau;
  vector< double > pPlus;
  int k;
  double r0_norm = norm_2(g);
  vector< double > Hd;
  double rr;
  double rPlus2;
  
  double eps = 1e-16;
  
  noalias(p) = zero_vector< double >(setup_->n);
  
  if(r0_norm < eps)
  {
    p = zero_vector< double >(setup_->n);
    return p;
  }
  
  r = g;
  d = -r;
  
  k = 0;
  while(k < setup_->n)
  {
    Hd = prod(H, d);
    dHd = inner_prod(d, Hd);
    if(dHd <= 0.0)
    {
      tau = computeTau_(d, p);
      pPlus = p + tau*d;
      break;
    }
    rr = inner_prod(r, r);
    alpha = rr / dHd;
    pPlus = p + alpha*d;
    if(norm_2(pPlus) >= delta)
    {
      tau = computeTau_(d, p);
      pPlus = p + tau*d;
      break;
    }
    rPlus = r + alpha * Hd;
    rPlus2 = inner_prod(rPlus, rPlus);
    if(rPlus2 < (eps*r0_norm)*(eps*r0_norm))
      break;
    beta = rPlus2 / rr;
    d = -rPlus + beta*d;
    r = rPlus;
    p = pPlus;
    k++;
  }
    
  p = pPlus;
  
  return p;
}

NativeSolver::IterationStatus DSQA::iterate_()
{
  double fXPlus;
  double distSq;
  vector< double > dx(setup_->n);
  bool fImproved = false;
  int j;
  double maxDistSq = 0.0;
  double ratio;
  int t = 0;
  double eta = 0.0;
  
  computeTrsRegStep_(state_.model.getG(), state_.model.getH(), state_.delta, p_);
  xPlus_ = state_.x + p_;
  fXPlus = setup_->f(xPlus_);
  ratio = computeReduction_(state_.x, xPlus_, state_.fx, fXPlus, p_);
  
  t = -1;
  maxDistSq = 0.0;
  for(j = 0; j < m_; j++)
  {
    if(j == state_.model.getLowestIndex())
      continue;
    dx = state_.model.getLowestX() - column(state_.model.getX(), j);
    distSq = inner_prod(dx, dx);
    if(distSq > maxDistSq)
    {
      maxDistSq = distSq;
      t = j;
    }
  }
  
  if(t != -1 && (fXPlus < state_.model.getLowestF() || 
     inner_prod(p_, p_) <= maxDistSq))
    fImproved = state_.model.updatePoint(xPlus_, fXPlus, t);
  
  if(ratio > eta)
    state_.delta *= 1.5;
  else
    state_.delta *= 0.75;

  if(fImproved)
    state_.model.setOrigin(state_.model.getLowestIndex());
  
  state_.x  = state_.model.getLowestX();
  state_.fx = state_.model.getLowestF();
  
  //std::cout<<"niter: "<<nIter_<<" x: "<<x_<<" dx: "<<sqrt(maxDistSq)<<std::endl;
  /*if(nIter_ % 100 == 0)
  {
    model_ = QuadInterp(objFunc_, x_, 1e-3);
    //delta_ = 1e-6;
  }*/
  
  if(state_.nIter < 10000 && state_.delta > 1e-12)
    return NativeSolver::ITERATION_CONTINUE;
  else
    return NativeSolver::ITERATION_SUCCESS;
}

void DSQA::doSetup_(const Function &objFunc,
                    const vector< double > &x0,
                    const Solver::Setup &solverSetup,
                    const Constraints &C)
{
  NativeSolver::doSetup_(objFunc, x0, solverSetup, C);
  
  state_.delta = 1e-2;
  m_ = (setup_->n+1)*(setup_->n+2)/2;
  
  state_.x = x0;
  state_.fx = setup_->f(x0);
  
  state_.model = QuadInterp(objFunc, x0, state_.delta);
  state_.model.setOrigin(state_.model.getLowestIndex());
  p_.resize(setup_->n);
  xPlus_.resize(setup_->n);
  
  p_.resize(setup_->n);
  xPlus_.resize(setup_->n);
}
