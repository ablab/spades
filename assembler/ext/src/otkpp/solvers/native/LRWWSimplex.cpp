
#include "LRWWSimplex.h"

#include <boost/numeric/ublas/matrix_proxy.hpp>

using namespace boost::numeric::ublas;

static const double gamma_ = 0.5;
static const double chi    = 2.0;
static const double rho    = 1.0;
static const double sigma  = 0.5;

LRWWSimplex::Simplex::Simplex(int n)
{
  n_ = n + 1;
  X_.resize(n, n + 1);
  fx_.resize(n + 1);
  double fac = n;
  for(int ni = n-1; ni >= 1; ni--)
    fac *= ni;
  vol_ = 1.0 / fac;
}

void LRWWSimplex::Simplex::computeCentroid()
{
  xc_ = column(X_, 0);
  for(int i = 1; i < n_-1; i++)
    xc_ += column(X_, i);
  xc_ /= n_ - 1;
}

void LRWWSimplex::Simplex::computeExpansion(vector< double > &xe)
{
  xe = (1.0+rho*chi) * xc_ - rho*chi*column(X_, n_-1);
}

void LRWWSimplex::Simplex::computeInsideContraction(vector< double > &xic)
{
  xic = (1.0 - gamma_) * xc_ + gamma_*column(X_, n_-1);
}

void LRWWSimplex::Simplex::computeOutsideContraction(vector< double > &xoc)
{
  xoc = (1.0 + rho*gamma_) * xc_ - rho*gamma_*column(X_, n_-1);
}

void LRWWSimplex::Simplex::computeReflection(vector< double > &xr)
{
  xr = (1.0 + rho) * xc_ - rho*column(X_, n_-1);
}

const vector< double > &LRWWSimplex::Simplex::getCentroidPoint() const
{
  return xc_;
}

double LRWWSimplex::Simplex::getfx(int i) const
{
  return fx_[i];
}

const matrix_column< const matrix< double > > LRWWSimplex::Simplex::getx(int i) const
{
  return column(X_, i);
}

void LRWWSimplex::Simplex::improveHighestVertex(const vector< double > &x, double fx, LRWWSimplex::Simplex::OpType opType)
{
  int i = n_ - 2;
  
  while(i >= 0 && fx < fx_[i])
  {
    fx_[i + 1] = fx_[i];
    column(X_, i+1) = column(X_, i);
    i--;
  }
  
  column(X_, i+1) = x;
  fx_[i+1] = fx;
  
  switch(opType)
  {
    case LRWWSimplex::Simplex::EXPANSION:
    vol_ *= rho * chi;
    break;
    case LRWWSimplex::Simplex::INSIDE_CONTRACTION:
    vol_ *= gamma_;
    break;
    case LRWWSimplex::Simplex::OUTSIDE_CONTRACTION:
    vol_ *= rho * gamma_;
    break;
    case LRWWSimplex::Simplex::REFLECTION:
    vol_ *= rho;
    break;
  }
}

void LRWWSimplex::Simplex::setVertex(int i, const vector< double > &x, double fx)
{
  column(X_, i) = x;
  fx_[i] = fx;
}

void LRWWSimplex::Simplex::shrink(Function &f)
{
  for(int i = 1; i < n_; i++)
  {
    column(X_, i) = column(X_, 0) + sigma * (column(X_, i) - column(X_, 0));
    fx_[i] = f(column(X_, i));
    
    vol_ *= sigma;
  }
}

void LRWWSimplex::Simplex::sortVertices()
{
  int i, j;
  double fxi;
  vector< double > tmpv;
  
  for(i = 1; i < n_; i++)
  {
    fxi = fx_[i];
    tmpv = column(X_, i);
    j = i - 1;
    
    while(j >= 0 && fx_[j] > fxi)
    {
      fx_[j + 1] = fx_[j];
      column(X_, j + 1) = column(X_, j);
      j--;
    }
    
    fx_[j + 1] = fxi;
    column(X_, j + 1) = tmpv;
  }
}

unsigned int LRWWSimplex::getM() const
{
  return setup_->n + 1;
}

std::string LRWWSimplex::getName() const
{
  return "LRWWSimplex";
}

const matrix< double > LRWWSimplex::getXArray() const
{
  matrix< double > X(setup_->n, setup_->n + 1);
  
  for(int i = 0; i < setup_->n; i++)
    for(int j = 0; j < setup_->n + 1; j++)
      X(i, j) = state_.S.getx(j)[i];
  
  return X;
}

bool LRWWSimplex::usesGradient() const
{
  return false;
}

bool LRWWSimplex::usesHessian() const
{
  return false;
}

NativeSolver::IterationStatus LRWWSimplex::iterate_()
{
  double fr;
  double fe;
  double fic;
  double foc;
  
  // Step 1: Compute centroid point.
  state_.S.computeCentroid();
  
  // Step 2: Compute reflection point.
  state_.S.computeReflection(xr_);
  fr = setup_->f(xr_);
  if(state_.S.getfx(0) <= fr && fr < state_.S.getfx(setup_->n - 1))
  {
    state_.S.improveHighestVertex(xr_, fr, LRWWSimplex::Simplex::REFLECTION);
    goto end;
  }
  
  // Step 3: Compute expansion point.
  if(fr < state_.S.getfx(0))
  {
    state_.S.computeExpansion(xe_);
    fe = setup_->f(xe_);
    
    if(fe < fr)
      state_.S.improveHighestVertex(xe_, fe, LRWWSimplex::Simplex::EXPANSION);
    else
      state_.S.improveHighestVertex(xr_, fr, LRWWSimplex::Simplex::REFLECTION);
    
    goto end;
  }
  
  // Step 4: Compute contraction points.
  if(fr >= state_.S.getfx(setup_->n - 1))
  {
    // Case 1: outside
    if(state_.S.getfx(setup_->n - 1) <= fr && fr < state_.S.getfx(setup_->n))
    {
      state_.S.computeOutsideContraction(xoc_);
      foc = setup_->f(xoc_);
      
      if(foc <= fr)
      {
        state_.S.improveHighestVertex(xoc_, foc, LRWWSimplex::Simplex::OUTSIDE_CONTRACTION);
        goto end;
      }
    }
    // Case 2: inside
    else if(fr >= state_.S.getfx(setup_->n))
    {
      state_.S.computeInsideContraction(xic_);
      fic = setup_->f(xic_);
      
      if(fic < state_.S.getfx(setup_->n))
      {
        state_.S.improveHighestVertex(xic_, fic, LRWWSimplex::Simplex::INSIDE_CONTRACTION);
        goto end;
      }
    }
  }
  
  // Step 5: Shrink if the previous steps did not yield any improvement.
  state_.S.shrink(setup_->f);
  state_.S.sortVertices();
  
  end:
  state_.x  = state_.S.getx(0);
  state_.fx = state_.S.getfx(0);
  
  return NativeSolver::ITERATION_CONTINUE;
}

void LRWWSimplex::doSetup_(const Function &objFunc,
                           const vector< double > &x0,
                           const Solver::Setup &solverSetup,
                           const Constraints &C)
{
  const int n = objFunc.getN();
  double fx;
  double fxi;
  int i;
  vector< double > xi(n);
  double stepSize = 1.0;  // TODO: read this from solverSetup
  
  NativeSolver::doSetup_(objFunc, x0, solverSetup, C);
  
  state_.S = Simplex(n);
  
  fx = objFunc(x0);
  state_.S.setVertex(0, x0, fx);
  for(i = 1; i < n + 1; i++)
  {
    xi = x0;
    xi[i - 1] += stepSize;
    fxi = objFunc(xi);
    state_.S.setVertex(i, xi, fxi);
  }
  
  state_.S.sortVertices();
  
  xe_.resize(n);
  xic_.resize(n);
  xoc_.resize(n);
  xr_.resize(n);
}
