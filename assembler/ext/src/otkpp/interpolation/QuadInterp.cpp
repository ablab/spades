
#include "QuadInterp.h"

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <iostream>

using namespace boost::numeric::ublas;

QuadInterp::QuadInterp() : n_(0), m_(0) { }

QuadInterp::QuadInterp(const Function &f, const boost::numeric::ublas::vector< double > &xb, double delta)
{
  int i;
  
  f_ = &f;
  n_ = f.getN();
  m_ = (n_+1)*(n_+2)/2;
  
  X_.resize(n_, m_);
  F_.resize(m_);
  
  g_.resize(n_);
  H_.resize(n_, n_);
  
  cl_.resize(m_);
  gl_.resize(m_);
  for(i = 0; i < m_; i++)
    gl_[i].resize(n_);
  Hl_.resize(m_);
  for(i = 0; i < m_; i++)
    Hl_[i].resize(n_, n_);
  
  cl_hat_.resize(2*n_ + 1);
  gl_hat_.resize(2*n_ + 1);
  for(i = 0; i < 2*n_+1; i++)
    gl_hat_[i].resize(n_);
  Hl_hat_.resize(2*n_ + 1);
  for(i = 0; i < 2*n_+1; i++)
    Hl_hat_[i].resize(n_, n_);
  
  initialize_(xb, delta);
}

double QuadInterp::eval(const boost::numeric::ublas::vector< double > &d)
{
  return c_ + inner_prod(g_, d) + 0.5*inner_prod(d, prod(H_, d));
}

double QuadInterp::evalLagrangian(int j, const boost::numeric::ublas::vector< double > &d)
{
  return cl_[j] + inner_prod(gl_[j], d) + 0.5*inner_prod(d, prod(Hl_[j], d));
}

double QuadInterp::getC() const
{
  return c_;
}

double QuadInterp::getCl(int i) const
{
  return cl_[i];
}

const std::vector< double > &QuadInterp::getF() const
{
  return F_;
}

const boost::numeric::ublas::vector< double > &QuadInterp::getG() const
{
  return g_;
}

const boost::numeric::ublas::vector< double > &QuadInterp::getGl(int i) const
{
  return gl_[i];
}

const matrix< double > &QuadInterp::getH() const
{
  return H_;
}

const matrix< double > &QuadInterp::getHl(int i) const
{
  return Hl_[i];
}

double QuadInterp::getLowestF() const
{
  return F_[xiLowest_];
}

int QuadInterp::getLowestIndex() const
{
  return xiLowest_;
}

const matrix_column< const matrix< double > > QuadInterp::getLowestX() const
{
  return column(X_, xiLowest_);
}

const boost::numeric::ublas::vector< double > &QuadInterp::getOrigin() const
{
  return xb_;
}

const matrix< double > &QuadInterp::getX() const
{
  return X_;
}

void QuadInterp::setOrigin(int xbi)
{
  const boost::numeric::ublas::vector< double > &xb = column(X_, xbi);
  boost::numeric::ublas::vector< double > dx = xb - xb_;
  
  for(int i = 0; i < m_; i++)
  {
    cl_[i] += 0.5*inner_prod(dx, prod(Hl_[i], dx)) + inner_prod(dx, gl_[i]);
    gl_[i] += prod(Hl_[i], dx);
  }
  
  //c_ += 0.5*inner_prod(dx, prod(H_, dx)) + inner_prod(dx, g_);
  c_ = F_[xbi];
  g_ += prod(H_, dx);
  
  xb_ = xb;
}

void QuadInterp::test()
{
  std::cout<<"Initial model:"<<std::endl;
  printInfo_();
  
  std::cout<<"Model after update:"<<std::endl;
  boost::numeric::ublas::vector< double > xNew = column(X_, 2);
  xNew[0] += 0.1;
  xNew[1] += -0.2;
  updatePoint(xNew, 2);
  printInfo_();
  
  std::cout<<"Model after origin shift:"<<std::endl;
  xNew = xb_;
  xNew[0] -= 0.1;
  xNew[1] += 0.2;
  setOrigin(getLowestIndex());
  printInfo_();
}

void QuadInterp::testInvariants()
{
  int i;
  double lix, Lx = 0.0;
  boost::numeric::ublas::vector< double > x;
  
  for(i = 0; i < m_; i++)
  {
    if((*f_)(column(X_, i)) != F_[i])
      throw std::runtime_error("invalid function value");
    if((*f_)(column(X_, i)) < F_[xiLowest_])
      throw std::runtime_error("invalid best point");
    if(fabs(eval(column(X_, i) - xb_) - F_[i]) > 1e-3)
    {
      std::cout<<"Q(x): "<<eval(column(X_, i) - xb_)<<std::endl;
      std::cout<<"f(x): "<<F_[i]<<std::endl;
      throw std::runtime_error("invalid interpolation");
    }
  }
  
  x = xb_;
  for(i = 0; i < n_; i++)
    x[i] += (rand() % 100) / 100.0;
  for(i = 0; i < m_; i++)
  {
    lix = evalLagrangian(i, x-xb_);
    Lx += F_[i] * lix;
  }
  if(fabs(Lx - eval(x-xb_)) > 1e-3)
  {
    std::cout<<"L(x): "<<Lx<<std::endl;
    std::cout<<"Q(x): "<<eval(x-xb_)<<std::endl;
    throw std::runtime_error("mismatching lagrange functions and quadratic model");
  }
}

bool QuadInterp::updatePoint(const boost::numeric::ublas::vector< double > &x, int j)
{
  return updatePoint(x, NAN, j);
}

bool QuadInterp::updatePoint(const boost::numeric::ublas::vector< double > &x, double fx, int j)
{
  // Update the Lagrange coefficients.
  double c1, c2;
  boost::numeric::ublas::vector< double > dx = x - xb_;
  int i;
  bool improved = false;
  double m;
  
  if(std::isnan(fx))
    fx = (*f_)(x);
  
  c2 = evalLagrangian(j, dx);
  
  cl_[j] /= c2;
  gl_[j] /= c2;
  Hl_[j] /= c2;
  
  for(i = 0; i < m_; i++)
  {
    if(i == j)
      continue;
    
    c1 = evalLagrangian(i, dx);
    
    cl_[i] -= c1*cl_[j];
    gl_[i] -= c1*gl_[j];
    Hl_[i] -= c1*Hl_[j];
  }
  
  if(fx < F_[xiLowest_])
  {
    xiLowest_ = j;
    improved = true;
  }
  
  // Update the quadratic model.
  m = fx - eval(dx);
  c_ += m*cl_[j];
  g_ += m*gl_[j];
  H_ += m*Hl_[j];
  
  column(X_, j) = x;
  F_[j] = fx;
  
  return improved;
}

void QuadInterp::initialize_(const boost::numeric::ublas::vector< double > &xb, double delta)
{
  double D;
  double gammaj, gammak;
  int i, j, k;
  double m1, m2, m3, m4;
  double r1, r2;
  boost::numeric::ublas::vector< double > alpha(n_);
  boost::numeric::ublas::vector< double > beta(n_);
  
  xb_ = xb;
  
  c_ = (*f_)(xb);
  
  column(X_, 0) = xb;
  for(j = 0; j < n_; j++)
  {
    column(X_, 2*j+1) = xb;
    alpha[j] = delta;
    column(X_, 2*j+1)[j] += alpha[j];
    
    if((*f_)(column(X_, 2*j+1)) < c_)
      beta[j] = 2.0*delta;
    else
      beta[j] = -2.0*delta;
    
    column(X_, 2*j+2) = xb;
    column(X_, 2*j+2)[j] += beta[j];
    
    m1 = 0.5*delta*delta;
    m2 = delta;
    m3 = 0.5*beta[j]*beta[j];
    m4 = beta[j];
    
    r1 = (*f_)(column(X_, 2*j+1)) - c_;
    r2 = (*f_)(column(X_, 2*j+2)) - c_;
    D = m1*m4 - m2*m3;
    H_(j, j) = (m4*r1 - m2*r2) / D;
    g_(j) = (-m3*r1 + m1*r2) / D;
  }
  
  for(j = 0; j < n_; j++)
  {
    for(k = j+1; k < n_; k++)
    {
      i = 2*n_+1+j+1+k*(k-1)/2-1;
      column(X_, i) = xb;
      if((*f_)(column(X_, 2*j+1)) < c_)
        gammaj = delta;
      else
        gammaj = -delta;
      if((*f_)(column(X_, 2*k+1)) < c_)
        gammak = delta;
      else
        gammak = -delta;
      
      column(X_, i)[j] += gammaj;
      column(X_, i)[k] += gammak;
      
      H_(j, k) = ((*f_)(column(X_, i)) - 0.5*gammak*gammak*H_(k, k) - 
                  0.5*gammaj*gammaj*H_(j, j) - gammaj*g_(j) - gammak*g_(k) - c_) / 
                  (gammaj*gammak);
      H_(k, j) = H_(j, k);
      
      computeLagrangeCoeff_last_(xb, gammaj*gammak, i, j, k);
    }
  }
  
  double denom1, denom2;
  for(j = 0; j < n_; j++)
  {
    denom1 = alpha[j]*alpha[j] - alpha[j]*beta[j];
    denom2 = beta[j]*beta[j] - alpha[j]*beta[j];
    
    computeLagrangeCoeff_hat_(xb, denom1, denom2, j);
  }
  
  for(i = 1; i < 2*n_+1; i++)
    computeLagrangeCoeff_first_(i, xb);
  computeLagrangeCoeff_first_(0, xb);
  
  for(i = 0; i < m_; i++)
    F_[i] = (*f_)(column(X_, i));
  
  double fLowest = F_[0];
  xiLowest_ = 0;
  for(i = 1; i < m_; i++)
  {
    if(F_[i] < fLowest)
    {
      xiLowest_ = i;
      fLowest = F_[i];
    }
  }
}

void QuadInterp::computeLagrangeCoeff_first_(int i, const boost::numeric::ublas::vector< double > &xb)
{
  if(i > 0)
  {
    double coeff;
    boost::numeric::ublas::vector< double > dx;
    
    cl_[i] = cl_hat_[i];
    gl_[i] = gl_hat_[i];
    Hl_[i] = Hl_hat_[i];
    
    for(int t = 2*n_+1; t < m_; t++)
    {
      dx = column(X_, t) - xb;
      coeff = 0.5*inner_prod(prod(dx, Hl_hat_[i]), dx) + 
          inner_prod(gl_hat_[i], dx) + 
          cl_hat_[i];
      cl_[i] -= coeff * cl_[t];
      gl_[i] -= coeff * gl_[t];
      Hl_[i] -= coeff * Hl_[t];
    }
  }
  else
  {
    cl_[0] = 1.0;
    gl_[0] = zero_vector< double >(n_);
    Hl_[0] = zero_matrix< double >(n_, n_);
    
    for(int t = 1; t < m_; t++)
    {
      cl_[0] -= cl_[t];
      gl_[0] -= gl_[t];
      Hl_[0] -= Hl_[t];
    }
  }
}

void QuadInterp::computeLagrangeCoeff_last_(const boost::numeric::ublas::vector< double > &xb, double denom, int i, int j, int k)
{
  Hl_[i] = zero_matrix< double >(n_, n_);
  Hl_[i](j, k) = Hl_[i](k, j) = 1.0 / denom;
  gl_[i] = zero_vector< double >(n_);
  cl_[i] = 0.0;
}

void QuadInterp::computeLagrangeCoeff_hat_(const boost::numeric::ublas::vector< double > &xb, double denom1, double denom2, int j)
{
  Hl_hat_[2*j+1] = zero_matrix< double >(n_, n_);
  Hl_hat_[2*j+1](j, j) = 2.0 / denom1;
  gl_hat_[2*j+1] = zero_vector< double >(n_);
  gl_hat_[2*j+1][j] = (-column(X_, 2*j+2)[j] + xb[j]) / denom1;
  cl_hat_[2*j+1] = 0.0;
  
  Hl_hat_[2*j+2] = zero_matrix< double >(n_, n_);
  Hl_hat_[2*j+2](j, j) = 2.0 / denom2;
  gl_hat_[2*j+2] = zero_vector< double >(n_);
  gl_hat_[2*j+2][j] = (-column(X_, 2*j+1)[j] + xb[j]) / denom2;
  cl_hat_[2*j+2] = 0.0;
}

void QuadInterp::printInfo_()
{
  double Lx;
  double ljx;
  boost::numeric::ublas::vector< double > x;
  
  for(int i = 0; i < m_; i++)
  {
    Lx = 0;
    x = column(X_, i);
    std::cout<<"X["<<i+1<<"]:"<<std::endl;
    std::cout<<"f(x)="<<(*f_)(x)<<std::endl;
    std::cout<<"Q(x)="<<c_ + inner_prod(g_, x-xb_) + 0.5 * inner_prod(x-xb_, prod(H_, x-xb_))<<std::endl;
    for(int j = 0; j < m_; j++)
    {
      /*std::cout<<"clj: "<<cl_[j]<<std::endl;
      std::cout<<"glj: "<<gl_[j]<<std::endl;
      std::cout<<"Hlj: "<<Hl_[j]<<std::endl;*/
      ljx = cl_[j] + inner_prod(gl_[j], x-xb_) + 0.5*inner_prod(x-xb_, prod(Hl_[j], x-xb_));
      std::cout<<"lj[x]: "<<ljx<<std::endl;
      Lx += (*f_)(column(X_, j)) * ljx;
    }
    std::cout<<"L(x)="<<Lx<<std::endl;
  }
}
