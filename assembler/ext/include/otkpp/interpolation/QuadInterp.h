
#ifndef QUADINTERP_H

#include <otkpp/objfunc/Function.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <vector>

/// Represents a quadratic model interpolating the objective function.
class QuadInterp
{
 public:
  QuadInterp();
  QuadInterp(const Function &f, const boost::numeric::ublas::vector< double > &xb, double delta);
  
  double eval(const boost::numeric::ublas::vector< double > &d);
  double evalLagrangian(int j, const boost::numeric::ublas::vector< double > &d);
  double getC() const;
  double getCl(int i) const;
  const std::vector< double > &getF() const;
  const boost::numeric::ublas::vector< double > &getG() const;
  const boost::numeric::ublas::vector< double > &getGl(int i) const;
  const boost::numeric::ublas::matrix< double > &getH() const;
  const boost::numeric::ublas::matrix< double > &getHl(int i) const;
  double getLowestF() const;
  int getLowestIndex() const;
  const boost::numeric::ublas::matrix_column< const boost::numeric::ublas::matrix< double > > getLowestX() const;
  const boost::numeric::ublas::vector< double > &getOrigin() const;
  const boost::numeric::ublas::matrix< double > &getX() const;
  void setOrigin(int xbi);
  void test();
  void testInvariants();
  bool updatePoint(const boost::numeric::ublas::vector< double > &x, int j);
  bool updatePoint(const boost::numeric::ublas::vector< double > &x, double fx, int j);
 private:
  const Function *f_;
  int n_, m_;
  
  boost::numeric::ublas::matrix< double > X_;
  std::vector< double > F_;
  
  boost::numeric::ublas::vector< double > xb_;
  
  double c_;
  boost::numeric::ublas::vector< double > g_;
  boost::numeric::ublas::matrix< double > H_;
  
  std::vector< double > cl_;
  std::vector< boost::numeric::ublas::vector< double > > gl_;
  std::vector< boost::numeric::ublas::matrix< double > > Hl_;

  std::vector< double > cl_hat_;
  std::vector< boost::numeric::ublas::vector< double > > gl_hat_;
  std::vector< boost::numeric::ublas::matrix< double > > Hl_hat_;
  
  int xiLowest_;
  
  void initialize_(const boost::numeric::ublas::vector< double > &xb, double delta);
  void computeLagrangeCoeff_first_(int i, const boost::numeric::ublas::vector< double > &x0);
  void computeLagrangeCoeff_last_(const boost::numeric::ublas::vector< double > &x0, double denom, int i, int j, int k);
  void computeLagrangeCoeff_hat_(const boost::numeric::ublas::vector< double > &x0, double denom1, double denom2, int j);
  void printInfo_();
};

#define QUADINTERP_H

#endif
