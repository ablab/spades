//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "quake_enhanced/quake.hpp"
#include "mathfunctions.hpp"
#include <fstream>
using std::ofstream;
using std::endl;
using std::string;
using std::vector;
using quake_enhanced::Quake;

long double Quake::GetLimit(uint32_t x, long double threshold) {
  long double l = 0.5;
  long double r = 1;
  const long double eps = 1e-10;
  while (r - l > eps) {
    long double m = l + (r - l) / 2;
    long double err_prob = 0;
    long double k_prob = 0;
    for (uint32_t i = x; i < trusted_hist_.size(); ++i) {
      long double de = kK * Bernoulli(x, i, (1 - m) / kK) * trusted_hist_[i];
      long double dk = Bernoulli(x, i, m) * trusted_hist_[i];
      if (de != de || dk != dk) {
        continue;
      }
      err_prob += de;
      k_prob += dk;
    }
    //    k_prob = trusted_hist_[x];
    if (k_prob < threshold * err_prob) {
      l = m;
    } else {
      r = m;
    }
  }
  if (l + 2 * eps > 1) {
    return 1;
  }
  return l;
}

void Quake::CountLimits(long double threshold) {
  limits_.resize(25); 
  // ToDo: move to options or find a way to
  // count limits for greater numbers. Problem is in
  // precision of arithmetic operations.
  for (uint32_t i = 0; i < 25; ++i) {
    limits_[i] = GetLimit(i, threshold);
  }
}

void Quake::PrintLimits(string limits_file) {
  ofstream out(limits_file.c_str());
  for (uint32_t i = 0; i < limits_.size(); ++i) {
    out << i << " " << limits_[i] << endl;
  }
}

void Quake::PrepareLimits(long double threshold, string limits_file) {
  assert(cur_state_ >=  kTrustedHistPrepared);
  if (cur_state_ < kLimitsCounted) {
    CountLimits(threshold);
  }
  if (cur_state_ < kLimitsPrinted) {
    PrintLimits(limits_file);
  }
}
