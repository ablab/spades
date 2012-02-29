///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
// author: Filipe Ribeiro      03/2011
// 
//
//


#ifndef _MATH__INT_FREQUENCIES_H
#define _MATH__INT_FREQUENCIES_H

#include "MainTools.h"

#include "math/IntFunction.h"

class IntFrequencies : public IntFunction<size_t>
{
public:
  size_t freq(const int x) const { return (*this)[x]; }

  void to_text_file(const String & head) const
  {
    const String fn = head + ".freq";

    ofstream os;
    os.open(fn.c_str());
    os << "# 1:x  2:freq(x) 3:cum(x) 4:freq_norm(x) 5:cum_norm(x)" << endl;

    const int x0 = x_min();
    const int x1 = x_max();

    os << fixed;
    size_t total = 0;
    size_t cum = 0;
    for (int x = x0; x <= x1; x++) cum += freq(x);
    double cum_double = cum;

    cum = 0;
    for (int x = x0; x != x1; x++) {
      const size_t fx = freq(x);
      cum += fx;
      os << setw(10) << x << " "
         << setw(16) << fx << " "
         << setw(16) << cum << " " 
         << setw(16) << setprecision(12) << double(fx) / cum_double << " "
         << setw(16) << setprecision(12) << double(cum) / cum_double << endl;
    }

    os.close();
  }
};



SELF_SERIALIZABLE(IntFrequencies);




#endif
