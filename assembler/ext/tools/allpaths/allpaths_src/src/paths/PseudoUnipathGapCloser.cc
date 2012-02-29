/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/PseudoUnipathGapCloser.h"


/// Returns the mean size of the gap sequences for this closer
int gapcloser::meanGapSize() const {
  if (solo())
    return gaps[0].size;
  else {
    int total = 0;
    for (int i = 0; i < gaps.isize(); ++i)
      total += gaps[0].size;
    return (total/gaps.isize());
  }
}

ostream& operator<<(ostream& out, const gapcloser& gc) {
  out << "uid1 = " << gc.uid1 << ", uid2 = " << gc.uid2 << "\n";
  for (int i = 0; i < gc.gaps.isize(); ++i)
    out << gc.gaps[i] << "\n";
  return out;
}

