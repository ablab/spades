/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "system/Utils.h"
#include <iostream>

void StrongWarning( String msg ) {
  using std::cout;
  using std::endl;
  cout << "\n\n*************************************************************\n";
  cout <<     "*** WARNING: " << msg << endl;
  cout << "*************************************************************\n\n";
}
