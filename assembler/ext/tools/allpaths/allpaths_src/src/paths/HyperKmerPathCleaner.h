/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_HYPERKMERPATHCLEANER_H
#define PATHS_HYPERKMERPATHCLEANER_H

#include "paths/HyperKmerPath.h"

class HyperKmerPathCleaner {
public:
  void CleanUpGraph( HyperKmerPath& ans ) const;

  void Zip( HyperKmerPath& ans ) const;

  void ZipLeft( HyperKmerPath& ans ) const;
  void ZipRight( HyperKmerPath& ans ) const;
  // Maybe this should be a HyperKmerPath member function
  bool ZipVertexRight( int vx, HyperKmerPath& ans ) const;
};


#endif
