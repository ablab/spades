/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "graphics/Eps.h"

void PrintEpsHeader( ostream& out, 
                     const float horizSize, 
                     const float vertSize,
                     const float border )
{
  const int leftBoundbox = 0;
  const int rightBoundbox = border + horizSize + border;
  const int bottomBoundbox = 0;
  const int topBoundbox = border + vertSize + border;

  out.setf( ios::fixed );
  
  out << "%!PS-Adobe-3.0 EPSF-3.0\n"
      << "%%BoundingBox: "
      << leftBoundbox << " " << bottomBoundbox << " "
      << rightBoundbox << " " << topBoundbox
      << "\n";

  // Translate boundbox to origin.
  out << border << " " <<  border << " translate\n";
}
