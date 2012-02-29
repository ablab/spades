///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include "Misc.h"

char PrintBool( Bool b ) {
  if ( b != False )
    return 'R';
  else
    return 'F';
}

int BoolToInt( Bool b ) {
  return ( b ? 1 : 0 );
}
