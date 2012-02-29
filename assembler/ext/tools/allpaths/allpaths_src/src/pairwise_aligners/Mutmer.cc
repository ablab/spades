///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "pairwise_aligners/Mutmer.h"
#include "PrintAlignment.h"

void mutmer::Print( ostream &out, const basevector &rd1, const basevector &rd2 )
{
  align a;
  a.Setpos1( pos1_ );
  a.Setpos2( pos2_ );
  a.SetNblocks( 1 );
  a.SetGap( 0, 0 );
  a.SetLength( 0, len_ );

  out << *this << endl;
  PrintVisualAlignment( False, out, rd1, rd2, a );
}

