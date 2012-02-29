///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

const char *DOC =
  "A quick and simple program to view the contents of a PairsManager "
  "binary file (.pairs).";
 
#include "MainTools.h"
#include "PairsManager.h"


int main( int argc, char *argv[] )
{
  RunTime();
  
  BeginCommandArguments;
  CommandDoc(DOC);
  CommandArgument_String_Doc( PAIRS, "PairsManager file (ends in .pairs)" );
  EndCommandArguments;
  
  // Be a little flexible with the filename
  PAIRS = PAIRS.SafeBefore(".fastb");
  PAIRS = PAIRS.SafeBefore(".pairs");
  if (PAIRS.EndsWith("."))
    PAIRS = PAIRS.RevBefore(".");

  PairsManager pairs( PAIRS + ".pairs" );
  pairs.print( cout );
}
