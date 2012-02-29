///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef QUERY_LOOKUP_TABLE_CORE_H
#define QUERY_LOOKUP_TABLE_CORE_H

#include "CoreTools.h"
#include "TokenizeString.h"

void QueryLookupTableCore( int argc, char *argv[] );

inline void QueryLookupTableCore( const String& args )
{    vec<String> tokens;
     Tokenize( "QueryLookupTable " + args, tokens );
     int argc = tokens.size( );
     vec<char*> argv(argc);
     for ( int i = 0; i < argc; i++ )
       argv[i] = const_cast<char*>(tokens[i].c_str());

     QueryLookupTableCore( argc, &argv[0] );    }

#endif
