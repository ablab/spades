// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// Qualb2Quala: convert a file reads.qualb into a human-readable reads.qual.

#include "MainTools.h"
#include "math/Functions.h"
#include "Qualvector.h"

int main( int argc, char *argv[] )
{
     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     EndCommandArguments;

     vecqualvector Q;
     Q.ReadAll( IN );

     Ofstream( out, OUT);
     for ( vecqvec::size_type id = 0; id < Q.size( ); id++ ) {
          out << ">sequence_" << id << "\n";
          for ( qvec::size_type j = 0; j < Q[id].size( ); ++j ) {
               if ( j > 0 ) out << " ";
               out << (int) Q[id][j];
          }
          out << "\n";
     }
     out.close( );
}
