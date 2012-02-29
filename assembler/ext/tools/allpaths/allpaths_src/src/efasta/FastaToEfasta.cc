///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FastaToEfasta: convert a file of nucleotides from fasta to efasta format, by
//
// (a) raising acgt to ACGT
// (b) expanding fasta ambiguity codes, including N
// (c) raise n to N and treat as gap character
// (d) sanitizing header lines to match the requirements of efasta.
//
// This code also folds lines to a max of 80 characters.  Long header lines are
// truncated.
//
// Note the very confusing role of N and n.

#include "FastIfstream.h"
#include "MainTools.h"
#include "TokenizeString.h"
#include "efasta/EfastaTools.h"

#define Err(message)                                      \
{    cout << message << endl << "\nInvalid.\n" << endl;   \
     return 1;    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(IN);
     CommandArgument_String(OUT);
     EndCommandArguments;

     Ofstream( out, OUT );
     fast_ifstream in(IN);
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) 
          {    if ( line.size( ) == 0 ) Err( "Illegal empty file." );
               break;    }
          if ( line.size( ) == 0 ) continue;
          if ( line[0] != '>' )
          {    Err( "See line = '" << line << "', which was expected "
                    << "to start with >." );    }
          for ( int i = 0; i < line.isize( ); i++ )
          {    char x = line[i];
               if ( x < 33 && x > 126 && x != ' ' && x != '\t' ) line[i] = '%';    }
          if ( line.size( ) < 2 || line[1] == ' ' || line[1] == '\t' )
               Err( "Illegal header line '" << line << "'." );
          if ( line.size( ) > 80 )
          {    for ( int j = 0; j < 76; j++ )
                    out << line[j];
               out << " ...\n";    }
          else out << line << "\n";
          vec<String> lines;
          Bool eof = False;
          while(1)
          {    char c;
               in.peek(c);
               if ( in.fail( ) ) { eof = True; break; }
               if ( c == '>' ) break;
               getline( in, line );
               lines.push_back(line);    }
          if ( lines.empty( ) ) Err( "Illegal record of empty length." );
          vec<char> all;
          int64_t all_size = 0;
          for ( size_t i = 0; i < lines.size( ); i++ )
               all_size += lines[i].size( );
          all.reserve(all_size);
          for ( size_t i = 0; i < lines.size( ); i++ )
          {    for ( size_t j = 0; j < lines[i].size( ); j++ )
               {    String s = ExpandAmbCode( lines[i][j] );
                    for ( int k = 0; k < s.isize( ); k++ )
                         all.push_back( s[k] );    }
               all.push_back( '\n' );    }
          int count = 0;
          for ( size_t i = 0; i < all.size( ); i++ )
          {    if ( all[i] == '\n' )
               {    if ( count > 0 )
                    {    out << all[i];
                         count = 0;    }    }
               else
               {    out << all[i];
                    count++;
                    if ( count == 80 )
                    {    out << "\n";
                         count = 0;    }    }    }
          if (eof) break;    }    }
