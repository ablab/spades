/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// FetchReadsAmb: like FetchReads, but only determines if bases are ambiguous,
// and stores ambiguity information in a vecbitvector.

// If lower_case_instead=True, check instead if bases are lower case, as would be
// the case if RepeatMasker had marked the bases as repeats.

#include <ctype.h>

#include "Bitvector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"

void FetchReadsAmb( vecbitvector& b, String fasta_file, AMB_STYLE amb_style )
{
     // Scan the file to determine the total number (nseq) of bitvectors to be
     // generated, and a good upper bound (totalbases) for the total number of bases.

     int nseq = 0;
     longlong totalbases = 0;
     {    fast_ifstream in(fasta_file);
          if ( in.fail( ) ) FatalErr( "Trouble opening " << fasta_file << "." );
          String line, gt( ">" );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( gt, 0 ) ) ++nseq;
               else totalbases += line.size( );    }    }
     b.Reserve( totalbases/32 + nseq, nseq );

     if ( nseq == 0 )
     {    cout << fasta_file << " has size 0 -- I hope this is OK\n";
          return;    }

     // Now generate them.

     ifstream text( fasta_file.c_str( ) );
     char c;
     text.get(c);

     if ( c != '>' ) FatalErr( "File " << fasta_file << " is supposed to be in "
          << "fasta format.  In particular, each sequence should be prefaced by a "
          << "line which starts with >." );
     text.putback(c);
     vector<char> read;
     read.reserve(1000);
     bitvector bx;
     unsigned int i, ia = 0;
     for ( i = 0; ; i++ )
     {    if ( !text ) break;

          // Skip over comment line.

          text.get(c);
          ForceAssert( c == '>' );
          do
          {    ForceAssert(! text.fail() );
               text.get(c);    }
          while( c != '\n' );

          int read_ptr = 0;
          read.resize(0);
          while(text)
          {    text.get(c);
               if ( isspace(c) ) continue;
               if ( c == '>' )
               {    text.putback(c);
                    break;    }
               if ( !GeneralizedBase::isGeneralizedBase(c) )
                    FatalErr( "FetchReadsAmb: unrecognized character " << c << " in " << fasta_file );
               read.push_back(c);
               ++read_ptr;    }

          // Convert it to a bitvector.

          bx.resize(read_ptr);
          for ( int j = 0; j < read_ptr; j++ )
          {    if ( amb_style == AMB_EQ_ANY_AMBIGUOUS )
                  bx.Set(j,GeneralizedBase::fromChar(read[j]).isAmbiguous());
               else if ( amb_style == AMB_EQ_LOWER_CASE )
                  bx.Set(j,islower(read[j]));
               else if ( amb_style == AMB_EQ_Nn )
                  bx.Set(j, read[j] == 'N' || read[j] == 'n');
               else if ( amb_style == AMB_EQ_n )
                  bx.Set(j, read[j] == 'n');
          }
          b.push_back( bx );
          ++ia;    }

     b.resize(ia);    }
