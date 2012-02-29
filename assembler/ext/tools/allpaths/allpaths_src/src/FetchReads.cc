///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/// Read n sequences from a file in fasta format.
/// \fn FetchReads
/// \file FetchReads.h
/// The file is supposed to consist of  >= n records of the following form:
/// > ... (newline), followed by a string from the alphabet {A,C,G,T,a,c,g,t},
/// possibly interspersed with white space.  We also allow various ambiguous
/// bases (as defined by ambiguous_base in FetchReads.h); these are converted
/// to random bases.
///
/// If n = 0, return all reads in the file.
///
/// The vector q is set to a vector of quality scores for the bases, computed
/// according to the following rule: if the base in the file is ambiguous,
/// the score is 1, and otherwise it is 25.  (This step is ignored if
/// no_q = True.)
///
/// If min_size is set to a positive value, ignore reads of length exceeding
/// it, and log correlation information.
///
/// If the "ids_to_read" argument is provided, it is to be a sorted list of
/// the sequence ids to be returned.

#include "Basevector.h"
#include "CoreTools.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "Qualvector.h"
#include "random/Random.h"

// Heuristic constants:

namespace
{
const int HighQuality = 25;
const int LowQuality = 1;
}

void FetchReads( vecbasevector& b, vecqualvector& q, unsigned int n,
                 String fasta_file, int amb_break, int min_size,
                 ostream& out, Bool no_q,
                 const vec<int>* ids_to_read )
{
  FetchReads(b, q, 0, n, fasta_file, amb_break, min_size, out,  no_q,
             ids_to_read );
}

void FetchReads( vecbasevector& b, vecqualvector& q, vecString * names,
                 unsigned int n,  String fasta_file,
                 int amb_break, int min_size, ostream& out, Bool no_q,
                 const vec<int>* ids_to_read )
{
  // Check for unsupported arguments.

  if ( amb_break > 0 )
    FatalErr( "amb_break option for FetchReads is no longer supported." );

  // Initialize random number generator.  The reason for doing this is that if
  // FetchReads is called twice on the same file, and that file has ambiguous bases,
  // we want to get the same answer.

  srandomx(1643295582);

  // Scan the file to determine the total number (nseq) of basevectors to be
  // generated, and a good upper bound (totalbases) for the total number of bases.

  int nseq = 0;
  longlong totalbases = 0;
  {
    fast_ifstream in(fasta_file);
    if ( in.fail( ) ) FatalErr( "Trouble opening " << fasta_file << "." );
    String line, gt( ">" );
    while(1) {
      getline( in, line );
      if ( line.Contains( gt, 0 ) ) {
	if (names && ( ids_to_read == 0 || BinMember( *ids_to_read, nseq ) )) {
	  names->push_back_reserve(line.After(gt));
	}
	++nseq;
      }
      else totalbases += line.size( );
      if ( in.fail( ) ) break;
      if ( n > 0 && nseq == (int) n ) break;
    }
  }
  b.clear( );
  b.Reserve( totalbases/16 + nseq, nseq );
  if ( !no_q ) {
    q.clear();
    q.Reserve( totalbases, nseq );
  }
  if ( nseq == 0 ) {
    cout << fasta_file << " has size 0 -- I hope this is OK\n";
    return;
  }

  // Now generate them.

  int count = 0;
  ifstream text( fasta_file.c_str( ) );
  int line = 1;
  char c;
  text.get(c);

  if ( c != '>' ) FatalErr( "File " << fasta_file << " is supposed to be in "
                            << "fasta format.  In particular, each sequence should be prefaced by a "
                            << "line which starts with >.  Problem at line 1." );
  text.putback(c);
  vector<char> read;
  read.reserve(1000);
  basevector bx;
  qualvector qx;
  unsigned int i, ia = 0;
  for ( i = 0; n == 0 || i < n; i++ ) {
    if ( !text ) {
      if ( n > 0 ) {
	cerr << "FetchReads: failed on read " << i+1 << "\n";
	cerr << "Check n (argument 2) " << "\n";
	exit(1);
      }
      break;
    }

    // Skip over comment line.
    text.get(c);
    ForceAssert( c == '>' );
    do {
      ForceAssert(! text.fail() );
      text.get(c);
    } while( c != '\n' );
    ++line;

    int read_ptr = 0;
    read.resize(0);
    while ( true )
    {
        text.get(c);
        if ( text.fail() )
            break;
        if ( c == '\n' )
            ++line;
        if ( isspace(c) )
            continue;
        if ( c == '>' )
        {
            text.putback(c);
            break;
        }
        if ( !GeneralizedBase::isGeneralizedBase(c) )
            FatalErr( "FetchReads: unrecognized character " << c << " in " << fasta_file << ", at line " << line << "." );
        read.push_back(c);
        ++read_ptr;
    }

    // Convert it to a basevector.

    if ( read_ptr >= min_size )
    {
        if ( min_size > 0 )
            out << i << " --> " << ia << "\n";

        bx.Setsize(read_ptr);
        if ( !no_q )
            qx.resize(read_ptr);

        for ( int j = 0; j < read_ptr; j++ )
        {
            qual_t qVal = HighQuality;
            GeneralizedBase const& base = GeneralizedBase::fromChar(read[j]);
            if ( !base.isAmbiguous() )
                bx.Set(j,static_cast<Base const&>(base).val());
            else
            {
                bx.Set(j,base.random());
                qVal = LowQuality;
            }
            if ( !no_q )
                qx[j] = qVal;
        }
        if ( ids_to_read == 0 || BinMember(*ids_to_read, count) )
        {
            b.push_back(bx);
            if ( !no_q )
                q.push_back(qx);
            ++ia;
        }
    }
    ++count;
    }

  b.resize(ia);
}

void FetchReads( vecbasevector& b, unsigned int n, String fasta_file,
                 int amb_break, int min_size, ostream& out,
		 const vec<int>* ids_to_read )
{
  vecqualvector q;
  FetchReads( b, q, n, fasta_file, amb_break, min_size, out, True,
	      ids_to_read );
}

void FetchReads( VecFloatVec& f, const String& fasta_file )
{    f.clear( );
     String line;
     FloatVec x;
     fast_ifstream in(fasta_file);
     getline( in, line );
     ForceAssert( line.Contains( ">", 0 ) );
     while(1)
     {    getline( in, line );
          if ( in.fail( ) )
          {    f.push_back_reserve(x);
               break;    }
          if ( line.Contains( ">", 0 ) )
          {    f.push_back_reserve(x);
               x.clear( );    }
          else
          {   istringstream iline( line.c_str( ) );
              float a;
              while ( iline >> a )
                  x.push_back(a);
          }
     }
}
