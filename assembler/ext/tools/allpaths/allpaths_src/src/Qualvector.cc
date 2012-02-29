///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "FastIfstream.h"
#include "Qualvector.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< qual_t, MempoolAllocator<qual_t> >;
template class OuterVec<qvec>;
template class OuterVec<qvec,qvec::alloc_type>;
template class OuterVec< OuterVec<qvec,qvec::alloc_type>,
                  MempoolOwner<qual_t> >;

/// \file
/// \ingroup grp_quals
/// \copydoc Qualvector.h

void Print( ostream &out, const qualvector &q, const String &name,
            const int scores_per_line )
{
    out << '>' << name;
    for ( qvec::size_type i = 0; i < q.size(); ++i )
    {
        if (i % scores_per_line)
        {
            out << ' ';
        }
        else
        {
            out << '\n';
        }
        out << static_cast<unsigned int>(q[i]);
    }
    out << '\n';
}

void ReadFastaQuals( const String& fn, vecqualvector& qual,
     const vec<int>* ids_to_read )
{    int total_seqs = 0;
     longlong total_bases = 0;
     int count = 0;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) qual.Reserve( total_bases, total_seqs );
          fast_ifstream quals(fn);
          String line;
          qualvector q;
          Bool first = True;
          while(1)
          {    getline( quals, line );
               if ( quals.fail( ) )
               {    if ( pass == 1 )
                    {    total_bases += q.size( );
                         ++total_seqs;    }
                    if ( pass == 2 ) qual.push_back(q);
                    q.clear( );
                    break;    }
               ForceAssert( line.size( ) > 0 );
               if ( line[0] == '>' )
               {    if ( !first )
                    {    if ( pass == 1 )
                         {    total_bases += q.size( );
                              ++total_seqs;    }
                         if ( pass == 2 ) qual.push_back(q);    }
                    first = False;
                    q.clear( );
                    while(1)
                    {    char c;
                         quals.peek(c);
                         if ( quals.fail( ) || c == '>' ) break;
                         getline( quals, line );
                         for ( int j = 0; j < (int) line.size( ); j++ )
                              ForceAssert( isspace(line[j]) || isdigit(line[j]) );
                         istrstream i( line.c_str( ) );
                         while(1)
                         {    int n;
                              i >> n;
                              if ( i.fail( ) ) break;
                              q.push_back(n);    }    }    }
               if ( quals.fail( ) )
               {    if ( pass == 1 )
                    {    total_bases += q.size( );
                         ++total_seqs;    }
                    if ( pass == 2 )
                    {    if ( ids_to_read == 0 || BinMember( *ids_to_read, count ) )
                              qual.push_back(q);
                         ++count;    }
                    q.clear( );
                    break;    }    }    }    }
