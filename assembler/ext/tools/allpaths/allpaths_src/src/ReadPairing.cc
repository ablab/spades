///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <ctype.h>
#include "CoreTools.h"
#include "ReadPairing.h"
#include "system/file/FileReader.h"

void read_pairing::ChangeId( int old_id, int new_id )
{    if ( id1 == old_id ) id1 = new_id;
     if ( id2 == old_id ) id2 = new_id;    }

bool read_pairing::InvolvesId( int id ) const
{    return id1 >= 0 && id2 >= 0 && ( id1 == id || id2 == id );    }

int read_pairing::Partner( int id ) const
{    if ( id1 != id && id2 != id )
     {    PRINT3( id, id1, id2 );
          ForceAssert( id1 == id || id2 == id );    }
     if ( id1 == id ) return id2;
     else return id1;    }

void read_pairing::Swap()
{   swap(id1, id2);    }

void read_pairing::Kill( ) { id1 = -1; }

ostream& operator<<( ostream& o, const read_pairing& p )
{    return o << p.id1 << " " << p.sep << " " << p.id2 <<
          " " << p.sd << " " << int(p.t) << " " << p.weight;    }

istream& operator>>( istream& i, read_pairing& p )
{    int t;
     i >> p.id1 >> p.sep >> p.id2 >> p.sd >> t >> p.weight;
     p.t = (read_type) t;
     return i;    }

istream& operator>>( istream& i, vec<read_pairing>& v )
{    int length;
     i >> length;
     v.resize(length);
     for ( int j = 0; j < length; j++ )
          i >> v[j];
     return i;    }

void ReadBinaryPairs( istream &in, vec<read_pairing> &vecPairs )
{   int size;
    in >> size;
    vecPairs.resize( size );
    in.get();
    if ( vecPairs.nonempty( ) )
         in.read( (char *) &vecPairs[0], 
               sizeof( read_pairing ) * vecPairs.size() );    }

void ReadPairsFile( const String &strTextFilename, vec<read_pairing> &vecPairs )
{    String strBinaryFilename( strTextFilename  + "b" );
     if ( IsRegularFile( strBinaryFilename ) &&
          IsOlder( strTextFilename, strBinaryFilename ) )
     {    BinaryRead0( strBinaryFilename, vecPairs );    }
     else READX( strTextFilename, vecPairs );    }

int PairsFileSize( const String &strTextFilename ) {
     String strBinaryFilename( strTextFilename  + "b" );
     if ( IsRegularFile( strBinaryFilename ) &&
          IsOlder( strTextFilename, strBinaryFilename ) )
       return AsciiOrBinary0VecSize( strBinaryFilename );
     else
       return AsciiOrBinary0VecSize( strTextFilename );
}

void WritePairs( const String& dir, const vec<read_pairing>& pairs, size_t nreads,
     Bool extras_only, const String& reads_name  )
{    String file_head = dir + "/" + reads_name;
     WritePairs(pairs, nreads, file_head, extras_only);   }

void WritePairs( const vec<read_pairing>& pairs, size_t nreads, const String& file_head,
		 Bool extras_only )
{    if ( !extras_only )
     {    Ofstream( outp, file_head + ".pairto" );
          outp << pairs.size( ) << "\n";
          for ( size_t i = 0; i < pairs.size(); i++ ) 
               outp << pairs[i] << "\n";    }
     BinaryWrite0( file_head + ".pairtob", pairs );
     vec<int> pairs_index( nreads, -1 );
     for ( size_t i = 0; i < pairs.size(); i++ )
     {    ForceAssertLt( static_cast<size_t>(pairs[i].id1), nreads );
          ForceAssertLt( static_cast<size_t>(pairs[i].id2), nreads );
          pairs_index[ pairs[i].id1 ] = i;
          pairs_index[ pairs[i].id2 ] = i;    }
     BinaryWrite( file_head + ".pairto_index", pairs_index );    }

void ReadsToPairs( const String& dir, const vec<int>& ids,
     vec<read_pairing>& pairs, const String& reads_name )
{    pairs.resize( ids.size( ) );
     String file_head = dir + "/" + reads_name;
     vec<int> ps;
     BinaryReadSubset( file_head + ".pairto_index", ids, ps );
     String npairs;
     {    Ifstream( in, file_head + ".pairtob" );
          in >> npairs;    }
     int head = npairs.size( ) + 1;
     FileReader fr( (file_head + ".pairtob").c_str() );
     for ( size_t i = 0; i < ids.size(); i++ )
     {    if ( ps[i] < 0 ) pairs[i].Kill( );
          else
          {    fr.seek( head + ps[i] * sizeof(read_pairing) );
               fr.read( &pairs[i], sizeof(read_pairing) );
               ForceAssert( ids[i] == pairs[i].id1 
                    || ids[i] == pairs[i].id2 );    }    }    }
