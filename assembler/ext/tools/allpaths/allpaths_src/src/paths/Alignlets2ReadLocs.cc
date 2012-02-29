///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include "paths/Alignlets2ReadLocs.h"
#include "paths/ReadLoc.h"

/**
 * Alignlets2ReadLocs
 */
void Alignlets2ReadLocs( const PairsManager &pairs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 vec<read_loc> &locs,
			 ostream *log )
{
  // WARNING - reads will be all tagged as type 1 (jump).

  locs.clear( );
  
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  // Map reads onto target.
  out << Date( ) << ": Alignlets2ReadLocs - mapping reads to target" << endl;
  vec< triple<int,int,int64_t> > helper;   // cg_id, start_on_cg, read_id
  helper.reserve( aligns.size( ) );
  for (int64_t read_id=0; read_id<(int64_t)index.size( ); read_id++) {
    if ( index[read_id] < 0 ) continue;
    int cid = aligns[ index[read_id] ].TargetId( );
    int start = aligns[ index[read_id] ].pos2( );
    helper.push_back( triple<int,int,int64_t>( cid, start, read_id ) );
  }
  sort( helper.begin( ), helper.end( ) );
  
  // Reserve memory.
  locs.reserve( helper.size( ) );

  // Generate vector of read_locs.
  out << Date( ) << ": Alignlets2ReadLocs - generating read_locs" << endl;
  for (size_t ii=0; ii<helper.size( ); ii++) {
    const alignlet &al = aligns[ index[ helper[ii].third ] ];
    int64_t rid = helper[ii].third;
    int64_t pairid = pairs.getPairID( rid );
    int64_t p_id =
      pairs.ID1( pairid ) == rid
      ? pairs.ID2( pairid )
      : pairs.ID1( pairid );
    if ( index[p_id] < 0 ) continue;   // partner not aligned
    const alignlet &p_al = aligns[ index[p_id] ];
    int32_t cid = al.TargetId( );
    int32_t p_cid = p_al.TargetId( );
    int cpos = al.pos2( );
    int p_cpos = p_al.pos2( );
    Bool fw = al.Fw1( );
    Bool p_fw = p_al.Fw1( );
    int8_t rclass = 1;
    uint8_t lib_id = pairs.libraryID( pairid );
    uint16_t rlen = al.Pos2( ) - al.pos2( );
    uint16_t p_rlen = p_al.Pos2( ) - p_al.pos2( );
    uint8_t band = 0;    // bandwidth - undefined!
    uint8_t p_band = 0;    // bandwidth - undefined!
    bool p_placed = true;

    locs.push_back( read_loc( rid, p_id, cid, p_cid, cpos, p_cpos,
			      fw, p_fw, rclass, lib_id, rlen, p_rlen,
			      band, p_band, p_placed ) );
  }
  
  out << Date( ) << ": Alignlets2ReadLocs - sorting read_locs" << endl;  
  sort( locs.begin( ), locs.end( ) );

  // Done.
  out << Date( ) << ": Alignlets2ReadLocs - done" << endl;
  
}

/**
 * LoadReadLocs
 */
void LoadReadLocs( const String &pairs_file,
		   const String &head_file,
		   vec<read_loc> &locs,
		   ostream *log )
{
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  // File names.
  String aligns_file = head_file;
  String index_file = head_file + ".index";

  // Load data.
  out << Date( ) << ": LoadReadLocs - loading pairs" << endl;
  PairsManager pairs( pairs_file );

  out << Date( ) << ": LoadReadLocs - loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryRead3( aligns_file, aligns );
  
  vec<int> index;
  BinaryRead3( index_file, index );
  
  // Convert (logging within).
  Alignlets2ReadLocs( pairs, aligns, index, locs, log );
  
  // Done.
  out << Date( ) << ": LoadReadLocs - done" << endl;

}
