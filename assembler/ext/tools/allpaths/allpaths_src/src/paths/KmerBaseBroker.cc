/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "paths/KmerBaseBroker.h"
#include "paths/PathEmbedding.h"

#include <set>






// Check that the indicated n bases from two basevectors are identical
bool CheckEq( basevector b1, int start1, basevector b2, int start2, int n ) {
  for( int i=0; i<n; i++ )
    if( b1[start1+i] != b2[start2+i] ) return false;
  return true;
}

template <class TAG>
const basevector& KmerBaseBrokerTemplate<TAG>::Bases(kmer_id_t k) const {

  const vecKmerPath& paths = *pathsp;
  const vecKmerPath& paths_rc = *paths_rcp;
  const vec<TAG>& pathsdb = *pathsdbp;

  basevector& kmer_bases = bases_cache[k];
  
  if( ! kmer_bases.Initialized() ) {

    kmer_bases.Setsize(0);

    path_interval_id_t match_index = Instance( pathsdb, k );
    if( match_index == -1 ) {
      cout << "Error: Bases(" << k << "): no such kmer in database!" << endl;
      return kmer_bases;
    }
  
    TAG match = pathsdb[match_index];
    longlong read_id = match.PathId();
    bool rc = (read_id < 0);
    if( rc ) read_id = -read_id-1;
    const KmerPath& the_path = ( rc ? paths_rc[read_id] : paths[read_id] );

    int offset = k - match.Start();
    for(int i=0; i<match.PathPos(); i++)
      offset += the_path.Length(i);

    if( verbose ) 
      cout << "kmer " << k << " from "
	   << (rc ? "rc of read " : "read ") << read_id
	   << "(len=" << bases[read_id].size() << ")"
	   << ", bases " << offset << "-" << offset+K-1 << endl;
    
    if( offset+K-1 > (int)bases[read_id].size() ) {
      cout << "Error: not enough bases to match path!\n";
      PRINT(read_id);
      PRINT(paths[read_id]);
      cout << "bases[read_id] = ";
      bases[read_id].Print(cout);
      cout << endl;
      return kmer_bases;
    }

    if(rc) {
      kmer_bases.SetToSubOf( bases[read_id], bases[read_id].size()-offset-K , K );
      kmer_bases.ReverseComplement();
    } else {
      kmer_bases.SetToSubOf( bases[read_id], offset, K);
    }
  }

  return kmer_bases;
}


// This asserts if a kmer in the interval is not in the database.
// (That would indicate a serious problem somewhere upstream!)
template <class TAG>
basevector KmerBaseBrokerTemplate<TAG>::Bases(KmerPathInterval rpi) const {
  
  const vecKmerPath& paths = *pathsp;
  const vecKmerPath& paths_rc = *paths_rcp;
  const vec<TAG>& pathsdb = *pathsdbp;
  ForceAssert( rpi.isSeq() );
  basevector rpi_bases( rpi.Length()+K-1 );
  int bases_done = 0;

  basevector extracted_bases;

  for(kmer_id_t k = rpi.Start(); k <= rpi.Stop(); ) {
    // Find this kmer
    path_interval_id_t match_index = Instance( pathsdb, k );
    if( match_index == -1 ) {
      cout << "\nERROR: KmerPathInterval " << rpi
	   << " contains kmer " << k
	   << ", not from any path!" << endl;
      cout << "*** Make sure you use the correct run directory.\n" <<
	" For assemblies it is under ASSEMBLIES/assemblyName/run.\n" <<
	" The reads.pathsdb.kN kmer database there has additional kmers\n" <<
	" (added to represent merges of reads overlapping by <K)\n" <<
	" that were not in the reads.pathsdb.kN in the original run directory." << endl;
      ForceAssert( match_index != -1 );
    }

    TAG match = pathsdb[match_index];
    longlong read_id = match.PathId();
    bool rc = (read_id < 0);
    if( rc ) read_id = -read_id-1;
    if ( verbose ) { PRINT2( rpi, rc ); }
    
    // Some versions of the KmerBaseBroker constructor do not load paths_rc.
    // These versions are incompatible with this function.
    if ( rc && paths_rcp == NULL )
      FatalErr( "Called KmerBaseBroker::Bases(KmerPathInterval) using a KmerBaseBroker\nwith non-reversible paths (i.e. not unipaths) and without the paths_rc loaded." );
    
    const KmerPath& the_path = ( rc ? paths_rc[read_id] : paths[read_id] );
    int offset = k - match.Start();
    for(int i=0; i<match.PathPos(); i++)
      offset += the_path.Length(i);

    // Get as many bases as we can 
    kmer_id_t kstop = min( rpi.Stop(), match.Stop() );
    int bcount = K + (kstop - k);
    
    if(rc) {
      extracted_bases.SetToSubOf(bases[read_id], 
			  bases[read_id].size()-offset-bcount, bcount );
      extracted_bases.ReverseComplement();
    } else {
      extracted_bases.SetToSubOf( bases[read_id], offset, bcount );
    }

    if( k == rpi.Start() ) {
      // beginning of interval -- copy all the bases.
      CopyBases( extracted_bases, 0, rpi_bases, 0, bcount );
      bases_done = bcount;
    } else {
      // check that the overlapping K-1 bases match
      if( ! CheckEq( extracted_bases, 0, rpi_bases, bases_done-K+1, K-1 ) ) {
	cout << "\nERROR: bad overlap of kmers " << k-1 << " and " << k << endl;
      cout << "*** Make sure you use the correct run directory.\n" <<
	" For assemblies it is under ASSEMBLIES/assemblyName/run.\n" <<
	" The reads.pathsdb.kN kmer database there has additional kmers\n" <<
	" (added to represent merges of reads overlapping by <K)\n" <<
	" that were not in the reads.pathsdb.kN in the original run directory." << endl;
      TracebackThisProcess( );
      }
      // copy the remaining ones.
      CopyBases( extracted_bases, K-1, rpi_bases, bases_done, bcount-K+1 );
      bases_done += bcount - K + 1;
    }
    
    k = kstop + 1;
  }
  return rpi_bases;
}


template <class TAG>
int KmerBaseBrokerTemplate<TAG>::MinOffset(kmer_id_t k1, kmer_id_t k2, int d_min) const {
  // Save stupid users some time:
  if( d_min >= K ) return d_min;

  pair< typename map<cache_key,int>::iterator,bool > lookup =
    min_offset_cache.insert( make_pair(cache_key(k1, k2, d_min),K) );
  
  if( !lookup.second ) { // answer is in the map
    return lookup.first->second;
  }

  // Otherwise, we need to do the search.

  const basevector& kmer1 = Bases(k1);
  const basevector& kmer2 = Bases(k2);

  // If we don't know one of the kmers, return K:
  if( kmer1.size()==0 || kmer2.size()==0 )
    return K;

  int d,i;

  for(d=max(d_min,0); d<K; d++) {
    for(i=0; i<K-d; i++)
      if (kmer1[i+d] != kmer2[i]) break;
    if (i == K-d) {
      lookup.first->second = d;
      return d;
    }
  }

  return K;
}


template <class TAG>
int KmerBaseBrokerTemplate<TAG>::MaxOffset(kmer_id_t k1, kmer_id_t k2, int d_max) const {
  if( d_max >= K ) return d_max;

  pair< typename map<cache_key,int>::iterator,bool > lookup =
    max_offset_cache.insert( make_pair(cache_key(k1, k2, d_max),-1) );
  
  if( !lookup.second ) { // answer is in the map
    return lookup.first->second;
  }

  // Otherwise, we need to do the search.

  const basevector& kmer1 = Bases(k1);
  const basevector& kmer2 = Bases(k2);

  // If we don't know one of the kmers, return the given d_max:
  if( kmer1.size()==0 || kmer2.size()==0 ) {
    lookup.first->second = d_max;
    return d_max;
  }

  int d,i;

  for(d=d_max; d>=0; d--) {
    for(i=0; i<K-d; i++)
      if (kmer1[i+d] != kmer2[i]) break;
    if (i == K-d) {
      lookup.first->second = d;
      return d;
    }
  }

  // Hmm, we found no possible overlap at all.  Let's return -1.
  return -1;
}


// Is this offset possible?
template <class TAG>
bool KmerBaseBrokerTemplate<TAG>::PossibleOffset( kmer_id_t k1, kmer_id_t k2, int d ) const {
  if( d >= K ) return true;
  const basevector& kmer1 = Bases(k1);
  const basevector& kmer2 = Bases(k2);
  for(int i=0; i<K-d; i++)
    if( kmer1[i+d] != kmer2[i] )
      return false;
  return true;
}




// Convert a KmerPath, in kmer space, to a representation of it in seq space.
template <class TAG>
SuperBaseVector KmerBaseBrokerTemplate<TAG>::ToSequence( const KmerPath& path, 
						    String name ) const {
  SuperBaseVector full_seq;
  full_seq.name = name;

  String current_seq_str, new_seq_str;

  for( int seg = 0; seg < path.NSegments(); seg++ ) {
    if( path.isGap(seg) && (seg==0 || seg==path.NSegments()-1) )
      continue; // ignore illegal initial or final gaps
    if( path.isSeq(seg) ) {
      // Handle a segment of sequence:
      if( current_seq_str.empty() ) {
	current_seq_str = Bases(path.Segment(seg)).ToString();
      } else {
	new_seq_str = Bases(path.Segment(seg)).ToString();
	if( current_seq_str.substr( current_seq_str.size()-K+1, K-1) !=
	    new_seq_str.substr( 0, K-1 ) ) {
	  cout << "\nERROR: bad overlap of adjacent intervals " 
	       << path.Segment(seg-1) << path.Segment(seg) << endl;
	  PRINT(current_seq_str);
	  PRINT(current_seq_str.substr( current_seq_str.size()-K+1, K-1));
	  PRINT(new_seq_str);
	  PRINT(new_seq_str.substr( 0, K-1 ));
	
	  cout << "*** Make sure you use the correct run directory.\n" <<
	    " For assemblies it is under ASSEMBLIES/assemblyName/run.\n" <<
	    " The reads.pathsdb.kN kmer database there has additional kmers\n" <<
	    " (added to represent merges of reads overlapping by <K)\n" <<
	    " that were not in the reads.pathsdb.kN in the original run directory." << endl;
	  TracebackThisProcess( );
	  // ERROR out!
	}
	current_seq_str += (new_seq_str.c_str() + K - 1);
      }
    }
    else { // path.isGap(seg)
      if ( path.isSeq(seg+1)          // should always be true
	   && path.Maximum(seg) < K-1 // K-1 means seq abutting, not overlaping
	   && path.Minimum(seg) == path.Maximum(seg) ) {
	// Walk right through the negative gap (if overlapping seq matches)
	int overlap = K-1 - path.Minimum(seg);
	new_seq_str = Bases(path.Segment(seg+1)).ToString();
	if( current_seq_str.substr( current_seq_str.size()-overlap, overlap) !=
	    new_seq_str.substr( 0, overlap ) ) {
	  cout << "\nERROR: bad overlap of fixed-sized gap " 
	       << path.Segment(seg-1) << path.Segment(seg)
	       << path.Segment(seg+1) << endl;
	  PRINT(overlap);
	  PRINT(current_seq_str);
	  PRINT(current_seq_str.substr( current_seq_str.size()-overlap, overlap));
	  PRINT(new_seq_str);
	  PRINT(new_seq_str.substr( 0, overlap ));
	  
	  cout << "*** Make sure you use the correct run directory.\n" <<
	    " For assemblies it is under ASSEMBLIES/assemblyName/run.\n" <<
	    " The reads.pathsdb.kN kmer database there has additional kmers\n" <<
	    " (added to represent merges of reads overlapping by <K)\n" <<
	    " that were not in the reads.pathsdb.kN in the original run directory." << endl;
	  
	  TracebackThisProcess( );
	  // ERROR out!
	}
	current_seq_str += (new_seq_str.c_str() + overlap);
	seg++;
      }
      else { // we have a gap in base space.
	basevector current_seq;
	current_seq.SetFromString( current_seq_str );
	full_seq.PushSeq( current_seq );
	full_seq.PushGap( path.Minimum(seg)-K+1, path.Maximum(seg)-K+1 );
	current_seq_str.resize(0);
      }
    }
  }
  // Final sequence:
  basevector current_seq;
  current_seq.SetFromString( current_seq_str );
  full_seq.PushSeq( current_seq );
      
  return full_seq;
  
}



	


template <class TAG>
bool KmerBaseBrokerTemplate<TAG>::KmersBetween( kmer_id_t k1, 
						kmer_id_t k2, 
						int gapsize,
						KmerPath& ans ) const {
  ans.Clear();

  const vecKmerPath& paths = *pathsp;
  const vecKmerPath& paths_rc = *paths_rcp;
  const vec<TAG>& pathsdb = *pathsdbp;

  // Only gaps of fixed size <K are determined
  if( gapsize >= K ) return false;
  // Don't bother with lookup for gaps of size zero
  if( gapsize == 0 ) return PossibleOffset( k1, k2, 1 );

  pair< typename map<cache_key,KmerPath>::iterator,bool > lookup =
    between_cache.insert( make_pair(cache_key(k1, k2, gapsize),KmerPath()) );

  if( !lookup.second ) { // answer is in the map
    ans = lookup.first->second;
    return(!ans.IsEmpty());
  }
  
  // Otherwise, we need to do the search.

  vec<path_interval_id_t> indices;
  Contains( pathsdb, k1, indices );

  for(uint i=0; i<indices.size(); i++) {
    const TAG& hit = pathsdb[indices[i]];
    const KmerPath& read = ( hit.PathId() < 0 
			     ? paths_rc[-hit.PathId()-1] 
			     : paths[hit.PathId()] );
    KmerPathLoc loc1( read, hit.PathPos() );
    loc1.SetKmer( k1 );

    KmerPathLoc loc2 = loc1;

    if( loc2.IncrementHaltAtGap( gapsize+1 ) && loc2.GetKmer() == k2 ) {
      // Hooray!  Found the stretch of kmers!
      loc1.IncrementHaltAtGap(+1);
      loc2.IncrementHaltAtGap(-1);
      read.CopySubpath(loc1, loc2, ans);
      // Store it in the map for next time
      lookup.first->second = ans;
      return true;
    }
  }

  return false;
}

template <class TAG>
basevector KmerBaseBrokerTemplate<TAG>::Seq( const KmerPath& path ) const
{   ForceAssert( path.GapFree( ) );
    return ToSequence(path).Seq(0);    }


template <class TAG>
Bool KmerBaseBrokerTemplate<TAG>::ToSequenceSafe( const KmerPath& path, 
						    String name ) const {
  SuperBaseVector full_seq;
  full_seq.name = name;

  String current_seq_str, new_seq_str;

  for( int seg = 0; seg < path.NSegments(); seg++ ) {
    if( path.isGap(seg) && (seg==0 || seg==path.NSegments()-1) )
      continue; // ignore illegal initial or final gaps
    if( path.isSeq(seg) ) {
      // Handle a segment of sequence:
      if( current_seq_str.empty() ) {
	current_seq_str = Bases(path.Segment(seg)).ToString();
      } else {
	new_seq_str = Bases(path.Segment(seg)).ToString();
	if( current_seq_str.substr( current_seq_str.size()-K+1, K-1) !=
	    new_seq_str.substr( 0, K-1 ) ) {
          return False;
	}
	current_seq_str += (new_seq_str.c_str() + K - 1);
      }
    }
    else { // path.isGap(seg)
      if ( path.isSeq(seg+1)          // should always be true
	   && path.Maximum(seg) < K-1 // K-1 means seq abutting, not overlaping
	   && path.Minimum(seg) == path.Maximum(seg) ) {
	// Walk right through the negative gap (if overlapping seq matches)
	int overlap = K-1 - path.Minimum(seg);
	new_seq_str = Bases(path.Segment(seg+1)).ToString();
	if( current_seq_str.substr( current_seq_str.size()-overlap, overlap) !=
	    new_seq_str.substr( 0, overlap ) ) {
          return False;
	}
	current_seq_str += (new_seq_str.c_str() + overlap);
	seg++;
      }
      else { // we have a gap in base space.
	basevector current_seq;
	current_seq.SetFromString( current_seq_str );
	full_seq.PushSeq( current_seq );
	full_seq.PushGap( path.Minimum(seg)-K+1, path.Maximum(seg)-K+1 );
	current_seq_str.resize(0);
      }
    }
  }
  return True;
}

template class KmerBaseBrokerTemplate<tagged_rpint>;
template class KmerBaseBrokerTemplate<big_tagged_rpint>;
