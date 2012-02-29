// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

// This file defines the "assembly" class, which was intended to encapsulate
// data from an assembly (and various indices in it) to facilitate analysis
// an operations on that assembly.  However, we don't like the structure much,
// and are likely to redesign it from scratch at some point in the future.

#ifndef AN_ASSEMBLY_CLASS
#define AN_ASSEMBLY_CLASS

#include <algorithm>
#include <map>
#include <set>

#include "Basevector.h"
#include "BriefAlign.h"
#include "CoreTools.h"
#include "Qualvector.h"
#include "ReadLocation.h"
#include "ReadPairing.h"

class assembly;

class super {

     public:

     vec<int> mtig;
     vec<int> gap;    // there are mtig.size( ) - 1 of these

     int Length( const assembly& A ) const;
     int ReducedLength( const assembly& A ) const;

};

class assembly {

     public:

     // core data (mutable):

     vecbasevector mtig;                     // the mtigs
     vecqualvector mtig_qual;                // and their quality scores
     vec<read_location> reads, reads_orig;   // read placement relative to mtigs
     vec<super> supers;                      // how the mtigs fit into supers
     vec<brief_align> all_aligns;            // good aligns between original reads
     vec<brief_align> mtig_aligns;           // known alignments between mtigs
     vec<Bool> mtig_aligns_lost;

     // core data (mostly const):

     int N;                                  // number of original reads
     vec<read_pairing> pairs;                // pairings between reads
     vec<int> orig_read_lengths;

     // indices to facilitate access:

     vec< vec<int> > reads_index, reads_orig_index;
     vec<int> simple_reads_orig_index; // read id --> index in reads_orig
     vec<int> pairs_index;
     map<int, int> mtigs_to_supers, mtigs_to_super_pos;
     vec<int> all_aligns_index;
     vec<int> mtig_aligns_index;

     // elogp should point to an open contig event log stream
     ostream *elogp;
     
     void SetMtigAligns( const vec<brief_align>& b )
     {    mtig_aligns = b;
          mtig_aligns_index.resize( mtig.size( ) );
          for ( size_t i = 0; i < mtig.size( ); i++ )
               mtig_aligns_index[i] = -1;
          for ( int i = (int) mtig_aligns.size( ) - 1; i >= 0; i-- )
               mtig_aligns_index[ mtig_aligns[i].id1 ] = i;
          mtig_aligns_lost.resize( mtig.size( ) );
          for ( size_t i = 0; i < mtig.size( ); i++ )
               mtig_aligns_lost[i] = False;    }

     // Constructor

     assembly() : elogp(0) { }

     assembly( String fasta_file, String qual_file,
               const vec<read_location>& reads, const vec<read_location>& reads_orig,
               const vec<super>& supers, const vec<brief_align>& all_aligns, 
               int N, const vec<read_pairing>& pairs, 
               const vec<int>& orig_read_lengths, const vec<int>& all_aligns_index,
               Bool store_contigs_one_by_one = False );

     // Number of mtigs in a supercontig

     int SuperSize( int s ) const { return supers[s].mtig.size( ); }

     // Length of a contig or a supercontig

     int Len( int m ) const { return mtig[m].size( ); }

     int SuperLen( int s ) const
     {    const super& sup = supers[s];
          int answer = 0;
          for ( int i = 0; i < (int) sup.mtig.size( ); i++ )
          {    answer += mtig[ sup.mtig[i] ].size( );
               if ( i < (int) sup.mtig.size( ) - 1 ) answer += sup.gap[i];    }
          return answer;   }

     int SuperLenWithoutGaps( int s ) const
     {    return supers[s].ReducedLength( *this );   };

     // Get the start of the given read location on its supercontig.

     int StartOnSuper( int read_loc_index ) const;

     // Delete the specified mtig from the specified super, adjusting
     // the various data structures as necessary.  Assumes that the
     // supercontig will not break if this mtig is removed.

     void DeleteMtigInSuper( int m, int s );

     // Schedule an mtig or super for annihilation

     void ClearMtig( int m ) 
     {    mtig[m].Reinitialize( );
          mtig_qual[m].Reinitialize( );    }

     void ClearSuper( int s ) { supers[s].mtig.clear( ); }

     // Return false if contig has been scheduled for annihilation; true if not.
     
     Bool Alive( int m ) const { return mtig[m].size( ) > 0; }

     // Kill supercontigs which have at most min_contigs contigs, and 
     // ( have less bases of sequence than min_len or less reads than min_reads ).

     void KillShortStandaloneMtigs( int min_len, int min_reads = 0, 
          int min_contigs = 0 );

     // Delete all supercontigs in a specified list:

     void KillSupers( vec<int> bads );

     // clean up by annihilating mtigs and renumbering (also supers)

     void CleanUp( );

     // write assembly files to the specified directory

     void Write( const String& out_dir, vec<read_location> *mappingLocs=0 );

};


#endif
