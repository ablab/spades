// Copyright (c) 2000-2003 Whitehead Institute for Biomedical Research

#include <math.h>

#include <algorithm>

#include "AnAssemblyClass.h"
#include "AnnotatedContig.h"
#include "ContigEvent.h"
#include "CoreTools.h"
#include "FindGaps.h"
#include "Quality.h"
#include "ReadLocationUtil.h"
#include "Set.h"
#include "Superb.h"

int super::Length( const assembly& A ) const
{    int answer = 0;
     for ( int i = 0; i < (int) mtig.size( ); i++ )
     {    answer += A.mtig[ mtig[i] ].size( );
          if ( i < (int) mtig.size( ) - 1 ) answer += gap[i];    }
     return answer;    }

int super::ReducedLength( const assembly& A ) const
{    int answer = 0;
     for ( int i = 0; i < (int) mtig.size( ); i++ )
          answer += A.mtig[ mtig[i] ].size( );
     return answer;    }



int assembly::StartOnSuper( int read_loc_index ) const
{
  const read_location& r = reads_orig[read_loc_index];
  ForceAssert( ! r.IsDead() );
  int mtig_id = r.Contig();
  map<int, int>::const_iterator iter;
  iter = mtigs_to_super_pos.find( mtig_id );
  ForceAssert( iter != mtigs_to_super_pos.end( ) );
  int mtig_pos = iter->second;
  iter = mtigs_to_supers.find( mtig_id );
  ForceAssert( iter != mtigs_to_supers.end( ) );
  int super_id = iter->second;
//   int mtig_pos = mtigs_to_super_pos[mtig_id];
//   int super_id = mtigs_to_supers[mtig_id];
  const super& sup = supers[super_id];

  int pos = 0;

  for ( int i = 0; i < mtig_pos; ++i )
  {
    pos += mtig[ sup.mtig[i] ].size( );
    pos += sup.gap[i];
  }

  pos += r.StartOnContig();

  return pos;
}

void assembly::DeleteMtigInSuper( int mtig_id, int super_id )
{
  ForceAssertEq( mtigs_to_supers[mtig_id], super_id );

  super& the_super = supers[super_id];

  int n = mtigs_to_super_pos[ mtig_id ];

  // If the mtig to be removed is not at either end, adjust the gap
  // for the previous contig to account for the removal of this contig.
  if ( n > 0 && n < (int) the_super.gap.size() )
    the_super.gap[n-1] = the_super.gap[n-1] + mtig[mtig_id].size() + the_super.gap[n] ;

  // Remove the gap between this mtig and the next from the gaps list.
  if ( the_super.gap.size() > 0 )
  {
    for ( int j = n + 1; j < (int) the_super.gap.size( ); j++ )
      the_super.gap[j - 1] = the_super.gap[j];
    the_super.gap.resize( the_super.gap.size( ) - 1 );
  }

  // Remove the index of this mtig from the super and adjust the mtig
  // to super pos map accordingly.
  for ( int j = n + 1; j < (int) the_super.mtig.size( ); j++ )
  {
    the_super.mtig[j - 1] = the_super.mtig[j];
    mtigs_to_super_pos[ the_super.mtig[j - 1] ] = j - 1;
  }
  the_super.mtig.resize( the_super.mtig.size( ) - 1 );

  mtigs_to_supers[mtig_id] = -2;
  mtigs_to_super_pos[mtig_id] = -2;

  // Schedule the mtig for removal.
  this->ClearMtig( mtig_id );

  // Kill the read locations associated with this mtig and adjust the
  // read to read location map accordingly.
  vec<int>& v = reads_orig_index[mtig_id];
  for ( int k = 0; k < (int) v.size( ); k++ )
  {
    read_location& r1 = reads_orig[ v[k] ];
    r1.Kill();
    simple_reads_orig_index[ r1.ReadId( ) ] = -1;
  }

  delete_contig_event e(mtig_id);
  if ( !elogp )
    cerr << "Unable to write to contig event log." << endl;
  e.BinaryWrite( *elogp );
}

void assembly::KillShortStandaloneMtigs( int min_len, int min_reads,
     int min_contigs )
{
     set<int> dead_mtigs;

     // Process deletions.
     if ( !elogp )
       cerr << "Unable to write to contig event log." << endl;

     for ( unsigned int i = 0; i < supers.size( ); i++ )
     {    super& s = supers[i];
          if ( (int) s.mtig.size( ) <= min_contigs )
          {    longlong total_bases = 0;
               int total_reads = 0;
               for ( int j = 0; j < (int) s.mtig.size( ); j++ )
               {    int m = s.mtig[j];
                    total_bases += mtig[m].size( );
                    total_reads += reads_orig_index[m].size( );    }
               Bool too_short = False;
               if ( total_bases < min_len ) too_short = True;
               if ( total_reads < min_reads ) too_short = True;
               if (too_short)
               {    for ( int j = 0; j < (int) s.mtig.size( ); j++ )
                    {    int m = s.mtig[j];
                         mtig[m].Setsize(0);
                         mtig_qual[m].resize(0);
                         dead_mtigs.insert(m);
                         delete_contig_event e(m);
                         e.BinaryWrite( *elogp );    }
                    s.mtig.clear( );    }    }    }

     // Clean up the other stuff.

     int current = 0;
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          if ( !Member( dead_mtigs, reads[i].Contig( ) ) )
               reads[current++] = reads[i];
     vec<int> dead_reads;
     reads.resize(current);
     current = 0;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
     {    if ( !Member( dead_mtigs, reads_orig[i].Contig( ) ) )
               reads_orig[current++] = reads_orig[i];
          else dead_reads.push_back( reads_orig[i].ReadId( ) );    }
     reads_orig.resize(current);
     reads_index.clear( );
     reads_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          reads_index[ reads[i].Contig( ) ].push_back(i);
     reads_orig_index.clear( );
     reads_orig_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          reads_orig_index[ reads_orig[i].Contig( ) ].push_back(i);
     for ( int i = 0; i < N; i++ )
          simple_reads_orig_index[i] = -1;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          simple_reads_orig_index[ reads_orig[i].ReadId( ) ] = i;
     current = 0;
     sort( dead_reads.begin( ), dead_reads.end( ) );
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
     {    if ( BinPosition( dead_reads, pairs[i].id1 ) < 0
               && BinPosition( dead_reads, pairs[i].id2 ) < 0 )
               pairs[current++] = pairs[i];    }
     pairs.resize(current);
     for ( int i = 0; i < N; i++ )
          pairs_index[i] = -1;
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
          if ( pairs[i].Alive( ) )
               pairs_index[ pairs[i].id1 ] = pairs_index[ pairs[i].id2 ] = i;    }

void assembly::KillSupers( vec<int> bads )
{    sort( bads.begin( ), bads.end( ) );
     vec<int> dead_mtigs;
     for ( int i = 0; i < (int) supers.size( ); i++ )
     {    super& s = supers[i];
          if ( !BinMember( bads, i ) ) continue;
          for ( int j = 0; j < (int) s.mtig.size( ); j++ )
          {    int m = s.mtig[j];
               mtig[m].Setsize(0);
               mtig_qual[m].resize(0);
               dead_mtigs.push_back(m);    }
          s.mtig.resize(0);    }
     sort( dead_mtigs.begin( ), dead_mtigs.end( ) );
     int current = 0;
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          if ( !BinMember( dead_mtigs, reads[i].Contig( ) ) )
               reads[current++] = reads[i];
     vec<int> dead_reads;
     reads.resize(current);
     current = 0;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
     {    if ( !BinMember( dead_mtigs, reads_orig[i].Contig( ) ) )
               reads_orig[current++] = reads_orig[i];
          else dead_reads.push_back( reads_orig[i].ReadId( ) );    }
     reads_orig.resize(current);
     reads_index.clear( );
     reads_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          reads_index[ reads[i].Contig( ) ].push_back(i);
     reads_orig_index.clear( );
     reads_orig_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          reads_orig_index[ reads_orig[i].Contig( ) ].push_back(i);
     for ( int i = 0; i < N; i++ )
          simple_reads_orig_index[i] = -1;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          simple_reads_orig_index[ reads_orig[i].ReadId( ) ] = i;
     current = 0;
     sort( dead_reads.begin( ), dead_reads.end( ) );
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
     {    if ( BinPosition( dead_reads, pairs[i].id1 ) < 0
               && BinPosition( dead_reads, pairs[i].id2 ) < 0 )
               pairs[current++] = pairs[i];    }
     pairs.resize(current);
     for ( int i = 0; i < N; i++ )
          pairs_index[i] = -1;
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
          if ( pairs[i].Alive( ) )
               pairs_index[ pairs[i].id1 ] = pairs_index[ pairs[i].id2 ] = i;    }

// clean up by annihilating mtigs and renumbering (also supers)

void assembly::CleanUp( )
{    map<size_t, size_t> new_mtig;
     size_t current = 0;
     if ( !elogp )
       cerr << "Unable to write to contig event log." << endl;
     for ( size_t i = 0; i < mtig.size( ); i++ )
          if ( mtig[i].size( ) != 0 )
          {
               // In effect, do:
               //      mtig[current] = mtig[i];
               //      mtig_qual[current] = mtig_qual[i];
               // but do it more efficiently, memory-wise.

               if ( current != i )
               {    mtig[current].Reinitialize( );
                    mtig[current].Swap( mtig[i] );
                    mtig_qual[current].Reinitialize( );
                    mtig_qual[current].Swap( mtig_qual[i] );    }

               rename_contig_event e( i, current );
               e.BinaryWrite( *elogp );
               reads_index[current] = reads_index[i];
               reads_orig_index[current] = reads_orig_index[i];
               mtig_aligns_lost[current] = mtig_aligns_lost[i];
               new_mtig[i] = current;
               current++;    }
     mtig.resize(current);
     mtig_qual.resize(current);
     reads_index.resize(current);
     reads_orig_index.resize(current);

     for ( size_t i = 0; i < mtig_aligns.size( ); i++ )
     {    mtig_aligns[i].id1 = new_mtig[ mtig_aligns[i].id1 ];
          mtig_aligns[i].id2 = new_mtig[ mtig_aligns[i].id2 ];    }
     for ( size_t i = 0; i < mtig.size( ); i++ )
          mtig_aligns_index[i] = -1;
     size_t i = mtig_aligns.size();
     while ( i-- )
        mtig_aligns_index[ mtig_aligns[i].id1 ] = i; 
       
     mtig_aligns_lost.resize(current);

     for ( size_t i = 0; i < reads.size( ); i++ )
     {    if ( reads[i].Contig( ) >= 0 )
               reads[i].SetContig( new_mtig[ reads[i].Contig( ) ] );    }
     for ( size_t i = 0; i < reads_orig.size( ); i++ )
     {    if ( reads_orig[i].Contig( ) >= 0 )
               reads_orig[i].SetContig( new_mtig[ reads_orig[i].Contig( ) ] );    }

     for ( size_t i = 0; i < supers.size( ); i++ )
          for ( size_t j = 0; j < supers[i].mtig.size( ); j++ )
               supers[i].mtig[j] = new_mtig[ supers[i].mtig[j] ];
     current = 0;
     for ( size_t i = 0; i < supers.size( ); i++ )
          if ( supers[i].mtig.size( ) > 0 ) supers[current++] = supers[i];
     supers.resize(current);

     mtigs_to_supers.clear( );
     mtigs_to_super_pos.clear( );
     for ( size_t i = 0; i < supers.size( ); i++ )
          for ( size_t j = 0; j < supers[i].mtig.size( ); j++ )
          {    mtigs_to_supers[ supers[i].mtig[j] ] = i;
               mtigs_to_super_pos[ supers[i].mtig[j] ] = j;    }    }

assembly::assembly( String fasta_file,
		    String qual_file,
		    const vec<read_location>& reads_arg,
		    const vec<read_location>& reads_orig_arg,
		    const vec<super>& supers_arg,
		    const vec<brief_align>& all_aligns_arg,
		    int N_arg,
		    const vec<read_pairing>& pairs_arg,
		    const vec<int>& orig_read_lengths_arg,
		    const vec<int>& all_aligns_index_arg,
		    Bool store_contigs_one_by_one )
          : reads(reads_arg),
	    reads_orig(reads_orig_arg),
	    supers(supers_arg),
	    all_aligns(all_aligns_arg),
	    N(N_arg),
	    pairs(pairs_arg),
	    orig_read_lengths(orig_read_lengths_arg),
	    all_aligns_index(all_aligns_index_arg),
	    elogp(0)
{

     // Read in contigs.

     if ( !store_contigs_one_by_one )
     {    mtig.ReadAll(fasta_file);
          if ( qual_file != "" ) mtig_qual.ReadAll(qual_file);
          else mtig_qual.resize( mtig.size( ) );    }
     else
     {    {    vecqualvector mtig_qual_temp;
               mtig_qual_temp.ReadAll(qual_file);
               mtig_qual.reserve( mtig_qual_temp.size( ) );
               for ( vecqvec::size_type i = 0; i < mtig_qual_temp.size( ); i++ )
                    mtig_qual.push_back_external( mtig_qual_temp[i] );    }
          {    vecbasevector mtig_temp;
               mtig_temp.ReadAll(fasta_file);
               mtig.reserve( mtig_temp.size( ) );
               for ( size_t i = 0; i < mtig_temp.size( ); i++ )
                    mtig.push_back_external( mtig_temp[i] );    }    }

     // Do everything else.

     reads_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads.size( ); i++ )
          reads_index[ reads[i].Contig( ) ].push_back(i);
     reads_orig_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          reads_orig_index[ reads_orig[i].Contig( ) ].push_back(i);
     pairs_index.resize(N);
     for ( int i = 0; i < N; i++ )
          pairs_index[i] = -1;
     for ( unsigned int i = 0; i < pairs.size( ); i++ )
          if ( pairs[i].Alive( ) )
               pairs_index[ pairs[i].id1 ] = pairs_index[ pairs[i].id2 ] = i;
     for ( unsigned int i = 0; i < supers.size( ); i++ )
          for ( unsigned int j = 0; j < supers[i].mtig.size( ); j++ )
          {    mtigs_to_supers[ supers[i].mtig[j] ] = i;
               mtigs_to_super_pos[ supers[i].mtig[j] ] = j;    }
     simple_reads_orig_index.resize(N);
     for ( int i = 0; i < N; i++ )
          simple_reads_orig_index[i] = -1;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          simple_reads_orig_index[ reads_orig[i].ReadId( ) ] = i;
     SetMtigAligns(mtig_aligns);    }

void assembly::Write( const String& out_dir, vec<read_location> *mapping_locs )
{
     // Eliminate dead reads.

     EraseIf( reads, &read_location::IsDead );
     EraseIf( reads_orig, &read_location::IsDead );

     // Recompute reads_orig_index and simple_reads_orig_index.

     reads_orig_index.clear( ), reads_orig_index.resize( mtig.size( ) );
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          reads_orig_index[ reads_orig[i].Contig( ) ].push_back(i);

     simple_reads_orig_index.clear(), simple_reads_orig_index.resize(N);
     for ( int i = 0; i < N; i++ )
          simple_reads_orig_index[i] = -1;
     for ( unsigned int i = 0; i < reads_orig.size( ); i++ )
          simple_reads_orig_index[ reads_orig[i].ReadId( ) ] = i;

     vec<superb> bs( supers.size( ) );
     for ( unsigned int i = 0; i < supers.size(); ++i )
     {    static vec<int> old_gaps, gaps_dev;
          old_gaps = supers[i].gap;
          FindGaps( *this, simple_reads_orig_index, i, gaps_dev, 0,
               supers[i].mtig.size( ) - 1 );
          super& s = supers[i];
          int n = s.mtig.size( );
          bs[i].SetNtigs(n);
          for ( int j = 0; j < n; j++ )
          {    int m = s.mtig[j];
               bs[i].SetTig( j, m );
               bs[i].SetLen( j, mtig[m].size( ) );
               if ( j < n - 1 )
               {    bs[i].SetGap( j, old_gaps[j] );
	            bs[i].SetDev( j, gaps_dev[j] );    }    }    }
     WriteSupercontigFiles( out_dir, bs );

  this->mtig.WriteAll( out_dir + "/mergedcontigs.fastb" );

  WriteQualityScores( out_dir + "/mergedcontigs.qualb", this->mtig_qual, False );

  Sort( this->reads );
  if ( this->reads.empty() )
    // No indexing.
    WriteLocs( out_dir + "/mergedcontigs.locs", this->reads );
  else
    // Indexed.
    WriteLocs( out_dir + "/mergedcontigs.locs", this->reads, this->mtig.size(), this->N );

  if ( mapping_locs != 0 )
  {
    vec<read_location> orig_locs;
    orig_locs.reserve(reads_orig.size()+mapping_locs->size());
    UnmergeLocations( reads_orig, *mapping_locs, orig_locs );
    Sort(orig_locs);

    reads_orig.swap(orig_locs);
  }

  Sort( this->reads_orig );
  if ( this->reads_orig.empty() )
    // No indexing.
    WriteLocs( out_dir + "/mergedcontigs_orig.locs", this->reads_orig );
  else
    // Indexed.
    WriteLocs( out_dir + "/mergedcontigs_orig.locs", this->reads_orig, this->mtig.size(), this->N );

}
