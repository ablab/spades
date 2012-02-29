///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Equiv.h"
#include "Intvector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "TaskTimer.h"
#include "Vec.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/OffsetDistributionInterface.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/reporting/CBundle.h"
#include "paths/reporting/CSuperLinks.h"

// # define FOLLOW_CONTIGS
#ifdef FOLLOW_CONTIGS   // for debugging only
const int CONTIG1 = 115;
const int CONTIG2 = 109;
#endif

/**
 * RegapSupersCore
 */
void RegapSupersCore( ostream &log,
		      vec<superb> &supers,
		      const int super_id,
		      const vec< pair<int,CBundle> > &bundles,
		      const bool VERBOSE,
		      const int MAX_OVERLAP,
		      const bool FIX_NEG_GAPS )
{
  // Nothing to do.
  const superb &sup = supers[super_id];
  const int ntigs = sup.Ntigs( );
  if ( ntigs < 2 ) return;
  
  // The super being regapped.
  superb regapped = sup;

  // Maps to original supers.
  vec<int> start( ntigs, 0 );
  int runner = 0;
  for (int ii=0; ii<ntigs; ii++) {
    start[ii] = runner;
    runner += sup.Len( ii );
    if ( ii < ntigs - 1 ) runner += sup.Gap( ii );
  }
  
  // Log start.
  log << Date( ) << ": regapping super_" << super_id << endl;

  // Find interval in bundles of links for the contigs of super_id.
  vec< pair<int,CBundle> >::const_iterator it;
  
  CBundle emptyb;
  pair<int,CBundle> sniffer = make_pair( super_id, emptyb );
  it = lower_bound( bundles.begin( ), bundles.end( ), sniffer);
  if ( it == bundles.end( ) || it->first != super_id ) {
    log << "ERROR! No links found for super_" << super_id << endl;
    return;
  }
  int begin_bid = distance( bundles.begin( ), it );
  int end_bid = -1;
  for (int ii=begin_bid+1; ii<bundles.isize( ); ii++) {
    if ( bundles[ii].first != super_id ) {
      end_bid = ii;
      break;
    }
  }
  if ( end_bid < 0 ) end_bid = bundles.isize( );
  
  // Sort links by weight.
  vec< pair<int,int> > weight2bid;
  weight2bid.reserve( end_bid - begin_bid );
  for (int ii=begin_bid; ii<end_bid; ii++)
    weight2bid.push_back( make_pair( bundles[ii].second.Weight( ), ii ) );
  sort( weight2bid.rbegin( ), weight2bid.rend( ) );
  
  // Components: connected sets of (cgpos, start)
  vec< vec< pair<int,int> > > comp;
  comp.reserve( ntigs );
  for (int ii=0; ii<ntigs; ii++) {
    vec< pair<int,int> > singlet = MkVec( make_pair( ii, 0 ) );
    comp.push_back( singlet );
  }

  // Place contigs (start with heaviest bundles).
  for (int ii=0; ii<weight2bid.isize( ); ii++) {
    const int bid = weight2bid[ii].second;
    const CBundle &bundle = bundles[bid].second;
    int offset = bundle.Offset( ).first;
    int p1 = bundle.Id1( );
    int p2 = bundle.Id2( );

     // Find components containing the two ids.
     int comp_id1 = -1;
     int comp_id2 = -1;
     int base1 = 0;
     int base2 = 0;
     for (int jj=0; jj<comp.isize( ); jj++) {
       if ( comp_id1 > -1 && comp_id2 > -1 ) break;
       for (int kk=0; kk<comp[jj].isize( ); kk++) {
	 if ( comp_id1 > -1 && comp_id2 > -1 ) break;
	 if ( comp_id1 < 0 && comp[jj][kk].first == p1 ) {
	   comp_id1 = jj;
	   base1 = comp[jj][kk].second;
	 }
	 if ( comp_id2 < 0 && comp[jj][kk].first == p2 ) {
	   comp_id2 = jj;
	   base2 = comp[jj][kk].second;
	 }
       }
     }
     ForceAssertGe( comp_id1, 0 );
     ForceAssertGe( comp_id2, 0 );

     // Already joined.
     if ( comp_id1 == comp_id2 ) continue;

     // Merge components (use offset).
     comp[comp_id1].reserve( comp[comp_id1].size( ) + comp[comp_id2].size( ) );
     for (int jj=0; jj<comp[comp_id2].isize( ); jj++) {
       int old_start = comp[comp_id2][jj].second;
       int new_start = base1 + offset - base2 + old_start;
       comp[comp_id2][jj].second = new_start;
       comp[comp_id1].push_back( comp[comp_id2][jj] );
     }
     comp[comp_id2].clear( );

   }

   // Should have one non empty component at the end.
   vec<int> comp_ids;
   for (int ii=0; ii<comp.isize( ); ii++)
     if ( comp[ii].size( ) > 0 )
       comp_ids.push_back( ii );
   if ( comp_ids.size( ) != 1 ) {
     log << "ERROR! There are " << comp_ids.size( ) << " components.\n";
     for (int ii=0; ii<comp_ids.isize( ); ii++) {
       int cid = comp_ids[ii];
       vec< pair<int,int> > start2pos = comp[cid];
       for (int jj=0; jj<start2pos.isize( ); jj++)
	 swap( start2pos[jj].first, start2pos[jj].second );
       sort( start2pos.begin( ), start2pos.end( ) );
       for (int jj=0; jj<start2pos.isize( ); jj++)
	 log << start2pos[jj].second << " (c"
	     << sup.Tig( start2pos[jj].second ) << ") ";
       log << endl;
     }

     // Skip super (maybe this should be an assert?)
     return;
   }

   // Sort contigs.
   int final_id = comp_ids[0];
   vec< pair<int,int> > fincomp;   // (start, pos_orig)
   fincomp.reserve( comp[final_id].size( ) );
   for (int ii=0; ii<comp[final_id].isize( ); ii++) {
     int start = comp[final_id][ii].second;
     int cpos = comp[final_id][ii].first;
     fincomp.push_back( make_pair( start, cpos ) );
   }
   ForceAssertEq( fincomp.isize( ), ntigs );
   sort( fincomp.begin( ), fincomp.end( ) );

   // Initialize regapped.
   regapped.SetNtigs( ntigs ); 
   for (int ii=0; ii<ntigs; ii++) {
     int pos_orig = fincomp[ii].second;
     int tig_len = sup.Len( pos_orig );

     regapped.SetTig( ii, sup.Tig( pos_orig ) );
     regapped.SetLen( ii, tig_len );
     if ( ii < ntigs - 1 ) {
       int gap_size = fincomp[ii+1].first - fincomp[ii].first - tig_len;
       regapped.SetGap( ii, gap_size );
       regapped.SetDev( ii, Max( Abs( gap_size ) / 4, 1 ) );
     }
   }

   // Get rid of large negative gaps (swap positions of contigs).
   if ( FIX_NEG_GAPS ) {
     for (int cpos=0; cpos<regapped.Ntigs( )-1; cpos++) {
       int start_cpos = 0;
       for (int ii=0; ii<cpos; ii++)
	 start_cpos += regapped.Len( ii ) + regapped.Gap( ii );
       int end_cpos = start_cpos + regapped.Len( cpos );
       int start_next = end_cpos + regapped.Gap( cpos );
       int end_next = start_next + regapped.Len( cpos + 1 );
       int start_next_next
	 = ( cpos + 2 < regapped.Ntigs( ) )
	 ? end_next + regapped.Gap( cpos + 1 )
	 : 0;
       
       // Positive gap, nothing to do.
       int current_gap = start_next - end_cpos;
       if ( current_gap >= 0 ) continue;
       
       // Negative gap, but no need to swap.
       int swapped_gap = start_cpos - end_next;
       if ( current_gap >= swapped_gap ) continue;
       
       // Swap cpos with next (contigs).
       int tig1 = regapped.Tig( cpos + 1 );
       int len1 = regapped.Len( cpos + 1 );
       int tig2 = regapped.Tig( cpos );
       int len2 = regapped.Len( cpos );
       regapped.SetTig( cpos, tig1 );
       regapped.SetLen( cpos, len1 );
       regapped.SetTig( cpos + 1, tig2 );
       regapped.SetLen( cpos + 1, len2 );
       
       // Swap cpos with next (gaps).
       if ( cpos >= 1 ) {
	 int new_gap = regapped.Gap( cpos - 1 ) + ( start_next - start_cpos );
	 int new_dev = Max( Abs( new_gap ) / 4, 1 );
	 regapped.SetGap( cpos - 1, new_gap );
	 regapped.SetDev( cpos - 1, new_dev );
       }
       regapped.SetGap( cpos, swapped_gap );
       regapped.SetDev( cpos, Max( Abs( swapped_gap ) / 4, 1 ) );
       if ( cpos + 2 < regapped.Ntigs( ) ) {
	 int new_gap = start_next_next - end_cpos;
	 int new_dev = Max( Abs( new_gap ) / 4, 1 );
	 regapped.SetGap( cpos + 1, new_gap );
	 regapped.SetDev( cpos + 1, new_dev );
       }
     }
   }     

   // A map orig_pos (in sup) to new_pos (in regapped).
   vec<int> orig2new( ntigs, -1 );
   for (int old_pos=0; old_pos<ntigs; old_pos++) {
     int contig_id = sup.Tig( old_pos );
     for (int jj=0; jj<ntigs; jj++) {
       if ( regapped.Tig( jj ) == contig_id ) {
	 orig2new[old_pos] = jj;
	 break;
       }
     }
     ForceAssertNe( orig2new[old_pos], -1 );
   }

   // Generate a local version of bundles (ids are new pos in regapped).
   vec<CBundle> rebundles;
   rebundles.reserve( end_bid - begin_bid );
   for(int bundle_id=begin_bid; bundle_id<end_bid; bundle_id++) {
     CBundle newb = bundles[bundle_id].second;
     int old1 = newb.Id1( );
     int old2 = newb.Id2( );
     newb.SetId1( orig2new[old1] );
     newb.SetId2( orig2new[old2] );
     rebundles.push_back( newb );
   }

   // For each gap, find all the bundles above it.
   vec< vec<int> > gap2bundles( ntigs - 1 );
   for (int ii=0; ii<rebundles.isize( ); ii++) {
     int id1 = rebundles[ii].Id1( );
     int id2 = rebundles[ii].Id2( );
     for (int jj=id1; jj<id2; jj++)
       gap2bundles[jj].push_back( ii );
   }

   // Reestimate gaps, one at a time (twice).
   for (int loop=0; loop<2; loop++) {
     for (int gap_id=0; gap_id<ntigs-1; gap_id++) {
       vec< pair<int,normal_distribution> > w2nds;
       for (int ii=0; ii<gap2bundles[gap_id].isize( ); ii++) {
	 const CBundle &theb = rebundles[ gap2bundles[gap_id][ii] ];
	   
	 int id1 = theb.Id1( );
	 int id2 = theb.Id2( );
	 int weight = theb.Weight( );
	 pair<int,int> offset_base = theb.Offset( );
	 int themu = offset_base.first;
	 for (int jj=id1; jj<id2; jj++) {
	   themu += -regapped.Len( jj );
	   if ( jj != gap_id ) themu += -regapped.Gap( jj );
	 }
	   
	 normal_distribution theND( themu, offset_base.second );
	 w2nds.push_back( make_pair( weight, theND ) );
       }
	 
       sort( w2nds.rbegin( ), w2nds.rend( ) );
	 
       // Want only heaviest weight.
       const double MAX_RATIO = 0.15;   // HEURISTICS
       int total_weight = 0;
       vec<normal_distribution> nds;
       for (int ii=0; ii<w2nds.isize( ); ii++) {
	 double weight = w2nds[ii].first;
	 if ( w2nds[ii].second.mu_ < - MAX_OVERLAP ) break;  // not continue
	 if ( weight < MAX_RATIO * double( total_weight ) ) break;
	 nds.push_back( w2nds[ii].second );
	 total_weight += weight;

#ifdef FOLLOW_CONTIGS
	 int cg1 = regapped.Tig( gap_id );
	 int cg2 = regapped.Tig( gap_id + 1 );
	 if ( cg1 == CONTIG1 && cg2 == CONTIG2 ) {
	   int gap_size = w2nds[ii].second.mu_;
	   int gap_dev = w2nds[ii].second.sigma_;
	   
	   cout << "Fixing gap between c" << cg1 << " and c" << cg2 << " "
		<< " before: " << gap_size
		<< " +/- " << gap_dev
		<< " (" << w2nds[ii].first << " links)"
		<< "\n";
	 }
#endif
       }
       
       // Accept "overlapping" bundles if there is nothing else to connect.
       if ( nds.size( ) < 1 ) {
	 if ( VERBOSE ) 
	   log << "WARNING! There are no valid links above gap "
	       << "c" << regapped.Tig( gap_id ) << " - " 
	       << "c" << regapped.Tig( gap_id + 1 ) << "\n"
	       << "Recycling the following bundles:\n";
	 for (int ii=0; ii<w2nds.isize( ); ii++) {
	   double weight = w2nds[ii].first;
	   if ( weight < MAX_RATIO * double( total_weight ) ) break;
	   if ( VERBOSE )
	     log << w2nds[ii].second.mu_ << " +/- "
		 << w2nds[ii].second.sigma_ << "\tw: "
		 << w2nds[ii].first << "\n";
	   nds.push_back( w2nds[ii].second );
	   total_weight += weight;
	 }
       }
	 
       normal_distribution fixed = CombineNormalDistributions( nds );
       regapped.SetGap( gap_id, (int)fixed.mu_ );
       regapped.SetDev( gap_id, Max( 1, (int)fixed.sigma_ ) );
	 
#ifdef FOLLOW_CONTIGS
       int cg1 = regapped.Tig( gap_id );
       int cg2 = regapped.Tig( gap_id + 1 );
       if ( cg1 == CONTIG1 && cg2 == CONTIG2 ) {
	 cout << "Fixing gap between c" << cg1 << " and c" << cg2 << " "
	      << "after: " << (int)fixed.mu_
	      << " +/- " << Max( 1, (int)fixed.sigma_ )
	      << "\n";
       }
#endif

     }
   }
   
   // Done.
   supers[super_id] = regapped;

 }

 /**
  * BreakUnlinked
  */
 int BreakUnlinked( ostream &log,
		    vec<superb> &supers,
		    const vec< pair<int,CBundle> > &bundles,
		    const bool VERBOSE )
 {
   // The extras from unconnected components.
  vec<superb> extras;

  // Loop over all supers.
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    if ( sup.Ntigs( ) < 2 ) continue;
    
    // Find clusters of linked contigs.
    vec< vec<int> > clusters;
    CBundle emptyb;
    pair<int,CBundle> sniffer = make_pair( super_id, emptyb );
    vec< pair<int,CBundle> >::const_iterator it
      = lower_bound( bundles.begin( ), bundles.end( ), sniffer );
    if ( it == bundles.end( ) || it->first != super_id ) {
      clusters.reserve( sup.Ntigs( ) );
      for (int ii=0; ii<sup.Ntigs( ); ii++) {
	vec<int> singleton = MkVec( ii );
	clusters.push_back( singleton );
      }
    }
    else {
      equiv_rel rel( sup.Ntigs( ) );
      int begin_bid = distance( bundles.begin( ), it );
      int end_bid = -1;
      for (int ii=begin_bid+1; ii<bundles.isize( ); ii++) {
	if ( bundles[ii].first != super_id ) {
	  end_bid = ii;
	  break;
	}
      }
      if ( end_bid < 0 ) end_bid = bundles.isize( );
      
      for (int bid=begin_bid; bid<end_bid; bid++) {
	ForceAssertEq( bundles[bid].first, super_id );
	const CBundle &bundle = bundles[bid].second;
	int p1 = bundle.Id1( );
	int p2 = bundle.Id2( );
	rel.Join( p1, p2 );
      }
      clusters.reserve( rel.OrbitCount( ) );
      vec<int> reps;
      rel.OrbitReps( reps );
      for (int ii=0; ii<rel.OrbitCount( ); ii++) {
	vec<int> orbit;
	rel.Orbit( reps[ii], orbit );
	sort( orbit.begin( ), orbit.end( ) );
	clusters.push_back( orbit );
      }
    }
    
    // Find segments.
    vec<int> tags;
    for (int clid=0; clid<clusters.isize( ); clid++) {
      tags.push_back( clusters[clid][0] );
      for (int ii=1; ii<clusters[clid].isize( ); ii++) {
	if ( clusters[clid][ii] != clusters[clid][ii-1] + 1 )
	  tags.push_back( clusters[clid][ii] );
      }
    }
    sort( tags.begin( ), tags.end( ) );
    
    // Break super.
    if ( tags.size( ) < 2 ) continue;
    
    vec<superb> chunks;
    log << "breaking s" << super_id << " into: ";
    for (int ii=0; ii<tags.isize( ); ii++) {
      vec<int> select;
      int beg = tags[ii];
      int end = ii < tags.isize( ) - 1 ? tags[ii+1] : sup.Ntigs( );
      log << "[" << beg << ", " << end << ") ";
      for (int jj=beg; jj<end; jj++)
	select.push_back( jj );
      chunks.push_back( sup.SubSuper( select ) );
    }
    log << "\n";
    
    supers[super_id] = chunks[0];
    for (int ii=1; ii<chunks.isize( ); ii++)
      extras.push_back( chunks[ii] );

  } // loop over all supers.
  
  // Add extra chunks.
  supers.reserve( supers.size( ) + extras.size( ) );
  for (int ii=0; ii<extras.isize( ); ii++)
    supers.push_back( extras[ii] );
  
  // Return.
  return extras.isize( );

}

/**
 * FindInSuperBundles
 */
void FindInSuperBundles( ostream &log,
			 vec< pair<int,CBundle> > &bundles,
			 const PairsManager &pairs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 const vec<superb> &supers,
			 const bool VERBOSE,
			 const int MIN_LINKS,
			 const int MAX_DISCREPANCY,
			 const UInt64VecVec *to_pairs,
			 const vec<IntDistribution> *distr )
{
  bundles.clear( );
  
  // Maps for the original supers.
  int n_contigs = -1;
  for (int ii=0; ii<supers.isize( ); ii++)
    for (int cpos=0; cpos<supers[ii].Ntigs( ); cpos++)
      n_contigs = Max( n_contigs, supers[ii].Tig( cpos ) );
  n_contigs++;

  vec<int> to_super( n_contigs, -1 );
  vec<int> to_pos( n_contigs, -1 );
  vec<int> cg_len( n_contigs, 0 );
  vec<int> start( n_contigs, 0 );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int pos = 0;
    for (int cpos=0; cpos<supers[ii].Ntigs( ); cpos++) {
      to_super[ supers[ii].Tig( cpos ) ] = ii;
      to_pos[ supers[ii].Tig( cpos ) ] = cpos;
      cg_len[ supers[ii].Tig( cpos ) ] = supers[ii].Len( cpos );
      start[ supers[ii].Tig( cpos ) ] = pos;
      pos += supers[ii].Len( cpos );
      if ( cpos < supers[ii].Ntigs( ) - 1 )
	pos += supers[ii].Gap( cpos );
    }
  }
  
  // Build singleton_supers (new super_id == old contig_id).
  vec<superb> singleton_supers;
  singleton_supers.resize( n_contigs );
  for (int tig=0; tig<n_contigs; tig++) {
    singleton_supers[tig].SetNtigs( 1 );
    singleton_supers[tig].SetTig( 0, tig );
    singleton_supers[tig].SetLen( 0, cg_len[tig] );
  }
  
  // The linker.
  CSuperLinks linker( &pairs, &singleton_supers, &aligns, &index );
  
  // Loop over all singleton_supers.
  for (int contig_id=0; contig_id<singleton_supers.isize( ); contig_id++) {
    float slop = -1.0;
    float stretch = 1.0;
    vec<COffset> all_offsets;
    linker.AllLinks( contig_id, all_offsets, slop, &stretch );
    
    // Loop over all other contigs linked to contig_id.
    for (int set_id=0; set_id<all_offsets.isize( ); set_id++) {
      const COffset &offset = all_offsets[set_id];
      if ( VERBOSE ) offset.Print( log );
      
      // Some contigs not in supers.
      int id1 = offset.Super1( );
      int id2 = offset.Super2( );
      bool rc2 = offset.Rc2( );
      if ( to_super[id1] < 0 || to_super[id2] < 0 ) continue;
      
      // Not fw-fw, or not in the same super.
      if ( rc2 || ( to_super[id1] != to_super[id2] ) ) continue;
      int pos1 = to_pos[id1];
      int pos2 = to_pos[id2];

      // All candidate bundles between id1[+] and id2[+] (scored).
      vec<int> weights;
      vec<CBundle> candidates; 
      size_t n_clusters = offset.NClusters( );
      for (size_t cluster_id=0; cluster_id<n_clusters; cluster_id++) {
	int weight = offset.NLinks( cluster_id );
	if ( weight < MIN_LINKS ) continue;
	float score = offset.Score( cluster_id );
	normal_distribution nd = offset.Offset( cluster_id );
	int offm = int( nd.mu_ );
	int offe = offset.DevErrorEstimate( cluster_id, pairs );
	int offs = Max( 1, int( nd.sigma_ ) );
	pair<int,int> off_pair = make_pair( offm, offs + offe );
	pair<int,int> h1 = offset.SpreadWinBases1( cluster_id );
	pair<int,int> h2 = offset.SpreadWinBases2( cluster_id );
	CBundle bundle( pos1, pos2, rc2, weight, score, off_pair, h1, h2 );
	
	// Filter implausible bundles.
	int boffset = bundle.Offset( ).first;
	int start1 = start[id1];
	int start2 = start[id2];
	if ( Abs( start2 - start1 - boffset ) > MAX_DISCREPANCY ) continue;
	
	// Add bundle.
	weights.push_back( weight );
	candidates.push_back( bundle );
      }

      // No valid links found.
      if ( candidates.size( ) < 1 ) continue;
      
      // Pick the winner bundle.
      ReverseSortSync( weights, candidates );
      bundles.push_back( make_pair( to_super[id1], candidates[0] ) );
    } // loop over all contigs linked to contig_id

  } // loop over all singleton_supers
  
  // Optionally, use OffsetDistribution to evaluate the offset.
  if ( to_pairs ) {
    TaskTimer timer;
    timer.Start( );
    size_t dotter = 10;
    log << Date( ) << ": starting OffsetDistribution" << endl;
    log << "Regapping " << bundles.size( )
	<< " gaps (.=" << dotter
	<< " gaps):" << endl;
    for (size_t ii=0; ii<bundles.size( ); ii++) {
      if ( ii % dotter == 0 ) Dot( log, ii / dotter );
      const superb &sup = supers[ bundles[ii].first ];

      // SANTEMP
      const CBundle &theB = bundles[ii].second;
      int super_id = bundles[ii].first;
      int cg1 = sup.Tig( theB.Id1( ) );
      int cg2 = sup.Tig( theB.Id2( ) );

      // if ( cg1 != 102 || cg2 != 4 ) continue;

      cout << "gap " << ii << "." << bundles.size( )
	   << ", between c" << cg1 << " and c" << cg2
	   << " (" << theB.Weight( ) << " links)"
	   << "   in: [" << theB.Offset( ).first - cg_len[cg1]
	   << ", " << theB.Offset( ).second
	   << ")" << flush;
      // SANTEMP

      OffsetFromDistributions( bundles[ii].second, sup, pairs, *to_pairs,
			       aligns, index, cg_len, *distr );
      
      // SANTEMP
      const CBundle &theB2 = bundles[ii].second;
      cout << "   out: [" << theB2.Offset( ).first - cg_len[cg1]
	   << ", " << theB2.Offset( ).second
	   << ")" << endl;
      // SANTEMP
      

    }
    timer.Stop( );
    log << "\n" << ToString( timer.UserSecs( ), 2 )
	<< "/" << ToString( timer.UserSecs( ), 2 )      
	<< " seconds (user/system time) spent in OffsetDistribution" << endl;
    log << Date( ) << ": done OffsetDistribution" << endl;
  }
  
  // Sort bundles.
  sort( bundles.begin( ), bundles.end( ) );
  
}

/**
 * RegapSupers
 */
void RegapSupers( ostream &log,
		  vec<superb> &supers,
		  const PairsManager &pairs,
		  const vec<alignlet> &aligns,
		  const vec<int> &index,
		  const bool VERBOSE,
		  const int MAX_OVERLAP,
		  const int MIN_LINKS ,
		  const int MAX_DISCREPANCY,
		  const bool FIX_NEG_GAPS,
		  const vec<IntDistribution> *lib_dist )
{
  log << Date( ) << ": RegapSupers - START" << endl;
  
  // A contig to pairs map.
  UInt64VecVec contig_to_pairs;
  if ( lib_dist ) {
    log << Date( ) << ": mapping contigs to pairs" << endl;
    MapContigsToPairs( contig_to_pairs, pairs, supers, aligns, index );
  }
  UInt64VecVec *to_pairs = lib_dist ? &contig_to_pairs : 0;

  // Find all bundles between all paired contigs.
  log << Date( ) << ": finding bundles" << endl;
  vec< pair<int,CBundle> > bundles;   // (super_id, bundle)
  FindInSuperBundles( log, bundles, pairs, aligns, index, supers,
		      VERBOSE, MIN_LINKS, MAX_DISCREPANCY,
		      to_pairs, lib_dist );

  // Break unlinked supers, and regenerate bundles (if needed).
  log << Date( ) << ": breaking unlinked... " << flush;
  int n_extras = BreakUnlinked( log, supers, bundles, VERBOSE );
  if ( n_extras > 0 ) {
    log << n_extras << " breaks performed" << endl;

    if ( lib_dist ) {
      log << Date( ) << ": remapping contigs to pairs" << endl;
      MapContigsToPairs( contig_to_pairs, pairs, supers, aligns, index );
    }
    UInt64VecVec *to_pairs2 = lib_dist ? &contig_to_pairs : 0;
        
    log << Date( ) << ": regenerating bundles" << endl;
    FindInSuperBundles( log, bundles, pairs, aligns, index, supers,
			VERBOSE, MIN_LINKS, MAX_DISCREPANCY,
			to_pairs2, lib_dist );
  }
  else log << " no breaks done" << endl;

  // Regap supers.
  log << Date( ) << ": regapping " << supers.size( ) << " supers" << endl;
  for (int super_id=0; super_id<supers.isize( ); super_id++)
    RegapSupersCore( log, supers, super_id, bundles,
		     VERBOSE, MAX_OVERLAP, FIX_NEG_GAPS );
  
  // Done.
  log << Date( ) << ": RegapSupers - done\n" << endl;
}

