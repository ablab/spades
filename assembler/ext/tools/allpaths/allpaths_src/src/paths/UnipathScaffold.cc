///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Basevector.h"
#include "CoreTools.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "paths/UnipathScaffold.h"
#include "graph/DigraphTemplate.h"
#include "graph/FindCells.h"

void UnipathScaffold(

     // input:

     const digraphE<linklet>& G,
     const vecbasevector& unibases,
     const int K,
     const vec<int>& to_rc,
     const vec<Bool>& exclude,
     const int min_kmers,

     // logging and debugging:

     const int FORCE_FIRST,
     const Bool SCAFFOLDING_VERBOSE,

     // output:

     vec<superb>& uscaffolds,
     vec<Bool>& circled ) 
{
     // Set up for building unipath scaffolds.  
     // Allow for the case of circular scaffolds.

     int nuni = unibases.size( );
     vec<Bool> used = exclude;
     vec< pair<int,int> > len_uni;
     for ( int u = 0; u < nuni; u++ )
     {    int nkmers = unibases[u].isize( ) - K + 1;
          if ( !exclude[u] && nkmers >= min_kmers ) len_uni.push( nkmers, u );    }
     ReverseSort(len_uni);

     // Process FORCE_FIRST.

     if ( FORCE_FIRST >= 0 )
     {    int jx;
          for ( jx = 0; jx < len_uni.isize( ); jx++ )
               if ( len_uni[jx].second == FORCE_FIRST ) break;
          ForceAssertLt( jx, len_uni.isize( ) );
          swap( len_uni[0], len_uni[jx] );    }

     // Build unipath scaffolds.

     if (SCAFFOLDING_VERBOSE) cout << "\nINITIAL SCAFFOLDS\n";
     for ( int ui = 0; ui < len_uni.isize( ); ui++ )
     {    int u = len_uni[ui].second;
          if ( used[u] ) continue;

          // Back up as far as possible.

          vec<int> seen;
          seen.push_back(u);
          int z = u;
          while(1)
          {    if ( G.To(z).empty( ) ) break;
               vec< pair<int,int> > links_uni;
               for ( int j = 0; j < G.To(z).isize( ); j++ )
               {    int z2 = G.To(z)[j];
                    if ( exclude[z2] ) continue;
                    const linklet& l = G.EdgeObjectByIndexTo( z, j );
                    links_uni.push( l.nlinks, z2 );    }
               ReverseSort( links_uni );
               z = links_uni[0].second;
               if ( Member( seen, z ) || used[z] ) break;
               seen.push_back(z);    }
          u = seen.back( );

          // Now make the scaffold.

          if (SCAFFOLDING_VERBOSE)
          {    cout << "\ninitial scaffold [" << uscaffolds.size( ) << "]\n";
               cout << "\nusing " << u << "\n";    }
          vec<int> s, sep, dev;
          s.push_back(u);
          used[u] = used[ to_rc[u] ] = True;
          int v = u;
          Bool circ = False;
          while(1)
          {    vec< triple<int,int,int> > links_uni;

               // Find the unipaths that are linked to.
               
               if (SCAFFOLDING_VERBOSE) cout << "\nlooking from " << v << "\n";
               vec<int> ahead = G.From(v);
               Sort(ahead);

               // Now go through more carefully.

               for ( int j = 0; j < G.From(v).isize( ); j++ )
               {    int u2 = G.From(v)[j];
                    if (SCAFFOLDING_VERBOSE) cout << "- see " << u2 << "\n";
                    const linklet& l = G.EdgeObjectByIndexFrom( v, j );

                    // If we've already seen u2 and it's not at the beginning,
                    // then we might have hit a heterozygous inversion.  Skip it.

                    if ( used[u2] && u2 != u ) 
                    {    if (SCAFFOLDING_VERBOSE)
                              cout << "- used and not at beginning, rejecting\n";
                         continue;    }

                    // If we link to both u2 and its reverse complement, then we
                    // might have hit a heterozygous inversion.  Skip it.

                    const int min_rc_ratio = 10;
                    int links_u2 = l.nlinks, links_u2_rc = 0;
                    for ( int jx = 0; jx < G.From(v).isize( ); jx++ )
                    {    if ( G.From(v)[jx] == to_rc[u2] ) 
                         {    links_u2_rc = G.EdgeObjectByIndexFrom( 
                                   v, jx ).nlinks;    }    }
                    if ( links_u2_rc * min_rc_ratio >= links_u2 )
                    {    if (SCAFFOLDING_VERBOSE)
                         {    cout << "- have link to its reverse complement too, "
                                   << "rejecting\n";    }
                         continue;    }

                    // Save link.

                    links_uni.push( l.nlinks, u2, j );    }
               if ( links_uni.empty( ) ) break;
               ReverseSort( links_uni );
               int w = links_uni[0].second;
               if (SCAFFOLDING_VERBOSE) cout << "\ninitially using " << w << "\n";

               // If there are unipaths between v and w, we want to use them.
               // In fact we should probably iterate on this, but don't at present.

               for ( int j = 1; j < links_uni.isize( ); j++ )
               {    int x = links_uni[j].second;
                    Bool have_xw = False;
                    for ( int k = 0; k < G.From(x).isize( ); k++ )
                    {    if ( G.From(x)[k] == w )
                         {    have_xw = True;
                              w = x;
                              if (SCAFFOLDING_VERBOSE) 
                                   cout << "replacing by " << w << "\n";
                              break;    }    }
                    if (have_xw) break;    }

               if ( used[w] )
               {    if ( w != u ) break;
                    circ = True;
                    if (SCAFFOLDING_VERBOSE) cout << "\ncircular\n";    }
               if (SCAFFOLDING_VERBOSE) cout << "\nusing " << w << "\n";
               s.push_back(w);
               const linklet& l = G.EdgeObjectByIndexFrom( v, links_uni[0].third );
               sep.push_back(l.sep), dev.push_back(l.dev);
               if (SCAFFOLDING_VERBOSE) PRINT2( l.sep, l.dev );
               used[w] = used[ to_rc[w] ] = True;
               if (circ) break;
               v = w;    }
          superb x;
          x.SetNtigs( s.size( ) );
          for ( int j = 0; j < s.isize( ); j++ )
          {    x.SetTig( j, s[j] );
               x.SetLen( j, unibases[ s[j] ].size( ) );
               if ( j < s.isize( ) - 1 )
               {    x.SetGap( j, sep[j] ), x.SetDev( j, dev[j] );    }    }
          uscaffolds.push_back(x);
          circled.push_back(circ);    }    }



// Alternative version


void GetNonTrivialEnds( const digraphE<linklet>& G, vec<int>& sources, vec<int>& sinks ){
  vec<int> all;
  G.Sources( all );
  sources.clear();
  for ( size_t is = 0; is < all.size(); is++ )
    if ( G.From(all[is]).size() > 0 )
      sources.push_back( all[is] );
  all.clear();
  G.Sinks( all);
  sinks.clear();
  for ( size_t is = 0; is < all.size(); is++ )
    if ( G.To(all[is]).size() > 0 )
      sinks.push_back( all[is] );
  return;
}



void UnipathScaffoldAlt(
			
			// input:
			
			const digraphE<linklet>& G_orig,
			const vecbasevector& unibases,
			const int K,
			const vec<int>& to_rc,
			const vec<Bool>& exclude,
			const int min_kmers,
			
			// logging and debugging:
			
			const int FORCE_FIRST,
			const Bool SCAFFOLDING_VERBOSE,
			
			// output:
			
			vec<superb>& uscaffolds,
			vec<Bool>& circled ) 
{
  // Set up for building unipath scaffolds.  
  // Allow for the case of circular scaffolds.
  
  digraphE<linklet> G(G_orig);

  // Simplify Graph
  
  cout << "Initially: "; PRINT( G.Edges().size() );
  
  cout << Date() << ": Finding cells" << endl;
  vec< vec<int> > cells;
  FindCells( G, 20, cells );
  cout << Date() << ": found " << cells.size() << " cells" << endl;
  vec<int> to_delete_all;
  for ( size_t ci = 0; ci < cells.size(); ci++ ){
    int v1 = cells[ci].front(), v2 = cells[ci].back();
    vec< vec<int> > paths;
    G.AllPaths( v1, v2, paths );
    vec<int> seqlens( paths.size(), 0 );
    for ( size_t pi = 0; pi < paths.size(); pi++ ){
      for ( size_t vi = 0; vi < paths[pi].size(); vi++ )
	seqlens[pi] += unibases[vi].size();   
    }
    ReverseSortSync( seqlens, paths );
    vec<int> all = cells[ci]; 
    UniqueSort(all);
    vec<int> goodEdges;
    for ( int i = 0; i < paths[0].isize() -1; i++ ){
      int j = i + 1;
      int v = paths[0][i], w = paths[0][j];
      vec<int> locEdges = G.EdgesBetween( v, w );
      goodEdges.append( locEdges );
    }
    UniqueSort( goodEdges );

    vec<int> allEdges = G.EdgesBetween( all );
    UniqueSort( allEdges );
    
    vec<int> to_delete;
    for ( size_t i = 0; i < allEdges.size(); i++ ){
      if ( ! BinMember( goodEdges, allEdges[i] ) )
	to_delete.push_back( allEdges[i] );
    }
    UniqueSort( to_delete );
    to_delete_all.append( to_delete );
  }

  UniqueSort( to_delete_all );
  G.DeleteEdges( to_delete_all );



  vec<int> sources, sinks;
  GetNonTrivialEnds( G, sources, sinks );
  
  vec< vec<int> > ocomps;
  G.Components(ocomps);
  
  cout << Date() << ": initial number of non-trivial ends: sources= " << sources.size() 
       << " sinks= " << sinks.size() << endl;
  vec<int> toRemove;
  for ( size_t i = 0; i < sources.size(); i++ ){
    int vs = sources[i];
    vec<int> vns = G.From(vs);
    for ( size_t j = 0; j < vns.size(); j++ ){
      int vn = vns[j];
      vec<int> to_vn;
      G.GetPredecessors1( vn, to_vn );
      UniqueSort( to_vn );
      
      vec< vec<int> > paths;
      G.AllPaths(vs, vn, paths );
      vec<int> pathsvtxs;
      for ( size_t k = 0; k < paths.size(); k++ )
	pathsvtxs.append( paths[k] );
      UniqueSort( pathsvtxs );
      
      Bool Remove = False;
      if ( to_vn.size() > 5 *  pathsvtxs.size() ){
	Remove = True;
	toRemove.push_back(vs);
      }
      if ( Remove ) break;
    }
  }
  
  for ( size_t i = 0; i < sinks.size(); i++ ){
    int vs = sinks[i];
    vec<int> vps = G.To(vs);
    for ( size_t j = 0; j < vps.size(); j++ ){
      int vp = vps[j];
      vec<int> from_vp;
      G.GetSuccessors1( vp, from_vp );
      UniqueSort( from_vp );
      
      vec< vec<int> > paths;
      G.AllPaths(vp, vs, paths );
      vec<int> pathsvtxs;
      for ( size_t k = 0; k < paths.size(); k++ )
	pathsvtxs.append( paths[k] );
      UniqueSort( pathsvtxs );
      
      Bool Remove = False;
      if ( from_vp.size() > 5 *  pathsvtxs.size() ){
	Remove = True;
	toRemove.push_back(vs);
      }
      if ( Remove ) break;
    }
  }
  UniqueSort( toRemove );
  PRINT( toRemove.size() );
  toRemove.clear();
  // marking affected components
  
  vec<Bool> affectedComps( ocomps.size(), False );
  for ( size_t ir = 0; ir < toRemove.size(); ir++ ){
    for ( size_t ic = 0; ic < ocomps.size(); ic++ ){
      if ( BinMember( ocomps[ic], toRemove[ir] ) )
	affectedComps[ic] = True;
    }
  }
  PRINT2( affectedComps.size(), Sum(affectedComps) );
  // Removing edges at spurious sinks and sources

  for ( size_t i = 0; i < toRemove.size(); i++ )
    G.DeleteEdgesAtVertex( toRemove[i] );

  sources.clear();
  sinks.clear();
  GetNonTrivialEnds( G, sources, sinks );

  cout << Date() << ": final number of non-trivial ends: sources= " << sources.size() 
       << " sinks= " << sinks.size() << endl;
  
  // create dots of affected components
  
  cout << Date() << ": writing dot files" << endl;
  String dot_dir = "dots";
  Mkdir777( dot_dir );
  {
    vec<double> tmp( G.N(), 1.0 );
    String fname1 = dot_dir + "/graph_o.dot";
    ofstream gout( fname1.c_str() );
    G_orig.PrettyDOT( gout, tmp, False, True, False, 0, False, 0, 0, 0 );
    gout.close();
    String fname2 = dot_dir + "/graph_m.dot";
    gout.open( fname2.c_str() );
    G.PrettyDOT( gout, tmp, False, True, False, 0, False, 0, 0, 0 );
    gout.close();
  }
  for ( size_t ic = 0; ic < affectedComps.size(); ic++ ){
    if ( ! affectedComps[ic] ) continue;
    String fname1 = dot_dir + "/graph_o_C" + ToString(ic) + ".dot";
    String fname2 = dot_dir + "/graph_m_C" + ToString(ic) + ".dot";
    vec<int> verts = ocomps[ic];
    vec<double> tmp( G.N(), 1.0 );
    ofstream gout( fname1.c_str() );
    G_orig.PrettyDOT( gout, tmp, False, True, False, 0, False, 0, 0, &verts );
    gout.close();
    //SystemSucceed("dot -Tpng -o " + fname1.SafeBefore(".dot") + ".png " + fname1 );
    gout.open( fname2.c_str() );
    G.PrettyDOT( gout, tmp, False, True, False, 0, False, 0, 0, &verts );
    gout.close();
    //SystemSucceed("dot -Tpng -o " + fname2.SafeBefore(".dot") + ".png " + fname2 );
  }


  // Preparing data for building scaffolds
  
  int nuni = unibases.size( );
  vec<Bool> used = exclude;
  vec< pair<int,int> > len_uni;
  for ( int u = 0; u < nuni; u++ ){    
    int nkmers = unibases[u].isize( ) - K + 1;
    if ( !exclude[u] && nkmers >= min_kmers ) len_uni.push( nkmers, u );    
  }
  ReverseSort(len_uni);

  // Process FORCE_FIRST.
  
  if ( FORCE_FIRST >= 0 ){    
    int jx;
    for ( jx = 0; jx < len_uni.isize( ); jx++ )
      if ( len_uni[jx].second == FORCE_FIRST ) break;
    ForceAssertLt( jx, len_uni.isize( ) );
    swap( len_uni[0], len_uni[jx] );    
  }

  // Build unipath scaffolds.

  cout << Date() << ": building initial scaffolds" << endl;
  if (SCAFFOLDING_VERBOSE) cout << "\nINITIAL SCAFFOLDS\n";
  for ( int ui = 0; ui < len_uni.isize( ); ui++ ){    
  //for ( int ui = len_uni.isize( ) -1; ui >= 0;  ui-- ){    
    int u0 = len_uni[ui].second;
    if ( used[u0] ) continue;
    
    
    // Now make the scaffold.
    
    if (SCAFFOLDING_VERBOSE){    
      cout << "\ninitial scaffold [" << uscaffolds.size( ) << "]\n";
      cout << "\nusing " << u0 << "\n";    
    }
    vec<int> s, sep, dev;
    s.push_back(u0);
    used[u0] = used[ to_rc[u0] ] = True;
    Bool circ = False;

    for ( int pass = 1; pass <= 2; pass++ ){
      
      int v = u0;
      while(1){    
	vec< triple<int,int,int> > links_uni;
	
	// Find the unipaths that are linked to.
	
	if (SCAFFOLDING_VERBOSE) cout << "\nlooking from " << v << "\n";
	vec<int> ahead = pass == 1 ? G.From(v) : G.To(v);
	Sort(ahead);
	
	// Now go through more carefully.
	
	for ( int j = 0; j < (pass == 1 ? G.From(v).isize( ) : G.To(v).isize()); j++ ){ 
	  int u2 = (pass == 1 ? G.From(v)[j] : G.To(v)[j] );
	  if (SCAFFOLDING_VERBOSE) cout << "- see " << u2 << "\n";
	  const linklet& l = (pass == 1 ? G.EdgeObjectByIndexFrom( v, j ) : G.EdgeObjectByIndexTo( v, j ) );
	  
	  // If we've already seen u2 and it's not at the beginning,
	  // then we might have hit a heterozygous inversion.  Skip it.
	  
	  if ( used[u2] && u2 != u0 ) {    
	    if (SCAFFOLDING_VERBOSE)
	      cout << "- used and not at beginning, rejecting\n";
	    continue;    
	  }
	  
	  // If we link to both u2 and its reverse complement, then we
	  // might have hit a heterozygous inversion.  Skip it.
	  
	  const int min_rc_ratio = 10;
	  int links_u2 = l.nlinks, links_u2_rc = 0;
	  for ( int jx = 0; jx < (pass == 1 ? G.From(v).isize( ) : G.To(v).isize() ); jx++ ){    
	    if ( (pass == 1 ? G.From(v)[jx] : G.To(v)[jx]) == to_rc[u2] ) {    
	      links_u2_rc = (pass == 1 ? G.EdgeObjectByIndexFrom( v, jx ).nlinks : G.EdgeObjectByIndexTo( v, jx ).nlinks);    
	    }    
	  }
	  if ( links_u2_rc * min_rc_ratio >= links_u2 ){    
	    if (SCAFFOLDING_VERBOSE){    
	      cout << "- have link to its reverse complement too, "
		   << "rejecting\n";    
	    }
	    continue;    
	  }
	  
	  // Save link.
	  
	  links_uni.push( l.nlinks, u2, j );    
	}
	if ( links_uni.empty( ) ) break;
	ReverseSort( links_uni );
	int w = links_uni[0].second;
	if (SCAFFOLDING_VERBOSE) cout << "\ninitially using " << w << "\n";
	
	// If there are unipaths between v and w, we want to use them.
	// In fact we should probably iterate on this, but don't at present.
	
	for ( int j = 1; j < links_uni.isize( ); j++ ){   
	  int x = links_uni[j].second;
	  Bool have_xw = False;
	  for ( int k = 0; k < (pass == 1 ? G.From(x).isize( ) : G.To(x).isize() ); k++ ){
	    if ( (pass == 1 ? G.From(x)[k] : G.To(x)[k] ) == w ) {    
	      have_xw = True;
	      w = x;
	      if (SCAFFOLDING_VERBOSE) 
		cout << "replacing by " << w << "\n";
	      break;    
	    }    
	  }
	  if (have_xw) break;    
	}
	
	if ( used[w] ){    
	  break;
	  //if ( w != u0 ) break;
	  //circ = True;
	  //if (SCAFFOLDING_VERBOSE) cout << "\ncircular\n";    
	}
	if (SCAFFOLDING_VERBOSE) cout << "\nusing " << w << "\n";
	pass == 1 ? s.push_back(w) : s.push_front(w);
	const linklet& l = 
	  (pass == 1 ? G.EdgeObjectByIndexFrom( v, links_uni[0].third ) : G.EdgeObjectByIndexTo( v, links_uni[0].third ) );
	pass == 1 ? sep.push_back(l.sep), dev.push_back(l.dev) : sep.push_front(l.sep), dev.push_front(l.dev);
	if (SCAFFOLDING_VERBOSE) PRINT2( l.sep, l.dev );
	used[w] = used[ to_rc[w] ] = True;
	if (circ) break;
	v = w;      
      }
      
    }

    Bool FoundSelf = False;
    superb x;
    x.SetNtigs( s.size( ) );
    for ( int j = 0; j < s.isize( ); j++ ){ 
      if ( s[j] == u0 ) FoundSelf = True;
      x.SetTig( j, s[j] );
      x.SetLen( j, unibases[ s[j] ].size( ) );
      if ( j < s.isize( ) - 1 ){    
	x.SetGap( j, sep[j] ), x.SetDev( j, dev[j] );    
      }    
    }
    if ( ! FoundSelf ){
      cout << "WARNING: starting unipath " << u0 << " not found in corresponding scaffold!!!!!!!\n";
    }
    //if ( x.Ntigs() == 1 && G.From(x.Tig(0)).size() == 0 && G.To(x.Tig(0)).size() == 0 )
    //continue;
    uscaffolds.push_back(x);
    circled.push_back(circ);    
  }

  int nCompleteEnds = 0, nWhole = 0, nComplete = 0, nTigs = 0, nTigsComplex = 0, nScaffComplex = 0;
  for ( size_t si = 0; si < uscaffolds.size(); si++ ){
    nTigs += uscaffolds[si].Ntigs();
    if ( uscaffolds[si].Ntigs() > 1 ){
      nTigsComplex += uscaffolds[si].Ntigs();
      nScaffComplex++;
    }
    int tigFirst = uscaffolds[si].Tig(0);
    int tigLast = uscaffolds[si].Tig( uscaffolds[si].Ntigs() -1 );
    if ( G.To( tigFirst ).size() == 0 && G.From( tigLast ).size() == 0 ){
      nCompleteEnds++;
      //ForceAssertNe( uscaffolds[si].Ntigs(), 1 );
    }
  }
  PRINT5( uscaffolds.size(), nTigs, nTigsComplex, nScaffComplex, nCompleteEnds );
}
