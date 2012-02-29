///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Fastavector.h"
#include "Set.h"
#include "Vec.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "graph/FindCells.h"
#include "paths/HyperFastavector.h"
#include "paths/CLinFastavec.h"


String error_string(const int error_code)
{
  
  if (error_code == 1) return "consensus.size() = 0";
  if (error_code == 2) return "too many ambiguities";
  return "undefined error";
}



/**
 * CLinFastavec
 * Constructor
 */
CLinFastavec::CLinFastavec( HyperFastavector &hyper, ostream *plog ) :
  hyper_ ( hyper ),
  base_inter_save_ ( "" ),
  log_ ( plog ),
  iter_ ( 0 ),
  new_algorithm_(False)
{ }

/**
 * CLinFastavec
 * SetBaseInterSave
 *
 * This will turn on saving all intermediate steps. Also, save the
 * initial step (unmodified).
 */
void CLinFastavec::SetBaseInterSave( const String &base_name )
{
  base_inter_save_ = base_name;
  this->SaveIntermediate( );
}

/**
 * CLinFastavec
 * FlattenCells
 *
 * Return the number of flattening events
 */
int CLinFastavec::FlattenCells( const int MAX_CELL_SIZE )
{
  Bool verbose = False;

  Bool PRINT_CELLS = False;

  vec< vec<int> > cells;
  FindCells( hyper_, MAX_CELL_SIZE, cells );
  vec<int> to_left;
  hyper_.ToLeft(to_left);

  // Fix cells, compensating for putative bug.  Some contain duplicates.
  // This fix probably isn't quite right.

  for ( int i = 0; i < cells.isize( ); i++ )
  {    vec<int> c;
       for ( int j = 0; j < cells[i].isize( ); j++ )
       {    if ( j == cells[i].isize( ) - 1 || !Member( c, cells[i][j] ) )
                 c.push_back( cells[i][j] );    }
       cells[i] = c;    }

  // Counter for the events.
  int n_events = 0;

  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log_ ? *log_ : devnull;
 
  out << "Flattening cells..." << endl;

  // Eliminate nonminimal cells.
  {
    vec<Bool> nonminimal( cells.size( ), False );
    vec< vec<int> > cell_index( hyper_.N( ) );
    for ( int i = 0; i < cells.isize( ); i++ )
      for ( int j = 0; j < cells[i].isize( ); j++ )
	cell_index[ cells[i][j] ].push_back(i);
      
    for ( int i = 0; i < cells.isize( ); i++ ) {
      vec<int> friends;
      for ( int j = 0; j < cells[i].isize( ); j++ )
	friends.append( cell_index[ cells[i][j] ] );
      UniqueSort(friends);
      friends.EraseValue(i);
      for ( int m = 0; m < friends.isize( ); m++ ) {
	if ( Subset( cells[ friends[m] ], cells[i] ) ) {
	  nonminimal[i] = True;
	  break;
	}
      }
    }
    
    EraseIf( cells, nonminimal );
  }

  if (PRINT_CELLS)
  {    for ( int i = 0; i < cells.isize( ); i++ )
       {    cout << "\ncell " << i << ":";
            for ( int j = 0; j < cells[i].isize( ) - 1; j++ )
            {    int v = cells[i][j];
                 for ( int k = 0; k < hyper_.From(v).isize( ); k++ )
                 {    cout << " " 
                           << BaseAlpha( hyper_.EdgeObjectIndexByIndexFrom( v, k ) );
                                   }    }
            cout << "\n";    }    }

  // Try to collapse the cells.

  for ( int i = 0; i < cells.isize( ); i++ ) {
    if (PRINT_CELLS) cout << "\ntrying cell " << i << "\n";
    vec< vec<int> > paths;
    hyper_.AllPaths( cells[i].front( ), cells[i].back( ), paths );

    // Don't use cells that contain vertices that are not between the first and last.
   
    vec<int> used;
    for ( int j = 0; j < paths.isize( ); j++ )
         used.append( paths[j] );
    UniqueSort(used);
    if ( used.size( ) < cells[i].size( ) ) 
    {    if (PRINT_CELLS) cout << "fails because of internal vertices\n";
         continue;    }

    // WHAT ABOUT CELLS THAT CONTAIN LOOPS?  OR MORE COMPLEX CYCLES?

    // SANTEMP
    if ( paths.size( ) < 1 ) {
//       cout << "ERROR - paths is empty" << endl;
//       cout << "cell: ";
//       for (int jj=0; jj<cells[i].isize( ); jj++)
// 	cout << cells[i][jj] << " ";
//       cout << "\n" << endl;
      if (PRINT_CELLS) cout << "fails because no paths\n";
      continue;
    }
    
    fastavector con = this->ConsensusPath( paths, True );
    int error_code = 0;
    if ( con.size( ) == 0 ) error_code = 1;
    if ( 10U * con.AmbCount( ) > con.size( ) ) error_code = 2;
    //if ( con.MaxWindowAmbCount(100) > 10 ) error_code = 2; 

    int cell1 = paths[0][0];
    int cell2 = paths[0].back( );
    if (verbose) out << "\n";
    out << "g." << iter_ << " - "
	<< "cell v[" << cell1 << "," << cell2 << "] "
	<< "(" << paths.size( ) << " paths)";
    if ( error_code == 0 ) out << "\n";
    else out << " - FAILED (" << error_string(error_code) << ")\n";
    if (PRINT_CELLS) PRINT(error_code);

    if (verbose)
    {    out << "vertices:";
         for ( int j = 0; j < cells[i].isize( ); j++ )
              out << " " << cells[i][j];
         out << endl;
         for ( int j = 0; j < cells[i].isize( ); j++ )
         {    int v = cells[i][j];
              for ( int k = 0; k < hyper_.From(v).isize( ); k++ )
              {    int w = hyper_.From(v)[k];
                   if ( !Member( cells[i], w ) ) continue;
                   out << "edge " << v << " --> " << w 
                       << " [" << hyper_.EdgeObjectByIndexFrom( v, k ).size( )
                       << " bases]" << endl;    }    }    }

  vec<fastavector> tails;
  for ( int l = 0; l < paths.isize( ); l++ ) {
    fastavector f;
    for ( int r = 0; r < paths[l].isize( ) - 1; r++ ) {
      int v = paths[l][r], w = paths[l][r+1];
      vec<fastavector> vw = hyper_.EdgeObjectsBetween( v, w );
      fastavector pop = this->CombineSlack( vw );
      // if ( pop.size( ) == 0 ) return fastavector( );   // popping failed
      // Originally the following 'if' tested for f.size( ) > 0, which probably
      // makes sense, but on rare occasions would cause horrible failures because
      // the resize would be to a negative number.  Changed to avoid this
      // catastrophe.
      if ( (int) f.size( ) - hyper_.K( ) + 1 >= 0 ) 
           f.resize( f.size( ) - hyper_.K( ) + 1 );
      f.append( pop.begin(), pop.end() );
    }
    tails.push_back(f);
  }

  if (verbose)
  {    out << "path sizes ";
       for ( int ii = 0; ii < tails.isize( ); ii++ )
            out << tails[ii].size( ) << " ";
       out << " --> consensus size " << con.size( ) << endl;    }

    if ( error_code > 0  ) continue;

    vec<int> to_delete;
    for ( int j = 0; j < cells[i].isize( ) - 1; j++ )
    { for ( int z = 0; z < hyper_.FromEdgeObj( cells[i][j] ).isize( ); z++ )
      {    int w = hyper_.FromEdgeObj( cells[i][j] )[z];
           if (verbose)
           {    out << "deleting edge of size " 
                    << hyper_.EdgeObject(w).size( ) << endl;    }    }
      to_delete.append( hyper_.FromEdgeObj( cells[i][j] ) );
      }
    Sort(to_delete);
    ForceAssertGt( to_delete.size(), 1u );
    hyper_.DeleteEdges( to_delete, to_left );
    hyper_.AddEdge( cells[i].front( ), cells[i].back( ), con );
    to_left.push_back(cells[i].front());
    n_events++;
  }

  out << "Flattened " << n_events << " cells." << endl;  

  this->CleanUp( );
  return n_events;
}

/**
 * CLinFastavec
 * FlattenFrayedEnds
 *
 * Return the number of flattening events.
 */
int CLinFastavec::FlattenFrayedEnds( )
{
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log_ ? *log_ : devnull;

  out << "Flattening frayed ends..." << endl;


  // Events counter.
  int n_events = 0;

  // Two passes.
  vec< vec<int> > compv;
  hyper_.ComponentsAlt(compv);
  for ( int pass = 1; pass <= 2; pass++ ) {
    hyper_.Reverse( );
    vec<int> to_left;
    hyper_.ToLeft(to_left);
    
    for ( int j = 0; j < compv.isize( ); j++ ) {
      vec<int> verts = compv[j], sinks;
      for ( int x = 0; x < verts.isize( ); x++ )
	if ( hyper_.Sink( verts[x] ) ) sinks.push_back( verts[x] );
      if ( sinks.size( ) <= 1 ) continue;
      
      // Find the vertices that are predecessors of every sink.  Then
      // amongst such vertices v, find those for which predecessors(v)
      // and successors(v) are connected to each other only through v.
      // Then amoungst these, determine if there is a single vertex
      // that is a successor of all the others.
      vec< vec<int> > predx( sinks.size( ) );
      for ( int x = 0; x < sinks.isize( ); x++ )
	hyper_.GetPredecessors1( sinks[x], predx[x] );
      vec<int> pred, goodpred, bestpred;
      Intersection( predx, pred );
      for ( int k = 0; k < pred.isize( ); k++ ) {
	vec<int> P, S;
	hyper_.GetPredecessors1( pred[k], P );
	hyper_.GetSuccessors1( pred[k], S );
	P.EraseValue( pred[k] ), S.EraseValue( pred[k] );
	
	Bool illegal_connection = False;
	for ( int l = 0; l < P.isize( ); l++ )
	  for ( int u = 0; u < hyper_.From( P[l] ).isize( ); u++ )
	    if ( BinMember( S, hyper_.From( P[l] )[u] ) )
	      illegal_connection = True;
	if ( !illegal_connection ) goodpred.push_back( pred[k] );
      }
      if ( goodpred.empty( ) ) continue;
      
      vec< vec<int> > suc( goodpred.isize( ) );
      for ( int k = 0; k < goodpred.isize( ); k++ )
	hyper_.GetSuccessors1( goodpred[k], suc[k] );
      vec<int> sucall, sucallx;
      Intersection( suc, sucall );
      Intersection( sucall, goodpred, sucallx );
      if ( !sucallx.solo( ) ) continue;
      int y = sucallx[0];
      
      // Define the tail.
      vec<int> s;
      hyper_.GetSuccessors1( y, s );
      
      // Check for loops in the tail.
      Bool looped = False;
      for ( int l = 0; l < s.isize( ); l++ ) {
	if ( hyper_.LoopAt( s[l] ) ) {
	  looped = True;
	  break;
	}
      }
      if (looped) continue;
      
      // Define the sequence for the tail.
      vec< vec<int> > paths;
      const int maxpaths = 1000;
      Bool AllPaths_OK = hyper_.AllPaths( y, -1, paths, maxpaths );
      if ( !AllPaths_OK ) continue;
      fastavector con = ConsensusPath( paths, False );

      int error_code = 0;
      if ( con.size( ) == 0 ) error_code = 1;
      if ( 10U * con.AmbCount( ) > con.size( ) ) error_code = 2;
      //if ( con.MaxWindowAmbCount(100) > 10 ) error_code = 2; 
     
      out << "g." << iter_ << " - "
	  << " frayed_" << pass << " v[" << y << "] "
	  << "( " << paths.size( ) << " paths)";
      if ( error_code == 0 ) out << "\n";
      else out << " - FAILED (" << error_string(error_code) << ")\n";
      
      if ( error_code > 0  ) continue;
      
      // Operate on the graph.
      vec<int> to_delete;
      for ( int k = 0; k < s.isize( ); k++ )
	to_delete.append( hyper_.FromEdgeObj( s[k] ) );
      Sort(to_delete);
      hyper_.DeleteEdges( to_delete, to_left );
      int nv = hyper_.N( );
      hyper_.AddVertices(1);
      hyper_.AddEdge( y, nv, con );
      to_left.push_back(y);
      n_events++;
    }
  }

  out << "Flattened " << n_events << " frayed ends." << endl;  

  this->CleanUp( true );
  return n_events;
}

/**
 * CLinFastavec
 * FlattenBubbles
 *
 * Return the number of flattening events.
 */
int CLinFastavec::FlattenBubbles( )
{
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log_ ? *log_ : devnull;

  out << "Flattening bubbles..." << endl;

  vec<int> to_left;
  hyper_.ToLeft(to_left);

  // Events counter.
  int n_events = 0;

  // Find and flatten bubbles.
  vec< pair<int,int> > bubbles;
  this->FindBubbles( bubbles );
  
  for (uint ii=0; ii<bubbles.size( ); ii++) {
    int v1 = bubbles[ii].first;
    int v2 = bubbles[ii].second;

    vec<fastavector> v1v2 = hyper_.EdgeObjectsBetween( v1, v2 );
    fastavector con = this->CombineSlack( v1v2 );

    int error_code = 0;
    if ( con.size( ) == 0 ) error_code = 1;
    if ( 10U * con.AmbCount( ) > con.size( ) ) error_code = 2;
    //if ( con.MaxWindowAmbCount(100) > 10 ) error_code = 2; 

    out << "g." << iter_ << " - "
	<< "bubble v[" << v1 << "," << v2 << "]";
    if ( error_code == 0 ) out << "\n";
    else out << " - FAILED (" << error_string(error_code) << ")\n";

    if ( error_code > 0  ) continue;
    
    vec<int> to_delete = hyper_.EdgesBetween( v1, v2 );
    Sort( to_delete );
    hyper_.DeleteEdges( to_delete, to_left );
    hyper_.AddEdge( v1, v2, con );
    to_left.push_back(v1);
    n_events++;
  }
  
  out << "Flattened " << n_events << " bubbles..." << endl;

  this->CleanUp( );
  return n_events;
}

/**
 * CLinFastavec
 * ReportBubbles
 */
void CLinFastavec::ReportBubbles( )
{
  if ( ! log_ ) return;
  ostream &out = *log_;

  vec< pair<int,int> > bubbles;
  this->FindBubbles( bubbles );

  out << "REPORTING EXISTING UNPOPPED BUBBLES:\n";
  for (int ii=0; ii<bubbles.isize( ); ii++) {
    int v1 = bubbles[ii].first;
    int v2 = bubbles[ii].second;
    vec<int> edges = hyper_.EdgesBetween( v1, v2 );
    
    out << "v" << v1 << "-v" << v2 << ":   ";
    for (int jj=0; jj<edges.isize( ); jj++)
      out << "e" << edges[jj] << " ("
	  << hyper_.EdgeObject( edges[jj] ).size( ) << "),  ";
    out << "\n";
  }
  out << endl;
  
}

/**
 * CLinFastavec
 * SaveIntermediate
 */
void CLinFastavec::SaveIntermediate( ) const
{
  if ( base_inter_save_ == "" ) return;
  
  String out_file = base_inter_save_ + "." + ToString( iter_ );
  if ( log_ ) *log_ << "Saving " << out_file << endl;
  int fd = OpenForWrite( out_file );
  BinaryWrite( fd, hyper_ );
}

/**
 * CLinFastavec
 * ConsensusPath
 */
fastavector CLinFastavec::ConsensusPath( vec< vec<int> > &paths,
					 bool same_size ) const
{
  vec<fastavector> tails;
  for ( int l = 0; l < paths.isize( ); l++ ) {
    fastavector f;
    for ( int r = 0; r < paths[l].isize( ) - 1; r++ ) {
      int v = paths[l][r], w = paths[l][r+1];
      vec<fastavector> vw = hyper_.EdgeObjectsBetween( v, w );
      fastavector pop = this->CombineSlack( vw );
      if ( pop.size( ) == 0 ) return fastavector( );   // popping failed
      if ( f.size( ) > 0 ) f.resize( f.size( ) - hyper_.K( ) + 1 );
      f.append( pop.begin(), pop.end() );
    }
    tails.push_back(f);
  }

  fastavector pop;
  if ( same_size ) {
    pop = this->CombineSlack( tails );
    if ( pop.size( ) == 0 ) pop.clear( );   // popping failed
    return pop;
  }
  pop = Combine( tails );

  // New algorithm: if there are twenty or more ambiguities, use a longest path.

  if ( NewAlgorithm( ) )
  {    const int max_ambiguities = 20;
       if ( pop.AmbCount( ) > max_ambiguities )
       {    int M = 0, best = -1;
            for ( int j = 0; j < tails.isize( ); j++ )
            {    if ( (int) tails[j].size( ) > M )
                 {    M = tails[j].size( );
                      best = j;    }    }
            ForceAssertGe( best, 0 );
            pop = tails[best];    }    }

  return pop;

}

/**
 * CLinFastavec
 * CombineSlack
 *
 * Wrapper around Combine( ) in Fastavector.h. If the fastavectors in
 * fvec are somewhat similar (ie their lengths are not exactly the
 * same, but are not too different either), then the "merge" is done
 * by simply taking the longest entry.
 */
fastavector CLinFastavec::CombineSlack( const vec<fastavector> &fvec ) const
{
  const int max_indels = 200;
  const double max_ratio = 0.1;

  // Nothing to do.
  if ( fvec.size( ) < 1 )
    return fastavector( );
  
  // Cannot merge if sizes are too different.
  vec< pair<int,int> > sizes2id;
  for (int ii=0; ii<fvec.isize( ); ii++)
    sizes2id.push_back( make_pair( fvec[ii].size( ), ii ) );
  sort( sizes2id.begin( ), sizes2id.end( ) );
  int smallest = sizes2id[0].first;
  int largest = sizes2id.back( ).first;
  double ratio = (double)( largest - smallest ) / (double)largest;
  if ( largest - smallest > max_indels && ratio > max_ratio )
    return fastavector( );
  
  // If all sizes match, return combined fastavector (matches old behavior).
  if ( largest == smallest )
    return Combine( fvec );

  // Return longest entry.
  return fvec[sizes2id.back( ).second];
  
}

/**
 * CLinFastavec
 * FindBubbles
 */
void CLinFastavec::FindBubbles( vec< pair<int,int> > &bubbles) const
{
  bubbles.clear( );
  uint n_bubbles = 0;

  // The first iteration counts the number of events.
  for (int iter=0; iter<2; iter++) {
    if ( iter == 1 ) bubbles.reserve( n_bubbles );

    // Loop over all vertices in the digraph.
    for (int v1=0; v1<hyper_.N( ); v1++) {
      
      // All the vertices connected to v1.
      vec<int> successors;
      hyper_.GetSuccessors1( v1, successors );
      
      // Loop over all successors of v1 (only if v1<v2).
      for (int succid=0; succid<successors.isize( ); succid++) {
	int v2 = successors[succid];
	// if ( v2 >= v1 ) continue;
        if ( v2 == v1 ) continue;
	
	// All edges between v1 and v2.
	vec<int> edges = hyper_.EdgesBetween( v1, v2 );
	if ( edges.size( ) != 2 ) continue;

	// A bubble.
	if ( iter == 0 ) n_bubbles++;
	else bubbles.push_back( make_pair( v1, v2 ) );
      }
    }
  }
}

/**
 * CLinFastavec
 * CleanUp
 *
 * Use alt for a second cleaning method. Both are stolen from David's
 * code, it may well be we do not need both. Notice that iter_ will be
 * updated as well.
 */
void CLinFastavec::CleanUp( bool alt )
{
  iter_++;

  if ( base_inter_save_ != "" ) this->SaveIntermediate( );

  if ( log_ ) *log_ << endl;

  if ( alt ) {
    hyper_.RemoveUnneededVertices( );
    hyper_.RemoveEdgelessVertices( );
    hyper_.RemoveDeadEdgeObjects( );
  }
  else {
    hyper_.RemoveUnneededVertices( );
    RemoveHangingEnds2( hyper_ );
    hyper_.RemoveUnneededVertices( );
    hyper_.RemoveDeadEdgeObjects( );
  }
  
}

