///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library PTHREAD
#include <strstream>

#include <pthread.h>

#include "CoreTools.h"
#include "Equiv.h"
#include "graph/DigraphTemplate.h"
#include "paths/HyperFastavector.h"
#include "system/Worklist.h"

// Template instantiations:

template void digraphE<fastavector>::ToLeft(vec<int>&) const;
template void digraphE<fastavector>::ToRight(vec<int>&) const;
template void digraphE<fastavector>::RemoveDeadEdgeObjects();
template void digraphE<fastavector>::ComponentEdges( vec< vec<int> >& edges ) const;
template digraphE<fastavector>::digraphE(digraphE<fastavector> const&, vec<int> const&);

// Order a pair of vec<int>s lexicographically using subset relation.

struct order_vecint_binsubset_pair
     : public binary_function< const pair< vec<int>, vec<int> >&,
          const pair< vec<int>, vec<int> >&, bool >
{    bool operator( ) ( const pair< vec<int> , vec<int> >& left,
          const pair< vec<int>, vec<int> >& right ) const
     {    if ( left.first == right.first )
          {    if ( left.second == right.second ) return false;
               else return BinSubset( left.second, right.second );    }
          else return BinSubset( left.first, right.first );    }
};

void HyperFastavector::SetToDisjointUnionOf( const vec<HyperFastavector>& v )
{    for ( size_t i = 0; i < v.size(); i++ )
          ForceAssertEq( K( ), v[i].K( ) );
     Clear( );
     int nvert = 0, nedge = 0;
     for ( size_t i = 0; i < v.size(); i++ )
     {    nvert += v[i].N( );
          nedge += v[i].EdgeObjectCount( );    }
     FromMutable( ).reserve(nvert), ToMutable( ).reserve(nvert);
     FromEdgeObjMutable( ).reserve(nvert), ToEdgeObjMutable( ).reserve(nvert);
     EdgesMutable( ).reserve(nedge);
     int vcount = 0, ecount = 0;
     for ( size_t i = 0; i < v.size(); i++ )
     {    const HyperFastavector& h = v[i];
          EdgesMutable( ).append( h.Edges( ) );
          vec< vec<int> > from = h.From( ), to = h.To( );
          vec< vec<int> > frome = h.FromEdgeObj( ), toe = h.ToEdgeObj( );
          for ( int j = 0; j < h.N( ); j++ )
          {    for ( size_t u = 0; u < from[j].size(); u++ )
               {    from[j][u] += vcount;
                    frome[j][u] += ecount;    }
               for ( size_t u = 0; u < to[j].size(); u++ )
               {    to[j][u] += vcount;
                    toe[j][u] += ecount;    }    }
          FromMutable( ).append(from), ToMutable( ).append(to);
          FromEdgeObjMutable( ).append(frome), ToEdgeObjMutable( ).append(toe);
          vcount += h.N( );
          ecount += h.EdgeObjectCount( );    }    }

void HyperFastavector::Reverse( )
{    for ( int i = 0; i < EdgeObjectCount( ); i++ )
          EdgeObjectMutable(i).ReverseComplement( );
     digraphE<fastavector>::Reverse( );    }

class ScaffoldingArgs {

     public:

     ScaffoldingArgs( ) { }

     vec< vec<int> > comp;
     vec<int> to_left, to_right;
     const HyperFastavector* H;
     vec<superb> superbs;
     vec< pair<int,int> > superb_indices;

};

// Attempt to linearize this HyperFastavector into a set of superb objects.
// This process involves many heuristics which are not thoroughly optimized.
  
// Define heuristic constants.

static int MIN_EDGE_TO_SAVE_GLOBAL = 1000;
static const unsigned int MAX_EDGE_TO_DISCARD = 10000;
static const int UNKNOWN_GAP = 1000;
static const int UNKNOWN_DEV = 1000;

static ScaffoldingArgs S;

static pthread_mutex_t scaffold_component_lock = PTHREAD_MUTEX_INITIALIZER;
static bool TEMP_DEBUG = false; // TEMP

void ScaffoldComponent( const int i )
{    
     Bool verbose = False; // for now, you have to manually toggle this

     // Maximium number of paths between edges - bail out if exceeded
     const int MAX_PATHS = 10000;
  
     vec<int> edges = S.comp[i];
     Sort(edges);
     int nedges = edges.size( );
     const HyperFastavector* H = S.H;

     // Find all edges that contain a source or a sink, or is a terminal
     // loop, or if removed from the graph, would disconnect it.

     vec<int> critical;
     for ( int ei = 0; ei < nedges; ei++ ) 
     {    int e = edges[ei];
          int v = S.to_left[e], w = S.to_right[e];
          if ( H->To(v).empty( ) ) critical.push_back(e);
          else if ( H->From(w).empty( ) ) critical.push_back(e);
          else if ( v == w && ( H->From(v).solo( ) || H->To(v).solo( ) ) )
	       critical.push_back(e);
          else if ( v == w ) continue;
          else
          {    vec<int> to_process;
               set<int> component;
               int ep = H->EdgeObjectIndexByIndexTo( v, 0 );
               to_process.push_back(ep), component.insert(ep);
               while( to_process.nonempty( ) )
               {    int ex = to_process.back( );
                    to_process.pop_back( );
	            int vx = S.to_left[ex], wx = S.to_right[ex];
	            for ( size_t l = 0; l < H->To(vx).size(); l++ ) 
                    {    int en = H->EdgeObjectIndexByIndexTo( vx, l );
                         if ( en == e ) continue;
                         if ( Member( component, en ) ) continue;
                         to_process.push_back(en);
                         component.insert(en);    }
	            for ( size_t l = 0; l < H->From(wx).size(); l++ ) 
                    {    int en = H->EdgeObjectIndexByIndexFrom( wx, l );
                         if ( en == e ) continue;
                         if ( Member( component, en ) ) continue;
                         to_process.push_back(en);
                         component.insert(en);    }    }
               if ( (int) component.size( ) < nedges - 1 )
                    critical.push_back(e);    }    }

     // Order the critical edges.

     vec< pair< vec<int>, vec<int> > > nexts( critical.size( ) );
     for ( size_t j = 0; j < critical.size(); j++ ) 
     {    H->GetPredecessors1( S.to_left[ critical[j] ], nexts[j].first );
          H->GetPredecessors1( S.to_right[ critical[j] ], nexts[j].second );    }
     SortSync( nexts, critical, order_vecint_binsubset_pair( ) );

     // Decide if linearization is accepted.

     Bool good = True;
     for ( size_t j = 0; j+1 < nexts.size(); j++ ) 
     {    if ( !BinSubset( nexts[j].first, nexts[j+1].first ) ) good = False;
          if ( !BinSubset( nexts[j].second, nexts[j+1].second ) ) good = False;
          if ( nexts[j] == nexts[j+1] ) good = False;    }
     for ( size_t j = 0; j < edges.size(); j++ ) 
     {    int e = edges[j];
          if ( !Member( critical, e ) 
               && H->EdgeObject(e).size( ) > MAX_EDGE_TO_DISCARD )
          {   good = False;    }    }

     // Take action.

     vec<superb> ss;
     if ( !good ) 
     {
          // Can't linearize.  Salvage by putting edges in separate scaffolds.

          for ( int ei = 0; ei < nedges; ei++ ) 
          {    fastavector e = H->EdgeObject( edges[ei] );
	       if ( (int) e.size( ) >= MIN_EDGE_TO_SAVE_GLOBAL )
               {    superb s;
                    s.SetNtigs(1);
	            s.SetTig( 0, edges[ei] );
	            s.SetLen( 0, e.size( ) );    
                    ss.push_back(s);    }    }    }
     else 
     {    vec<int> L(critical);
          if ( L.empty( ) ) return;
          superb s;
          s.SetNtigs( L.size( ) );
          for ( size_t j = 0; j < L.size(); j++ ) 
          {    fastavector e = H->EdgeObject( L[j] );
	       s.SetTig( j, L[j] );
	       s.SetLen( j, e.size( ) );
	       if ( j+1 < L.size() ) 
               {
	            // Now we need to figure out how big the gap is.  If there
	            // are no cycles, this is easy.  If there are cycles, at
	            // present we punt.  Either way, we should come back and
	            // use read pairing to compute this better.

	            int v = S.to_right[ L[j] ], w = S.to_left[ L[j+1] ];
	            vec< vec<int> > paths;
		    if ( TEMP_DEBUG ) cout << Date() << ": Before AllPaths" << endl;
		    // Set the maximum search steps to 1e8 to prevent AllPaths from running forever
	            // H->AllPaths( v, w, paths, MAX_PATHS );
	            H->AllPaths( v, w, paths, MAX_PATHS, False, 100*1000*1000 );
		    if ( TEMP_DEBUG ) cout << Date() << ": After AllPaths" << endl;
	            vec<int> verts_in_paths;
	            for ( size_t r = 0; r < paths.size(); r++ )
	                 verts_in_paths.append( paths[r] );
	            UniqueSort(verts_in_paths);
	            Bool have_loop = False;
	            for ( size_t r = 0; r < verts_in_paths.size(); r++ ) 
                    {    if ( H->LoopAt( verts_in_paths[r] ) ) 
                         {    have_loop = True;
	                      break;    }    }
	            if (have_loop || paths.empty()) 
                    {    s.SetGap( j, UNKNOWN_GAP );
	                 s.SetDev( j, UNKNOWN_DEV );    }
	            else 
                    {
	                 // ForceAssert( paths.nonempty( ) );

	                 // Find the shortest and longest paths from v to w.

	                 int min_path = 1000000000, max_path = 0;
	                 for ( size_t r = 0; r < paths.size(); r++ ) 
                         {    int path_low = 0, path_high = 0;
	                      for ( size_t u = 0; u+1 < paths[r].size(); u++ ) 
                              {    int x = paths[r][u], y = paths[r][u+1];
     
		                   // Find the min and max length of an edge
		                   // from x to y.
     
		                   unsigned int m = 100000000, M = 0;
		                   for (size_t q = 0; q < H->From(x).size(); q++) 
                                   {    if ( H->From(x)[q] == y ) 
                                        {    const fastavector& p =
		                                  H->EdgeObjectByIndexFrom( x, q );
		                             m = Min( m, p.size( ) - H->K( ) + 1 );
		                             M = Max( M, p.size( ) - H->K( ) + 1 );
		                                  }    }

		                   path_low += m, path_high += M;    }
	                      min_path = Min( min_path, path_low );
	                      max_path = Max( max_path, path_high );    }

	                 // Define the gap.

	                 int sep = ( min_path + max_path ) / 2 - H->K( ) + 1;
	                 int dev = ( max_path - min_path ) / 2;
	                 if ( paths.empty( ) ) sep = 1000, dev = 1000;
	                 s.SetGap( j, sep );
	                 s.SetDev( j, dev );    }    }    }
          ss.push_back(s);

     // Report what's happening.

     if (verbose)
     {    
          // Should be single threaded.

          {    cout << "\nscaffolding component " << i << endl;
               cout << "edges =";
               for ( int j = 0; j < nedges; j++ )
                    cout << " " << BaseAlpha( edges[j] );
               cout << "\n";
               cout << "critical edges =";
               for ( int j = 0; j < critical.isize( ); j++ )
                    cout << " " << BaseAlpha( critical[j] );
               cout << "\n";
               PRINT( int(good ) );    }    }

          // Recover any long edges that were not used.

          vec<Bool> used( nedges, False );
          for ( int j = 0; j < s.Ntigs( ); j++ )
               used[ Position( edges, s.Tig(j) ) ] = True;
          for ( int ei = 0; ei < nedges; ei++ )
          {    if ( !used[ei] )
               {    fastavector e = H->EdgeObject( edges[ei] );
	            if ( (int) e.size( ) >= MIN_EDGE_TO_SAVE_GLOBAL )
                    {    superb s;
                         s.SetNtigs(1);
	                 s.SetTig( 0, edges[ei] );
	                 s.SetLen( 0, e.size( ) );    
                         ss.push_back(s);    }    }    }    }

     pthread_mutex_lock( &scaffold_component_lock );
     for ( size_t j = 0; j < ss.size(); j++ )
     {    S.superbs.push_back( ss[j] );
          // cout << "component " << i << " yields " << ss[j].FullLength( ) // XXXXX
          //      << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          S.superb_indices.push_back( make_pair(i,j) );    }
     pthread_mutex_unlock( &scaffold_component_lock );    }

class ScaffoldComponentProcessor {
     public:
     void operator( )( int i )
     {    ScaffoldComponent(i);    }
};

void *ScaffoldComponentProcessorX(void *threadid)
{    int *tid = (int *) threadid;
     int i = *tid;
     ScaffoldComponent(i);
     return 0;    }

Bool cmp_size_reverse( const vec<int>& v1, const vec<int>& v2 )
{    if ( v1.size( ) > v2.size( ) ) return True;
     if ( v1.size( ) < v2.size( ) ) return False;
     return v1 > v2;    }

void HyperFastavector::ConvertToSuperbs( vec<superb> & superbs,
     const int max_threads, const int MIN_EDGE_TO_SAVE ) const
{
     MIN_EDGE_TO_SAVE_GLOBAL = MIN_EDGE_TO_SAVE;

     // Find components.

     S.H = this;
     ComponentsE(S.comp);
     sort( S.comp.begin( ), S.comp.end( ), cmp_size_reverse );
     ToLeft(S.to_left), ToRight(S.to_right);

     // Scaffold them.

     /*
     ScaffoldComponentProcessor p;
     Worklist<int,ScaffoldComponentProcessor> worklist( p, max_threads );
     for ( size_t i = 0; i < S.comp.size(); i++ )
          worklist.add(i);
     worklist.wait( ); 
     */

     if ( TEMP_DEBUG ) cout << Date() << ": A\t" << S.comp.size() << endl;
     for ( size_t l = 0; l < S.comp.size(); l++ )
     {    
          ScaffoldComponent(l);
          continue;
          int JOBS = Min( max_threads, int(S.comp.size() - l) );
          vec<pthread_t> threads(JOBS);
          vec<unsigned int> ids(JOBS);
          for ( int i = 0; i < JOBS; i++ )
          {    ids[i] = l+i;
               pthread_create( &threads[i], NULL,
                    ScaffoldComponentProcessorX, (void*) &ids[i] );    }
          for ( int i = 0; i < JOBS; i++ )
            {
               pthread_join( threads[i], NULL );
            }
	  if ( JOBS > 0 )
	    l += JOBS - 1;
     }
     if ( TEMP_DEBUG ) cout << Date() << ": Z" << endl;

     SortSync( S.superb_indices, S.superbs );
     /*
     for ( size_t i = 0; i < S.superbs.size(); i++ ) // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          PRINT2( i, S.superbs[i].FullLength( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     */
     superbs = S.superbs;    }





// Write a fasta-format file which contains the bases in this HyperFastavector
// scaffolded together (with gaps, etc.) as defined by the scaffolds.
// If supplied, rc defines which contigs are reverse-complement.
void
HyperFastavector::WriteScaffoldsFasta( const String & fasta_out,
				       const vec<superb> & scaffolds,
				       const vec<Bool> & rc ) const
{
  WriteScaffoldedFasta( fasta_out, this->Edges( ), scaffolds, rc );
}


void HyperFastavector::PrintSummaryDOT0w( ostream& out, Bool label_contigs,
     Bool label_vertices, Bool label_edges, const vec<int>* componentsToPrint,
     const Bool edge_labels_base_alpha, const vec<String> *label_edges_extra,
     const vec<String> *label_contigs_extra, const vec<int> *verticesToPrint ) const
{    vec<double> lengths( EdgeObjectCount( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
         lengths[i] = EdgeObject(i).size( );
     PrettyDOT( out, lengths, label_contigs, label_vertices, label_edges,
         componentsToPrint, edge_labels_base_alpha,
	 label_edges_extra, label_contigs_extra, verticesToPrint );    }


void BinaryWrite( int fd, const HyperFastavector& h )
{    WriteBytes( fd, &h.K_, sizeof(int) );
     BinaryWrite( fd, (const digraphE<fastavector>&) h );    }

void BinaryRead( int fd, HyperFastavector& h )
{    ReadBytes( fd, &h.K_, sizeof(int) );
     BinaryRead( fd, (digraphE<fastavector>&) h );    }

void HyperFastavector::LowerK( int newK )
{    ForceAssertLe( newK, K( ) );
     for ( int i = 0; i < EdgeObjectCount( ); i++ )
     {    ForceAssertGt( EdgeObject(i).size( ), 0u );
          EdgeObjectMutable(i).resize( EdgeObject(i).size( ) + newK - K( ) );    }
     K_ = newK;    }

Bool operator==( const HyperFastavector& h1, const HyperFastavector& h2 )
{    if ( h1.K( ) != h2.K( ) ) return False;
     return (const digraphE<fastavector>&) h1 
          == (const digraphE<fastavector>&) h2;    }

HyperFastavector::HyperFastavector( const HyperBasevector& h )
{    vec<fastavector> fastas( h.EdgeObjectCount( ) );
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
     {    fastas[i].SetFrom( h.EdgeObject(i) );    }
     (*this) = HyperFastavector( h.K( ), h, fastas );    }

HyperFastavector::HyperFastavector( const String& filename )
{    int fd = OpenForRead(filename);
     BinaryRead( fd, *this );
     close(fd);    }

void HyperFastavector::RemoveUnneededVertices( )
{    for ( int i = 0; i < N( ); i++ )
     {    if ( From(i).size( ) == 1 && To(i).size( ) == 1 && From(i)[0] != i )
          {    fastavector p = EdgeObjectByIndexTo( i, 0 );
               p.resize( p.size( ) - K( ) + 1 );
               p.Append( EdgeObjectByIndexFrom( i, 0 ) );
               JoinEdges( i, p );    }    }
     RemoveEdgelessVertices( );    }

String EdgeLabel( int v, int w, const fastavector& p, int K )
{    ostringstream out;
     out << v << " --- " << setprecision(9) << p.size( ) << " --> " << w;
     return out.str( );    }

void HyperFastavector::PrintSummaryPlus( ostream& out, const Bool show_bases ) const
{    
     // Define equivalence relation whose orbits are the components.

     equiv_rel e( N( ) );
     for ( int v = 0; v < N( ); v++ )
     {    for ( size_t j = 0; j < From(v).size(); j++ )
               e.Join( v, From(v)[j] );    }

     // Go through the components.  

     int count = 0;
     for ( int x = 0; x < N( ); x++ )
     {    if ( e.Representative(x) )
          {    out << "\ncomponent " << count << "\n";
               ++count;
               vec<int> o;
               e.Orbit( x, o );
               Sort(o);
               vec<float> pos( o.size( ) );
               vec<Bool> placed( o.size( ), False );
               pos[0] = 0.0, placed[0] = True;
               while( Sum(placed) < (int)o.size() )
               {    for ( size_t i1 = 0; i1 < o.size(); i1++ )
                    {    int v = o[i1];
                         for ( size_t j = 0; j < From(v).size(); j++ )
                         {    int w = From(v)[j];
                              int i2 = BinPosition( o, w );
                              if ( !( placed[i1] ^ placed[i2] ) ) continue;
                              const fastavector& p = EdgeObjectByIndexFrom( v, j );
                              if ( placed[i1] ) pos[i2] = pos[i1] + p.size( );
                              else pos[i1] = pos[i2] - p.size( );    
                              placed[i1] = placed[i2] = True;    }    }    }
               vec<float> opos( o.size( ) );
               for ( size_t i = 0; i < o.size(); i++ )
                    opos[i] = pos[i];
               SortSync( opos, o );
               int edgeid = 0;
               for ( size_t i = 0; i < o.size(); i++ )
               {    int v = o[i];
                    vec<int> f = From(v);
                    vec<int> ind( f.size( ) );
                    for ( size_t j = 0; j < ind.size(); j++ )
                         ind[j] = j;
                    vec<float> flen( f.size( ) );
                    for ( size_t j = 0; j < f.size(); j++ )
                         flen[j] = EdgeObjectByIndexFrom( v, j ).size( );
                    SortSync( flen, ind );
                    for ( size_t j = 0; j < ind.size(); j++ )
                    {    ++edgeid;
                         int w = f[ ind[j] ];
                         const fastavector& p = EdgeObjectByIndexFrom( v, ind[j] );

                         // Print edge.

                         int e = EdgeObjectIndexByIndexFrom( v, ind[j] );
                         out << "\n" << "[";
                         out << count-1 << "." << edgeid-1 
                             << " = " << BaseAlpha(e) << "] ";
                         out << EdgeLabel( v, w, p, K( ) ) << "\n";
                         if (show_bases) p.Print( out, ToString(e) );
                         if ( To(v).empty( ) ) out << "[" << v << " is source]\n";
                         if ( From(w).empty( ) ) out << "[" << w << " is sink]\n";
                                   }    }    }    }
     flush(out);    }
