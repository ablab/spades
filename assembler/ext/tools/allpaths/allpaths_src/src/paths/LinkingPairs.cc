///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "paths/LinkingPairs.h"
#include "math/Array.h"

void LinkingPairs::Init(int nlibs_, const vec<int>& lens )
{
  nlibs = nlibs_;
  ntigs = lens.size();
  tlens = lens;
  starts.assign( ntigs, LibCov_t( nlibs) );
  stops.assign( ntigs, LibCov_t( nlibs ) );
  dists.assign( nlibs, map<int,int>() );
  tig_pick.assign( nlibs, vec<bool> ( ntigs, false ) );
  all_links.clear();
  cout << "nlibs= " << nlibs << endl;
  cout << "ntigs= " << ntigs << endl;
  _empty.assign( nlibs, vec< pair<int,int> >() );
}


// Add a link that starts from pos x1 at contig u1 and stops at pos x2 at contig u2
void LinkingPairs::AddPair( int ilib, int u1, int start1, int u2, int stop2 ) { 
  ForceAssert( ilib >= 0 && ilib < nlibs );
  bool aligned1 = ( u1 >= 0 && u1 < ntigs );
  bool aligned2 = ( u2 >= 0 && u2 < ntigs );
  if ( aligned1 && aligned2 && u1 != u2 ) {
    all_links[ make_pair(u1,u2) ].resize(nlibs); 
    all_links[ make_pair(u1,u2) ][ilib].push( start1, stop2 );
  }
  if ( aligned1 ) 
    AddStart( ilib, u1, start1, u2, stop2 );
  if ( aligned2 ) 
    AddStop ( ilib, u2, stop2, u1, start1 );
}


void LinkingPairs::GetStartsSmoothed( int u1, vec< map<int,double> >& x1s, int delta ) const {
  const LibCov_t& x1 = this->GetStarts( u1);
  for ( int ilib = 0; ilib < nlibs; ilib++ ) {
    for ( multiset<LinkEnd_t> ::const_iterator it = x1[ilib].begin(); it != x1[ilib].end(); it++ )
      x1s[ilib][it->first] += 1;
    SmoothArrayGaussian( x1s[ilib], delta );
  }
}

void LinkingPairs::GetStopsSmoothed( int u1, vec< map<int,double> >& x1s, int delta ) const {
  const LibCov_t& x1 = this->GetStops( u1);
  for ( int ilib = 0; ilib < nlibs; ilib++ ) {
    for ( multiset<LinkEnd_t> ::const_iterator it = x1[ilib].begin(); it != x1[ilib].end(); it++ )
      x1s[ilib][it->first] += 1;
    SmoothArrayGaussian( x1s[ilib], delta );
  }
}

void LinkingPairs::GetStartsSampled( int u1, vec< map<int,double> >& x1s ) const {
  const LibCov_t& x1 = this->GetStarts( u1);
  const int len = tlens[u1];
  for ( int ilib = 0; ilib < nlibs; ilib++ ) {
    x1s[ilib].clear();
    for ( multiset<LinkEnd_t> ::const_iterator it = x1[ilib].begin(); it != x1[ilib].end(); it++ )
      x1s[ilib][it->first] += 1;
    for ( map<int,double>::iterator it = x1s[ilib].begin(); it != x1s[ilib].end(); it++ ) {
      double count = it->second;
      // find the prev and next position
      int prev = - it->first, next = len + it->first;
      if ( it != x1s[ilib].begin() ) {
	prev = (--it)->first;
	++it;
      }
      ++it;
      if ( it != x1s[ilib].end() ) {
	next = it->first;
      }
      --it;
      x1s[ilib][it->first] = count * 2.0 / ( next - prev );
    }
  }
}

void LinkingPairs::GetStopsSampled( int u1, vec< map<int,double> >& x1s ) const {
  const LibCov_t& x1 = this->GetStops( u1);
  const int len = tlens[u1];
  for ( int ilib = 0; ilib < nlibs; ilib++ ) {
    x1s[ilib].clear();
    for ( multiset<LinkEnd_t> ::const_iterator it = x1[ilib].begin(); it != x1[ilib].end(); it++ )
      x1s[ilib][it->first] += 1;
    for ( map<int,double>::iterator it = x1s[ilib].begin(); it != x1s[ilib].end(); it++ ) {
      double count = it->second;
      // find the prev and next position
      int prev = - it->first, next = len + it->first;
      if ( it != x1s[ilib].begin() ) {
	prev = (--it)->first;
	++it;
      }
      ++it;
      if ( it != x1s[ilib].end() ) {
	next = (++it)->first;
      }
      --it;
      x1s[ilib][it->first] = count * 2.0 / ( next - prev );
    }
  }
}

void LinkingPairs::DumpInfo( ostream& out, int u1, int u2, int delta ) const
{
  vec< map<int, double > > x1s( nlibs );        
  vec< map<int, double > > x2s( nlibs );       
  for ( int ilib = 0; ilib < nlibs; ilib++ ) {
       GetStartsSmoothed( u1, x1s, delta );
       GetStopsSmoothed( u1, x1s, delta );
       SmoothArrayGaussian( x1s[ilib], delta );
       SmoothArrayGaussian( x2s[ilib], delta );
  }
  DumpLinkingInfo( out, x1s, x2s, this->GetLinks(u1,u2) );
}

int LinkingPairs::GetNStarts ( int libid, int u ) const {
     ForceAssertLt( libid, nlibs );
     return starts[u][libid].size();
}

int LinkingPairs::GetNStops ( int libid, int u ) const {
     ForceAssertLt( libid, nlibs );
     return stops[u][libid].size();
}

void DumpLinkingInfo( ostream& out, 
   const vec< map<int,double> > & x1s, const vec< map<int,double> > & x2s,
   const vec< vec< pair<int,int> > > & links, int step, bool flip, int len1 ) 
{
  int nlibs = links.size();
  ForceAssertEq( nlibs, x1s.isize() );
  ForceAssertEq( nlibs, x2s.isize() );
  // output diagnostic information
  for ( int libid = 0; libid < nlibs; libid++ )
    for ( int i = 0; i < links[libid].isize(); i++ )
	out <<"link "<< libid  << " " << links[libid][i].first - (flip ? len1 : 0 ) 
	  << " "  <<  links[libid][i].second << endl;
  // coverage data
  for( int libId=0; libId < nlibs; libId++) {
    for ( map<int,double>::const_iterator it = x1s[libId].begin(); it!= x1s[libId].end(); it++ )
      if ( it->first % step == 0 )
	out << "counter1_"<<libId<<" " << it->first - (flip ? len1 : 0 )  << " " << it->second << endl;
    for ( map<int,double>::const_iterator it = x2s[libId].begin(); it!= x2s[libId].end(); it++ )
      if ( it->first % step == 0 )
        out << "counter2_"<<libId<<" " << it->first  << " " << it->second << endl;
  }
}

void DumpLinkingInfo( ostream& out, 
   const vec< vec<double> > & x1s, const vec< vec<double> > & x2s,
   const vec< vec< pair<int,int> > > & links, int step, bool flip, int len1 ) 
{
  int nlibs = links.size();
  ForceAssertEq( nlibs, x1s.isize() );
  ForceAssertEq( nlibs, x2s.isize() );
  // output diagnostic information
  for ( int libid = 0; libid < nlibs; libid++ )
    for ( int i = 0; i < links[libid].isize(); i++ )
      out <<"link "<< libid  << " " << links[libid][i].first - (flip ? len1 : 0 ) << " "  <<  links[libid][i].second << endl;
  // coverage data
  for( int libId=0; libId < nlibs; libId++) {
    for( int i = 0; i < x1s[libId].isize(); i++ )
      if ( i % step == 0 )
	out << "counter1_"<<libId<<" " << i - (flip ? len1 : 0 ) << " " << x1s[libId][i]<< endl;
    for( int i = 0; i < x2s[libId].isize(); i++ )
      if ( i % step == 0 )
	out << "counter1_"<<libId<<" " << i << " " << x2s[libId][i]<< endl;
  }
}

// Every read should have unique alignment.  
//
// The position of one end of the read (PER) is recorded ( usually the invariant, i.e., the
// sequencing starting end ). Of course you can choose any end, as long as the 
// choice is consistent when calculating the pair separation distributions.
// We make sure PER is in [0, tig_length ).
//
// Since RCs of the contigs are also included, each read is aligned to each contig in
// only one direction. I choose the direction so that the other reads are to be expected
// at position PER - s, where s is the distributuion of the pair spearation. In other
// words, only the 'stop' positions are included in the alignments. It
// corresponds to 'start' position of tig_len - PER in the RC of the contig.
//
// Diagram below illustrate the case of unfliped jump reads, where outies are jumps
// and innies are nonjumps. The symbol '@' marks the invarant read end position.
// Only the left cases are recorded.
//
// |<-x1->|<------y1------>                          |<------y1------>|<-x1->      
//        @----> r1                                              <----@ r1*
// ------------------------ tig1                     ------------------------ tig1*    
//                                      equiv to                                     
// |<-x2->|<------y2------>                          |<------y2------>|<-x2->      
//        @----> r2                                              <----@ r2*
// ------------------------ tig2                     ------------------------ tig2*    
//
// We record the starting and stopping positions on each contigs, counted from the
// beginning of the contig. The links are also counted in the same way. 
// 
void GatherLinks( const PairsManager& pairs, const vec< pair<int,int> > aligns, 
    const vec<int>& to_rc, 
     LinkingPairs &linking)
{
  int nlibs = pairs.nLibraries();
  int ntigs = to_rc.size();
  ForceAssertEq( ntigs, linking.NTigs() );
  ForceAssertEq( nlibs, linking.NLibs() );
  ForceAssertEq( pairs.nReads(), aligns.size() );

  //linking.Init( nlibs, ntigs ); // assuming already initialized to the right size

  const int cutoff_nonjump = 1000; // Maximium allowed nonjump distance
  const int cutoff_smalljump = 5;  // remove spurious peak of near-zero sized jumps
  // putative cutoffs 
  vec <int> cutoffs( nlibs, 0 );
  for ( int i = 0; i < nlibs; i++ ) 
    cutoffs[i] = pairs.getLibrarySep(i) + pairs.getLibrarySD(i) * 10;

  for ( size_t i = 0; i < pairs.nPairs( ); i++ ){    
    longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
    int libid = pairs.libraryID(i);
    int u1 = aligns[id1].first;
    int x1 = 0, y1 = 0;
    if ( u1 >= 0 ) {
      x1 = aligns[id1].second;  // the stop position on u1
      y1 = linking.GetLen( u1 )- x1 - 1;   // the start position on u1*
    }
    int u2 = aligns[id2].first;
    int x2 = 0, y2 = 0;
    if ( u2 >= 0 ) {
      x2 = aligns[id2].second;
      y2 = linking.GetLen( u2 ) - x2 - 1;
    }
    int ru1 = ( u1 >= 0 ? to_rc[u1] : -1 );  // the starting contig
    int ru2 = ( u2 >= 0 ? to_rc[u2] : -1 );
    linking.AddPair(libid, ru1, y1, u2, x2 );
    linking.AddPair(libid, ru2, y2, u1, x1 );
    if ( u1 >= 0 && u2 >= 0 && u1 == to_rc[u2] ) {
      // pairs aligned in the same contig
      // |<-x1->|<--------------- y1 ---------------------------->        
      //        @----> r1  
      // -------------------------------------------------------- tig1 
      // -------------------------------------------------------- tig2*  
      //                                             <-----@ r2*                                 
      // <--------------------- y2 ----------------------->|<-x2->      
      // The diagram above actually illustrates a nonjump (innie) case, where x1 < y2;                
      int d = x1 - y2; // remember that the pair starts from y2 and stops at x1
      if ( ( d > 0 && d < cutoffs[libid] && d > cutoff_smalljump )  // valid jumps
	  || ( d < 0 && -d < cutoff_nonjump && -d > cutoff_smalljump)  // valid nonjumps
	 ) 
	linking.AddDist( libid, d, u1, u2 );
    }
  }
}
