///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LINKING_PAIRS_H
#define LINKING_PAIRS_H

#include "CoreTools.h"
#include "MainTools.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "PairsManager.h"

/*
 * A centralized container to store all the links between a set of contigs or
 * unibases.  The links can be added one-by-one, and retrieved together by
 * libraries. 
 *
 * The class identify different contigs or unibase by indexes, hence logical
 * relation between the contigs  shall be maintained by other externally. For
 * example, by using different index, * fw and rc of the same contig are
 * treated as different one. 
 */

// -------------- Definitions------------------------------

//// A class to store all links between two contigs
//// LibLinks[libid][link_id] is a pair of ( dist1, dist2 ), which are
typedef vec< vec< pair<int, int> > > LibLink_t;
// LibCov_t[libid] is a set of triples (pos, partner_tig, partner_pos) to track the
// location of a pair of read.
typedef triple<int,int,int> LinkEnd_t;          // ( pos, partner_tig, partner_pos )
typedef vec< multiset< LinkEnd_t > > LibCov_t;

class LinkingPairs
{
public:

  // -------------- Constructors ------------------------------

  LinkingPairs () { };
  LinkingPairs (int nlibs_, const vec<int>& lens) { Init( nlibs_, lens); }
  void Init(int nlibs, const vec<int>& lens);

  // -------------- Retrieving methods ------------------------------

  int                  NLibs ( ) const { return nlibs; }
  int                  NTigs ( ) const { return ntigs; }
  const vec<int>&      GetLens( )  const { return tlens; }
  int                  GetLen ( int i ) const { return tlens[i]; }

  const LibCov_t&      GetStarts ( int itig ) const { return starts[itig]; }
  const LibCov_t&      GetStops ( int itig ) const { return stops[itig]; }

  void                 GetStartsSmoothed( int u1, vec< map<int,double> >& x1s, int delta=50 ) const;
  void                 GetStopsSmoothed( int u1, vec< map<int,double> >& x1s, int delta=50 ) const;
  void                 GetStartsSampled( int u1, vec< map<int,double> >& x1s ) const;
  void                 GetStopsSampled( int u1, vec< map<int,double> >& x1s ) const;

  const map<int,int>&  GetDists ( int ilib ) const { return dists[ilib]; }
  const vec<bool>&     GetTigPick( int ilib) const { return tig_pick[ilib]; }
  const LibLink_t&     GetLinks ( pair<int,int> u1u2 ) const {
                           map< pair<int,int>, LibLink_t >::const_iterator it = all_links.find( u1u2 );
                           if ( it != all_links.end() ) return it->second;
                           else return _empty;
  }
  const LibLink_t&     GetLinks ( int u1, int u2 ) const { return GetLinks( make_pair(u1,u2) ); }
  vec<pair<int,int> >  GetAllLinked( ) const { 
                           vec< pair<int,int> > results;
                           for( map< pair<int,int>, LibLink_t >::const_iterator it = all_links.begin(); 
                                     it != all_links.end(); it++ )
                                results.push_back( it->first );
                           return results;
  }

  // -------------- Modification methods ------------------------------

  // Add a link that starts from itig at position pos, and stops at partner_tig at posotion partner_pos 
  void AddStart( int ilib, int itig, const triple<int, int, int>& link  ) { starts[itig][ilib].insert(link); }
  void AddStart( int ilib, int itig, int pos, int partner_tig, int partner_pos ) 
                                    { AddStart(ilib, itig, triple<int,int,int>( pos, partner_tig, partner_pos ) ); }
  // Add a link that stops at  itig at position pos, and starts from partner_tig at position partner_pos 
  void AddStop( int ilib, int itig, const triple<int, int, int>& link  ) { stops[itig][ilib].insert(link); }
  void AddStop( int ilib, int itig, int pos, int partner_tig, int partner_pos ) 
                                    { AddStop(ilib, itig, triple<int,int,int>( pos, partner_tig, partner_pos ) ); }
  // Add a linking pair that starts at x1@u1 and stops at x2@u2. 
  void AddPair( int ilib, int u1, int x1, int u2, int x2 ) ; 
  // Add a pair separation 
  void AddDist ( int libid, int d, int u1, int u2 ) { dists[libid][d]++; tig_pick[libid][u1] = tig_pick[libid][u2] = true; }

  // -------------- Diagnostic methods ------------------------------
  void DumpInfo ( ostream& out, int u1, int u2, int delta = 50 ) const;
  int GetNStarts( int libid, int u ) const ; 
  int GetNStops ( int libid, int u ) const ; 

private:

  int nlibs, ntigs;
  vec <int> tlens;                    // lenghts of the contigs
  map< pair<int,int>, LibLink_t > all_links; // all_links[ make_pair(uid1, uid2) ][iLib][i] is a link (x1,x2)
                                             // that starts from pos x1 at u1 and end at pos x2 at u2.
  vec< LibCov_t > starts;            // starts[itig][ilib] is a set of links that starts from this contig
  vec< LibCov_t > stops;             // stops [itig][ilib] is a set of links that  stops from this contig
  // pair seperation statistics
  vec< map< int, int > > dists;      // dists[ilib][x] is the number of pairs at sep x
  vec< vec< bool > > tig_pick;      // unibase_pick[ilib][itig] mark if a contig is used for separation distribution

  LibLink_t _empty;                          // an empty array to be returned if request non-existing linking
};


void DumpLinkingInfo( ostream& out, 
   const vec< vec<double> > & x1s, const vec< vec<double> > & x2s,
   const vec< vec< pair<int,int> > > & links, int step = 10, bool flip = true, int len1 = 0 ) ;

void DumpLinkingInfo( ostream& out, 
   const vec< map<int,double> > & x1s, const vec< map<int,double> > & x2s,
   const vec< vec< pair<int,int> > > & links, int step = 10, bool flip = true, int len1 = 0 ) ;

void GatherLinks( const PairsManager& pairs, const vec< pair<int,int> > aligns, 
    const vec<int>& to_rc,  // end of input
     LinkingPairs &linking );

#endif
