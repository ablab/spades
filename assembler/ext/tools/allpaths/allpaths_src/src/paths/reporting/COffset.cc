///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "VecTemplate.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/reporting/COffset.h"

/**
 * COffset
 * Constructor
 */
COffset::COffset( ) :
  super1_ ( 0 ),
  super2_ ( 0 ),
  slen1_ ( 0 ),
  slen2_ ( 0 ) { }
  
/**
 * COffset
 * Constructor
 */
COffset::COffset( int s1, int s2, bool rc2, int slen1, int slen2 )
{
  this->SetSupers( s1, s2, rc2, slen1, slen2 );
}
  
/**
 * COffset
 * SetSupers
 */
void COffset::SetSupers( int s1, int s2, bool rc2, int slen1, int slen2 )
{
  super1_ = s1;
  super2_ = rc2 ? - s2 - 1 : s2;
  slen1_ = slen1;
  slen2_ = slen2;
}

/**
 * COffset
 * AddLink
 *
 * After all links have been added, run ClusterLinks to cluster them
 * into consistent sets. No new links can be added after clustering
 * (hence the assert if links_.size( ) > 1).
 */
void COffset::AddLink( const SLink &link )
{
  if ( links_.size( ) < 1 ) links_.resize( 1 );
  ForceAssert( links_.size( ) == 1 );
  links_[0].push_back( link );
}

/**
 * COffset
 * NClusters
 */
size_t COffset::NClusters( ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  return links_.size( );
}

/**
 * COffset
 * NLinksTotal
 */
size_t COffset::NLinksTotal( ) const
{
  // No need to run ClusterLinks.
  size_t total = 0;
  for (size_t ii=0; ii<links_.size( ); ii++)
    total += links_[ii].size( );
  return total;
}

/**
 * COffset
 * NLinks
 */
size_t COffset::NLinks( int cluster_id ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  return links_[cluster_id].size( );
}

/**
 * COffset
 * Score
 */
float COffset::Score( int cluster_id ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  return scores_[cluster_id];
}

/**
 * COffset
 * CombinedScore
 *
 * It mimics the CombinedScore in CLinkBundle: it returns an overall
 * combined score of Score( ) and NLinks( ). The returned value is in
 * the range (0, sqrt( NLinks( ) ], where higher is better.
 *
 * HEURISTICS: use CENTER to decide where the real distribution is
 * centered (in theory it should be 1.0).
 *
 * HEURISTICS: use DAMPING to damp the effect of the exponential
 * (at 1.0 it would decay too fast).
 */
double COffset::CombinedScore( int cluster_id ) const
{
  const double CENTER = 1.0;
  const double DAMPING = 100.0;
  double N = (double)this->NLinks( cluster_id );
  double score = (double)this->Score( cluster_id );
  double sig = ( score - CENTER ) * ( score - CENTER );
  return sqrt( N ) * exp( - N * sig / DAMPING );
}

/**
 * COffset
 * Offset
 */
normal_distribution COffset::Offset( int cluster_id ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  return offsets_[cluster_id];
}

/**
 * COffset
 * Offset
 *
 * This version of Offset breaks the bundle into libraries:
 * lib_nds[lib_id] is the combined normal distribution of all the
 * links in the bundle from library lib_id, or (-1,-1) if there are no
 * links from lib_id in bundle.
 */
normal_distribution COffset::Offset( int cluster_id,
				     const PairsManager &pairs,
				     vec<normal_distribution> &lib_nds ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );

  lib_nds.clear( );
  lib_nds.resize( pairs.nLibraries( ), normal_distribution( -1.0, -1.0 ) );
  
  vec< vec<normal_distribution> > lib_offsets( pairs.nLibraries( ) );
  for (int ii=0; ii<links_[cluster_id].isize( ); ii++) {
    const SLink &link = links_[cluster_id][ii];
    const int lib_id = pairs.libraryID( link.pair_id_ );
    lib_offsets[lib_id].push_back( link.offset_ );
  }

  for (int ii=0; ii<(int)pairs.nLibraries( ); ii++) {
    if ( lib_offsets.size( ) < 1 ) continue;
    lib_nds[ii] = CombineNormalDistributions( lib_offsets[ii] );
  }
  
  return offsets_[cluster_id];
}

/**
 * COffset
 * DevErrorEstimate
 *
 * It estimates the error in the measurement of the mean of the gap.
 * The error for a library is defined as lib_dev * RATIO, and the
 * error for a cluster as the min of the errors over all the libraries
 * involved in the cluster.
 *
 * HEURISTICS: use RATIO to define the error rate of a library.
 */
int COffset::DevErrorEstimate( int cluster_id,
			       const PairsManager &pairs ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  
  const double RATIO = 0.03;

  vec<int> involved( pairs.nLibraries( ), false );
  for (size_t ii=0; ii<links_[cluster_id].size( ); ii++)
    involved[ pairs.libraryID( links_[cluster_id][ii].pair_id_ ) ] += 1;

  int select = -1;
  int selected_size = 0;
  for (size_t ii=0; ii<involved.size( ); ii++) {
    if ( involved[ii] > selected_size ) {
      selected_size = involved[ii];
      select = ii;
    }
  }
  ForceAssertGt( select, -1 );
  
  return Max( 1, int( double( pairs.getLibrarySD( select ) ) * RATIO ) );
}

/**
 * COffset
 * SpreadWin1
 */
pair<double,double> COffset::SpreadWin1( int cluster_id ) const
{
  pair<int,int> win = this->SpreadWinBases1( cluster_id );
  double minr = 100.0 * SafeQuotient( win.first, slen1_ );
  double maxr = 100.0 * SafeQuotient( win.second, slen1_ );

  return make_pair( minr, maxr );
}

/**
 * COffset
 * SpreadWin2
 */
pair<double,double> COffset::SpreadWin2( int cluster_id ) const
{
  pair<int,int> win = this->SpreadWinBases2( cluster_id );
  double minr = 100.0 * SafeQuotient( win.first, slen2_ );
  double maxr = 100.0 * SafeQuotient( win.second, slen2_ );

  return make_pair( minr, maxr );
}

/**
 * COffset
 * SpreadWinBases1
 */
pair<int,int> COffset::SpreadWinBases1( int cluster_id ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  
  const vec<SLink> &clust = links_[cluster_id];
  vec<int> lens;
  lens.reserve( clust.size( ) );
  for (size_t ii=0; ii<clust.size( ); ii++)
    lens.push_back( clust[ii].win1_.Start( ) );
  
  return this->CoreSpread( lens );
}

/**
 * COffset
 * SpreadWinBases2
 */
pair<int,int> COffset::SpreadWinBases2( int cluster_id ) const
{
  if ( offsets_.size( ) < 1 ) this->ClusterLinks( );
  
  const vec<SLink> &clust = links_[cluster_id];
  vec<int> lens;
  lens.reserve( clust.size( ) );
  for (size_t ii=0; ii<clust.size( ); ii++)
    lens.push_back( clust[ii].win2_.Start( ) );
  
  return this->CoreSpread( lens );
}

/**
 * COffset
 * Separation
 */
normal_distribution COffset::Separation( int cluster_id ) const
{
  normal_distribution nd = this->Offset( cluster_id );
  if ( nd.mu_ > 0 ) nd.mu_ = nd.mu_ - slen1_;
  else nd.mu_ = - nd.mu_ - slen2_;
  return nd;
}

/**
 * COffset
 * MatchesSupersWith
 */
bool COffset::MatchesSupersWith( const COffset &other ) const
{
  if ( this->Super1( ) != other.Super1( ) ) return false;
  if ( this->Super2( ) != other.Super2( ) ) return false;
  return ( this->Rc2( ) == other.Rc2( ) );
}

/**
 * COffset
 * IsValid
 *
 * It runs validity tests on the given cluster. In this order:
 *
 * if MIN_LINKS: accept only clusters with these many links or more
 * if MAX_SPREAD: discard a cluster if either of its spreads is > MAX_SPREAD
 * if MIN_MAX_GAP: accept a cluster only if implied sep is in the given window
 */
bool COffset::IsValid( int cluster_id,
		       const size_t *MIN_LINKS,
		       const int *MAX_SPREAD,
		       const pair<int,int> *MIN_MAX_GAP ) const
{
  if ( MIN_LINKS ) {
    if ( this->NLinks( cluster_id ) < *MIN_LINKS )
      return false;
  }

  if ( MAX_SPREAD ) {
    pair<int,int> win1 = this->SpreadWinBases1( cluster_id );
    pair<int,int> win2 = this->SpreadWinBases2( cluster_id );
    int spread1 = win1.second - win1.first;
    int spread2 = win2.second - win2.first;
    if ( Min( spread1, spread2 ) > *MAX_SPREAD )
      return false;
  }

  if ( MIN_MAX_GAP ) {
    normal_distribution nd = this->Separation( cluster_id );
    if ( nd.mu_ < MIN_MAX_GAP->first || nd.mu_ > MIN_MAX_GAP->second )
      return false;
  }
      
  return true;  
}

/**
 * COffset
 * Print
 *
 * It will assert if ClusterLinks was not run. If pairs are given, all
 * links will be printed (many lines). If smart is true, run IsValid
 * above with a "reasonable" selection of arguments.
 */
void COffset::Print( ostream &out,
		     const PairsManager *pairs,
		     const bool smart ) const
{
  // Only used if smart = true.
  const size_t MIN_LINKS = 3;
  const int MAX_SPREAD = 20000;
  const pair<int,int> MIN_MAX_GAP = make_pair( -6000, 25000 );

  if ( this->NLinksTotal( ) < 1 ) return;
  
  ForceAssert( offsets_.size( ) > 0 );
  
  String str_or = this->Rc2( ) ? "[-]" : "[+]";
  out << "s" << super1_ << " (" << slen1_ << " bp)"
      << "  --->  s" << this->Super2( ) << str_or << " (" << slen2_ << " bp)"
      << " " << this->NLinksTotal( ) << " link(s)\n";

  int n_links = links_.size( );
  String str_sup = this->ToStringSupers( );
  
  for (size_t cluster_id=0; cluster_id<links_.size( ); cluster_id++) {
    if ( smart ) {
      if ( ! IsValid( cluster_id, &MIN_LINKS, &MAX_SPREAD, &MIN_MAX_GAP ) )
	continue;
    }
    
    out << "  " << this->ToStringCluster( cluster_id ) << "\n";

    if ( ! pairs ) continue;
    
    vec< vec<String> > table;
    vec<String> tline;
  
    String str_clu = "cl_" + ToString( cluster_id ) + "/" + ToString( n_links );
    int n_clinks = links_[cluster_id].size( );
    for (size_t ii=0; ii<links_[cluster_id].size( ); ii++) {
      tline.clear( );

      pair<int,int> valid1 = this->SpreadWinBases1( cluster_id );
      pair<int,int> valid2 = this->SpreadWinBases2( cluster_id );
      int pair_id = links_[cluster_id][ii].pair_id_;
      int lib_id = pairs->libraryID( pair_id );
      int lib_sep = pairs->getLibrarySep( lib_id );
      int lib_sd = pairs->getLibrarySD( lib_id );
      int start1 = links_[cluster_id][ii].win1_.Start( );
      int start2 = links_[cluster_id][ii].win2_.Start( );
      bool outlier1 = start1 < valid1.first || start1 >= valid1.second;
      bool outlier2 = start2 < valid2.first || start2 >= valid2.second;
      String str_id = ToString( pair_id ) + " (" + ToString( lib_id );
      String str_name = pairs->getLibraryName( lib_id );
      String str_sep = ToString( lib_sep ) + "+/-" + ToString( lib_sd );
      String str_pos1 = ( outlier1 ? "OUTLIER @" : "@" );
      String str_pos2 = ( outlier2 ? "OUTLIER @" : "@" );

      tline.push_back( "  " );
      tline.push_back( str_sup );
      tline.push_back( str_clu );
      tline.push_back( "link_" + ToString( ii ) + "." + ToString( n_clinks ) );
      tline.push_back( this->ToStringND( links_[cluster_id][ii].offset_ ) );
      tline.push_back( str_pos1 + ToString( start1 ) );
      tline.push_back( str_pos2 + ToString( start2 ) );
      tline.push_back( str_id + " = " + str_name + " " + str_sep + ")" );
      
      table.push_back( tline );
    }

    PrintTabular( out, table, 2, "rrrrrrrl" );
  }
  out << endl;
  
}
  
/**
 * COffset
 * ClusterLinks
 *
 * Sort all links, and then define two consecutive links to be
 * consistent with each other iff the two intervals with center the
 * offset and radius ( stretch * links_stdev ) overlap. ClusterLinks
 * will also generate offsets_, bases_, and scores_.
 */
void COffset::ClusterLinks( const float stretch ) const
{
  if ( links_.size( ) < 1 ) return;
  ForceAssert( links_.size( ) == 1 );
  vec<SLink> all_links = links_[0];
  links_.clear( );

  sort( all_links.begin( ), all_links.end( ) );
  links_.push_back( vec<SLink>( 1, all_links[0] ) );
  for (size_t ii=1; ii< all_links.size( ); ii++) {
    if ( ! all_links[ii-1].OverlapsWith( all_links[ii], stretch ) )
      links_.resize( links_.size( ) + 1 );
    links_[links_.size( )-1].push_back( all_links[ii] );
  }
  
  offsets_.resize( links_.size( ) );
  scores_.resize( links_.size( ), -1.0 );
  for (size_t ii=0; ii<links_.size( ); ii++) {
    vec<normal_distribution> nds;
    nds.reserve( links_[ii].size( ) );
    for (size_t jj=0; jj<links_[ii].size( ); jj++)
      nds.push_back( links_[ii][jj].offset_ );
    offsets_[ii] = CombineNormalDistributions( nds, &( scores_[ii] ) );
  }
  
}

/**
 * COffset
 * CoreSpread
 * private
 *
 * Look at the start positions of reads on a super, and determine if
 * there are any outliers, then determine the spread.  It returns
 * (-1,-1) on error.
 *
 * WARNING! This is not an optimal solution for two reasons:
 *   1. outliers should be detected by looking at the spreads on both
 *      supers, not just one super at a time.
 *   2. offsets_, scores_, and the NLinks( ) method (which provide the
 *      weight for a cluster) should reflect the fact that the cluster
 *      may be using less links that the whole set.
 */
pair<int,int> COffset::CoreSpread( const vec<int> &starts ) const
{
  // This could be made into an argument.
  const double MAX_STRETCH = 3.5;

  // Nothing to do, return error code (-1,-1).
  if ( starts.size( ) < 1 ) return make_pair( -1, -1 );

  // Remove outliers and find spread.
  int win_min = -1;
  int win_max = -1;
  double starts_mean = Mean( starts );
  double starts_dev = Max( 1.0, StdDev( starts, starts_mean ) );
  int range_low = int( starts_mean - ( starts_dev * MAX_STRETCH ) );
  int range_high = int( starts_mean + ( starts_dev * MAX_STRETCH ) );
  for (size_t ii=0; ii<starts.size( ); ii++) {
    if ( starts[ii] < range_low || starts[ii] > range_high ) continue;
    if ( win_min < 0 && win_max < 0 ) {
      win_min = starts[ii];
      win_max = starts[ii];
      continue;
    }
    win_min = Min( win_min, starts[ii] );
    win_max = Max( win_max, starts[ii] );
  }
  
  return make_pair( win_min, 1 + win_max );
}

/**
 * COffset
 * ToStringSupers
 * private
 */
String COffset::ToStringSupers( ) const
{
  return 
    "s" + ToString( super1_ )
    + " s" + ToString( this->Super2( ) )
    + "[" + ( this->Rc2( ) ? "-" : "+" ) + "]";
}

/**
 * COffset
 * ToStringCluster
 * private
 */
String COffset::ToStringCluster( const int cluster_id ) const  
{
  pair<int,int> win1 = this->SpreadWinBases1( cluster_id );
  pair<int,int> win2 = this->SpreadWinBases2( cluster_id );
  int spread1 = win1.second - win1.first;
  int spread2 = win2.second - win2.first;
  String str1
    = ToString( spread1 ) + " @[" + ToString( win1.first ) + 
    "," + ToString( win1.second ) + ")_" + ToString( slen1_ );
  String str2
    = ToString( spread2 ) + " @[" + ToString( win2.first ) + 
    "," + ToString( win2.second ) + ")_" + ToString( slen2_ );
  
  return 
    "cl=" + ToString( cluster_id ) + "." + ToString( links_.size( ) )
    + " w=" + ToString( links_[cluster_id].size( ) )
    + " score=" + ToString( scores_[cluster_id], 2 )
    + " spread1=" + str1
    + " spread2=" + str2
    + " " + this->ToStringND( offsets_[cluster_id] );
}

/**
 * COffset
 * ToStringND
 * private
 */
String COffset::ToStringND( const normal_distribution &nd ) const
{
  return "<" + ToString( nd.mu_, 0 ) + "," + ToString( nd.sigma_, 0 ) + ">";
}

