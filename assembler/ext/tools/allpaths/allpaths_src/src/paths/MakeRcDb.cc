/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "paths/KmerPath.h"
#include "paths/Unipath.h"

template<class TAG>
void BuildDbFile( vecKmerPath& paths,
		  vecKmerPath& paths_rc,
		  String pathsdb_file );

/**
 * MakeRcDb.cc
 *
 * Create standard paths_rc and pathsdb files.
 * All files are in <PRE>/<DATA>/<RUN>.
 *
 * Input:
 * <READS>.<PATHS>.k<K>
 *
 * Output:
 * <READS>.<PATHS>_rc.k<K>
 * <READS>.<PATHS>db.k<K>
 *
 * If NO_RC=True, skip the file <READS>.<PATHS>_rc.k<K>, and create a pathsdb
 * that is half as large.
 *
 * If FORCE_BIG=True, or if any of the paths contain KmerPathIntervals with
 * length greater than 65536, the pathsdb file will take the form of a
 * vec<big_tagged_rpint> instead of a vec<tagged_rpint>, and will be named
 * <READS>.<PATHS>db_big.k<K>.
 *
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault( PATHS, "paths" );
  CommandArgument_Bool_OrDefault( FORCE_BIG, False );
  CommandArgument_Bool_OrDefault( NO_RC, False );
  EndCommandArguments;

  String KS = ToString( K );
  String paths_head = PRE + "/" + DATA + "/" + RUN + "/" + READS + "." + PATHS;

  vecKmerPath paths;
  cout << Date() << ": Loading paths from file..." << endl;
  paths.ReadAll( paths_head + ".k" + KS );
  cout << "\tpaths:    Total " << paths.size() << " KmerPaths, containing " << paths.sumSizes() << " KmerPathIntervals." << endl;

  int max_segs=0;

  vecKmerPath paths_rc;
  if ( NO_RC )
    cout << Date() << ": [Not creating paths_rc because NO_RC=True.]" << endl;
  else {
    cout << Date() << ": Creating paths_rc" << endl;
    paths_rc.reserve( paths.size() );
    KmerPath temp;
    for(size_t i=0; i<paths.size(); i++) {
      temp = paths[i];
      max_segs = Max( max_segs, temp.NSegments() );
      temp.Reverse();
      max_segs = Max( max_segs, temp.NSegments() );
      paths_rc.push_back(temp);
    }
    
    String paths_rc_file = paths_head + "_rc.k" + KS;
    cout << "\tpaths_rc: Total " << paths_rc.size() << " KmerPaths, containing " << paths_rc.sumSizes() << " KmerPathIntervals." << endl;
    cout << Date() << ": Writing paths_rc to file" << endl;
    paths_rc.WriteAll( paths_rc_file );
  }

  bool need_big = uint(max_segs) >= tagged_rpint::POSITION_MAX;
  bool big = need_big || FORCE_BIG;
  String str_big = ( big ? "db_big.k" : "db.k") + KS;
  String pathsdb_file = paths_head + str_big;

  if ( need_big )
    cout << "   NOTE: the longest KmerPath has " << max_segs << " segments."
	 << "\n   This is larger than the " <<  tagged_rpint::POSITION_MAX
	 << " limit for a tagged_rpint,"
	 << "\n   so building a vec<big_tagged_rpint> instead" << endl;
  else if ( FORCE_BIG )
    cout << "FORCE_BIG=True, so building a vec<big_tagged_rpint>,"
	 << "\neven though there are no long paths which require it." << endl;

  if( big ) BuildDbFile<big_tagged_rpint>(paths, paths_rc, pathsdb_file );
  else BuildDbFile<tagged_rpint>(paths, paths_rc, pathsdb_file );

  cout << Date() << ": MakeRcDb is done!" << endl;
}

/**
 * BuildDbFile
 *
 * class TAG can be either a tagged_rpint or a big_tagged_rpint.
 */
template<class TAG>
void BuildDbFile( vecKmerPath& paths,
		  vecKmerPath& paths_rc,
		  String pathsdb_file )
{
  size_t paths_size_sum = paths.sumSizes();
  size_t paths_rc_size_sum = paths_rc.sumSizes();
  
  // Find total memory usage (converted to Mb, rounded up).
  size_t total_mem_usage = ( paths_size_sum + paths_rc_size_sum ) * sizeof(TAG);
  total_mem_usage >>= 20;
  total_mem_usage++;
  
  cout << Date() << ": Creating pathsdb (will be a vector with size = "
       << ( paths_size_sum + paths_rc_size_sum ) << " and mem usage "
       << total_mem_usage << " Mb)" << endl;
  vec<TAG> pathsdb;
  
  if ( paths_rc.empty( ) ) CreateDatabase( paths, pathsdb );
  else {
    
    // The following code is just like CreateDatabase( paths, paths_rc, pathsdb)
    // except that it destroys paths and paths_rc.
    pathsdb.reserve( paths_size_sum + paths_rc_size_sum );
    for ( size_t i = 0; i < paths.size(); i++ )
      paths[i].AppendToDatabase( pathsdb, i );
    paths.clear().shrink_to_fit();
    
    for ( size_t i = 0; i < paths_rc.size(); i++ )
      paths_rc[i].AppendToDatabase( pathsdb, ~i ); // ~i == -i-1
    paths_rc.clear().shrink_to_fit();
    
    cout << Date() << ": Preparing pathsdb" << endl;
    Prepare(pathsdb);
  }
  
  cout << Date() << ": Writing pathsdb to file" << endl;
  BinaryWrite3( pathsdb_file, pathsdb );
  return;
}

