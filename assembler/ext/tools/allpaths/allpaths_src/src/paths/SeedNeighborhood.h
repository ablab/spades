///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SEED_NEIGHBORHOOD
#define SEED_NEIGHBORHOOD

/*******************************************************************************
 *
 *        SEED NEIGHBORHOOD CLASS
 *
 * This class represents a "seed" unipath and its local "neighborhood".  Each
 * seed is a carefully-selected unipath from the global assembly, chosen to have
 * copy number 1 (CN-1) and be sufficiently close to other seeds so that their
 * neighborhoods will eventually overlap.
 *
 * The task of this class is to perform a local assembly of the neighborhood,
 * using only the reads and pairs believed to belong here.  The assembled
 * neighborhood will then be merged with other seeds' neighborhoods to create
 * the final assemble.
 *
 * The assembly of a seed into its local neighborhood takes the following steps:
 * 1. MakeUnipathCloud: Find all nearby CN-1 unipaths.  These are called the
 *    "unipath cloud."
 * 2. MakeReadCloud: Find all reads lying on unipaths in the unipath cloud.
 *    Also compile a set of pairs for which EITHER read lies on a cloud unipath,
 *    and put BOTH that pair's reads into the read cloud.  This "localizes" the
 *    reads.
 * 3. MakeLocalAssembly: Take the reads from the read cloud and make them into
 *    paths, then unipaths, then a HyperKmerPath.
 * 4. MakeInsertWalks: Attempt to close local pairs by walking them with local
 *    reads.  Create a new HyperBasevector accordingly.
 *
 * All member functions are defined in either SeedNeighborhood.cc or
 * SeedNeighborhoodAnnex.cc, and they are also documented there.
 *
 * Josh Burton
 * November 2008
 *
 ******************************************************************************/


#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "TaskTimer.h"
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h" // digraph
#include "math/NStatsTools.h" // PrintBasicNStats
#include "paths/HyperBasevector.h"
#include "paths/PdfEntry.h"
#include "paths/SeedNeighborhoodAnnex.h"
#include "paths/UnipathNhoodCommon.h" // fsepdev, ustart
#include "paths/simulation/Placement.h"


class SeedNeighborhood
{
public:
  
  // Constructor: immediately loads all global data structures (by reference)
  SeedNeighborhood( const int K,
		    const int ploidy,
		    const int mcn_other,
		    const vecKmerPath * paths,
		    const vecKmerPath * paths_rc,
		    const vec<tagged_rpint> * pathsdb,
		    const vec<int> * unipath_lengths,
		    const vec<int> * read_lengths,
		    const vec<int> * predicted_CNs,
		    const PairsManager * pairs,
		    const vec<ReadLocationLG> * unilocs,
		    const vec<longlong> * unilocs_index,
                    const vec<String>& LOCAL_DUMP = vec<String>( ) );
  
  ~SeedNeighborhood( );
  
  /* ACTION FUNCTIONS
   *
   * These functions work together to assemble a HyperKmerPath that is local to
   * this neighborhood.  It is best to run them in the order they appear here.
   *
   ****************************************************************************/
  
  // REQUIRED
  void SetSeedID( const int seed_ID ) { _seed_ID = seed_ID; }
  
  // Optional: Set up logging for this SeedNeighborhood
  // All diagnostic output goes to the file "logfile", except in two cases:
  // logfile = "stdout": output goes to stdout
  // logfile = "": output silenced (as if SetLog were never called)
  void SetLog( const String & logfile );
  
  // Optional: Activate eval_subdir, which will be used to dump
  // (large!) files for the evaluation of the nhood.
  void SetEvalSubdir( const String &evaldir );
  
  // Make the cloud of unipaths in the vicinity of this seed
  void MakeUnipathCloud( const digraphE<fsepdev> & LG,
			 const vec<int> & predicted_CNs,
			 const vec<Bool> & branches );
  
  // Find the reads and pairs in this neighborhood
  void MakeReadCloud( const vecbvec &global_reads_bases,
		      const vecbvec &unibases, const Bool& LOCAL_PRIMARY );
  
  // Assemble this seed's local reads into a unipath graph and HyperKmerPath
  void MakeLocalAssembly( const vecbasevector& unibases,
     const vec< vec<int> >& unibases_next, const vec<int>& to_rc, const int pass );
  
  // Experimental code from David (add global connections to _reads).
  void AddGlobalConnections( const vecbasevector& unibases,
     const vec< vec<int> >& unibases_next, const vec<int>& to_rc );

  // Walk inserts in this neighborhood with the help of fragments
  // returns 'false' on irregular exit (e.g. timeout)
  bool MakeInsertWalks( );
  
  // Write the HyperBasevector of this assembly to files
  void Write( size_t iSeed, BinaryWriter& hbvWriter,
              bool LOCAL_DOT, bool LOCAL_FASTA );
  
  
  /* EVALUATION FUNCTIONS
   *
   * These functions evaluate the local assembly in progress, writing useful
   * diagnostic output to the logfile.  They are the only functions that allow
   * the reference genome to be loaded in - and they are all const, guaranteeing
   * that the genome will not affect the assembly.
   * If you only want non-reference-based evaluation, you can call these
   * functions with empty unput objects.
   *
   ****************************************************************************/
  
  void EvalUnipathCloud( const vecbasevector & genome, 
			 const VecPlacementVec& unipath_POGs ) const;
  void EvalReadCloud( const vecbasevector & genome, 
		      const VecPlacementVec& unipath_POGs,
		      const vec<ReadLocationLG> & read_POGs ) const;
  void EvalLocalAssembly( ) const;
  void EvalInsertWalks( ) const;
  
  
  /* HELPER FUNCTIONS
   *
   * These functions are helpers for the action functions, above.
   * They perform much of the algorithmic work of assembling a neighborhood.
   * They are declared here but defined in SeedNeighborhoodAnnex.cc.
   *
   * Some of these functions modify SeedNeighborhood member variables and
   * return void.  Others are const and return an object, which is their output.
   *
   ****************************************************************************/
public:
  void FindUnipathCloud( const digraphE<fsepdev> & LG, const vec<int> & predicted_CNs, const vec<Bool> & branches, const int MAX_DEV );
  void ExpandUnipathCloud( const digraphE<fsepdev> & LG, const vec<int> & predicted_CNs, const int MAX_DEV );
  void FindPrimaryReadCloud( const vecbvec& global_reads_bases,
       const Bool& LOCAL_PRIMARY );
  void FindSecondaryReadCloud( );
  vec<longlong> SelectPairsToWalk( const int n_to_select, int & n_logical_pairs ) const;
  HyperKmerPath MakeAcyclicHKP( const vecKmerPath & new_unipaths, const vecbasevector & new_unibases, const int min_size ) const;
  HyperKmerPath WalkInserts( const vec<longlong> & selected_pair_IDs, const int STRETCH, TaskTimer & timer, int & n_to_walk, int & n_walked, vec<Bool>& walked ) const;
  void MergeInsertWalks( const vec<HyperKmerPath>& HKPs, const vecKmerPath& new_unipaths, const vecbasevector & new_unibases );
  
  
  /* QUERY FUNCTIONS
   *
   * These functions are all const and public.
   *
   ****************************************************************************/
  vec<ustart> GetCloudUnipaths() const { return _cloud_unipaths; }

  
  
  
  /* GLOBAL DATA STRUCTURES
   *
   * These objects have implicit indices/kmer_ids/etc. that apply to the general
   * assembly, not to this specific SeedNeighborhood.
   * They are loaded into the SeedNeighborhood object when it is initialized,
   * and cannot be modified thereafter.
   * Note that most of these are pointers - the SeedNeighborhood class does NOT
   * own the original objects!  If the originals are de-allocated or moved,
   * SeedNeighborhood will seg fault.
   *
   ****************************************************************************/
private:
  const int _K;
  const int _ploidy;
  const int _mcn_other;
  const vecKmerPath         * _global_paths;
  const vecKmerPath         * _global_paths_rc;
  const vec<tagged_rpint>   * _global_pathsdb;
  const vec<int>            * _global_unipath_lengths;
  const vec<int>            * _global_read_lengths;
  const vec<int>            * _global_predicted_CNs;
  const PairsManager        * _global_pairs;
  const vec<ReadLocationLG> * _global_unilocs;
  const vec<longlong>       * _global_unilocs_index;
  const vec<String>         _LOCAL_DUMP;
  
  
  /* LOCAL DATA STRUCTURES
   *
   * These data structures are created by the SeedNeighborhood object itself.
   * They are used in the creation of the neighborhood's assembly
   * Do not confuse these objects with global objects of similar names!
   * These local objects may (or may not) use local read numbering or kmer_ids.
   *
   ****************************************************************************/
private:
  // TODO: potentially dangerous truncation of index by these int members
  int _seed_ID; // in global space
  int _n_cloud_unipaths, _n_reads, _n_pairs, _n_local_unipaths;
  
  vec<ustart> _cloud_unipaths; // Cloud unipaths - with global unipath numbering
  
  vec< pair<longlong, Bool> > _read_IDs; // Local read IDs, with orientations
  vec<int> _read_locs; // (Estimated) locations of reads in local neighborhood
  vecbasevector _reads; // Local reads
  vecbasevector _reads_fw;
  
  vec<longlong> _local_to_global; // Read ID map
  
  vec<longlong> _pair_IDs; // Map from local to global pair numbering
  
  digraph _AG; // Local unipath adjacency graph
  vecKmerPath _local_unipaths;
  vec<Bool> _local_unipaths_fw;
  vecKmerPath _paths, _paths_rc; // Local paths
  vec<tagged_rpint> _pathsdb, _local_unipathsdb;
  KmerBaseBroker _kbb; // Local KmerBaseBroker
  HyperKmerPath _hkp; // Local HyperKmerPath
  // Final output: the local HyperBasevector containing the merged assembly
  HyperBasevector _hbv;
  
  
  
  /* Tallies and runtime counters */
public:
  int _n_inserts_walked, _n_inserts_total;
  double _T_pathing, _T_insert_walking, _T_insert_merging;
  
  
  
  
  /* Heuristic parameters */
  // These are all adapted from arguments to LocalizeReads (not LG)
public:
  
  // Neighborhood radius: Unipaths are only added to this neighborhood if their
  // distance to the seed unipath (expressed as a gap) is less than NHOD_RADIUS
  static int _NHOOD_RADIUS;
  static int _NHOOD_RADIUS_INTERNAL;
  static int _NHOOD_GLUEPERFECT_SIZE;
  static int _NHOOD_VERBOSITY;
  static Bool _DUMP_UNMERGED;
  static Bool _DUMP_LOCAL_UNIBASES;
  static Bool _USE_ACYCLIC;
  String _outhead;
  
private:
  ostream * _log; // logfile output stream
  String * _eval_subdir; // save here intermediate (possibly large) files
  void LogDate( const String& s, const int n_newlines = 0 ) const; // write to log

};

// Definitions of verbosity control bits
#define VERBOSITY_UNIPATH_CLOUD		1


#endif
