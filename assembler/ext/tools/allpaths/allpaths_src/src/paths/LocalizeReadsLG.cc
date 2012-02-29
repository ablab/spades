///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library PTHREAD

/* L      OOO   CCC    A   L      III  ZZZZZ EEEEE RRRR  EEEEE   A   DDDD   SSS
 * L     O   O C   C  A A  L       I       Z E     R   R E      A A   D  D S   S
 * L     O   O C     A   A L       I      Z  E     R   R E     A   A  D  D S
 * L     O   O C     A   A L       I      Z  E     R   R E     A   A  D  D S
 * L     O   O C     A   A L       I     Z   EEEE  RRRR  EEEE  A   A  D  D  SSS
 * L     O   O C     AAAAA L       I    Z    E     R R   E     AAAAA  D  D     S
 * L     O   O C     A   A L       I    Z    E     R  R  E     A   A  D  D     S
 * L     O   O C   C A   A L       I   Z     E     R   R E     A   A  D  D S   S
 * LLLLL  OOO   CCC  A   A LLLLL  III  ZZZZZ EEEEE R   R EEEEE A   A DDDD   SSS
 *
 * This is the LG (Large Genome) version of LocalizeReads.  It is a complete
 * overhaul of the original LocalizeReads, intended to be a clean-slate
 * re-envisioning of the central algorithm and an adaptation to the scenario of
 * the new 100-bp Illumina reads.  It belongs in the RunAllPathsLG pipeline.
 *
 * This module contains much of the high-level functionality of the original
 * LocalizeReads. The unipath seed selection at the beginning has been branched
 * off into UnipathSeedsLG, and the neighborhood merging at the end has been
 * branched off into MergeNeighborhoods.  However, the local neighborhood
 * assemblies - one for each seed unipath - are created here.
 *
 * The central subclass is the SeedNeighborhood class, which handles the
 * implementation of the local assemblies.  This module is mostly a client for
 * SeedNeighborhood.  Here files are loaded and the SeedNeighborhoods are
 * created and supplied with the necessary global data.
 *
 * You can parallelize LocalizeReadsLG by setting the NUM_THREADS parameter.
 * When LocalizeReadsLG is running, it will assemble NUM_THREADS local neighborhoods
 * in parallel and will thus use as many as NUM_THREADS processors; however, it
 * runs entirely within a single kernel process.  Large data structures (such
 * as the reads) are shared between neighborhoods, making the parallelization
 * memory-efficient.  Parallelization is implemented with the Worklist class.
 *
 * You can also run a subset of the local assemblies with the SEED_IDS option.
 * For example, if you set SEED_IDS="{0,1,3}", then only neighborhoods 0, 1, 3
 * will be assembled.  (Note that these are NOT unipath IDs - they are seed IDs,
 * corresponding to locations in the file seeds.ids.)  You can set SEED_IDS=all
 * to override the seed statuses and assemble all seeds.
 *
 * [Filipe Ribeiro 02.2009]
 * As each neighborhood is processed the file 'seeds.status' is updated.
 * This file stores completion information for each seed as a sequence of
 * keywords for each seed:
 *
 * 'done'   : neighborhood succesfully processed.
 * 'todo'   : neighborhood still unprocessed.
 * 'running': neighborhood is being processed or was being processed and there
 *            was an interruption.
 * 'later'  : neighborhood might have led to a crash and will be processed at
 *            the end.
 *
 * Josh Burton
 * November 2008
 *
 ******************************************************************************/

#include <pthread.h>
#include "MainTools.h"
#include "PairsManager.h" // PairsManager
#include "ParseSet.h"
#include "ReadLocationLG.h"
#include "TaskTimer.h" // TaskTimer
#include "feudal/BinaryStream.h"
#include "graph/Digraph.h" // digraph
#include "paths/GetNexts.h"
#include "paths/PdfEntry.h" // pdf_entry
#include "paths/SeedNeighborhood.h" // SeedNeighborhood
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h" // BuildUnipathAdjacencyGraph
#include "paths/reporting/ReportWalkingRatesCore.h"
#include "paths/simulation/Placement.h" // placement
#include "system/Worklist.h" // Worklist
#include "system/LockedData.h" // LockedData

static inline 
String Tag(String S = "LRLG") { return Date() + " (" + S + "): "; } 



/******************************************************************************
 *
 *
 *
 ******************************************************************************/


class SeedStore : public LockedData
{
private:
  static const size_t _SEEDS_PER_DIR = 1000;

  String      _sub_dir;

  String      _seed_status_file;
  vec<String> _seed_status;

  size_t      _num_seeds;
  vec<int>    _seed_IDs;
  vec<int>    _i_seeds_todo;

  size_t      _num_no_writes;
  size_t      _num_timedouts;

  String       mHBVFile;
  BinaryWriter mHBVWriter;

  SeedStore(const SeedStore & b); // unimplemented -- no copying
  SeedStore & operator=(const SeedStore & b); // unimplemented -- no copying

public:

  SeedStore(const String & sub_dir,
            const String SEED_IDS,
            bool makeSeedsDirs )
    : _sub_dir(sub_dir),
      _seed_status_file(sub_dir + "/seeds.status"),
      mHBVFile(sub_dir + "/localized.hbvs"),
      mHBVWriter(mHBVFile.c_str())
  {
    String seeds_file = sub_dir + "/seeds.ids";
    cout << Tag() << "Reading seed IDs from file " + seeds_file + "..." << endl;
    FstreamReadOrFatal(seeds_file, _seed_IDs);
    _num_seeds = _seed_IDs.size();
    mHBVWriter.write(_num_seeds);

    // Make the directories that will contain output from the local assemblies.
    // We do this now because we do NOT want to do it in parallel.
    if ( makeSeedsDirs )
        MakeSeedsDirs();

    ReadSeedStatus();

    if (SEED_IDS != "") {  // If the arg SEED_IDS was given, ignore seed_status entirely
      SeedsToDoFromChosen(SEED_IDS);
    }
    else {                // use seed_status = 'todo' to determine which seeds to process
      SeedsToDoFromStatus();
    }
    cout << Tag() << "Found " << _num_seeds << " seeds; listed " << GetNumSeedsToDo() << " as to-do." << endl;

    _num_no_writes = 0;
    _num_timedouts = 0;

  }

  void close() {
    mHBVWriter.write(~0UL);
    mHBVWriter.close();
    WriteSeedStatus(0);
  }

  BinaryWriter& getWriter() { return mHBVWriter; }

  size_t GetNumSeeds() { return _num_seeds; }
  size_t GetNumSeedsToDo() { return _i_seeds_todo.size(); }
  size_t GetSeedID(const size_t i_seed) { return _seed_IDs[i_seed]; }
  const vec<int> &GetSeedIDs( ) const { return _seed_IDs; }

  String GetSeedDir(const size_t i_seed)
  { return _sub_dir + "/seed/" + GetDirFromKey(i_seed, _SEEDS_PER_DIR); }

private:


  void MakeSeedsDirs()
  {
    Mkdir777(_sub_dir + "/seed");
    for (size_t i_seed = 0; i_seed < _num_seeds; i_seed += _SEEDS_PER_DIR)
      Mkdir777(GetSeedDir(i_seed));
  }


  void SeedsToDoFromChosen(const String & SEED_IDS)
  {
    cout << Tag() << "Only processing seeds in SEED_IDS=" << SEED_IDS << endl;
    ParseIntSet( SEED_IDS, _i_seeds_todo );
    UniqueSort( _i_seeds_todo );

    if ( (size_t)_i_seeds_todo.back() >= _num_seeds)
      FatalErr( "All integers in SEED_IDS must be less than n_seeds = "
                << _num_seeds << ". (You gave " << SEED_IDS << ".)" );
  }


  void SeedsToDoFromStatus()
  {
    _i_seeds_todo.clear();

    if (1) {  // hardwired choice between processing all not-done seeds in order ...

      for (size_t i = 0; i != _num_seeds; i++) // add "todo" seeds
        if (_seed_status[i] != "done")
          _i_seeds_todo.push_back(i);
    }
    else {  // ... or leaving "problematic" seeds for later.

      for (size_t i = 0; i != _num_seeds; i++) // add "todo" seeds
        if (_seed_status[i] == "todo")
          _i_seeds_todo.push_back(i);

      for (size_t i = 0; i != _num_seeds; i++) // add "later" seeds to be run at the end
        if (_seed_status[i] == "later")
          _i_seeds_todo.push_back(i);
    }
  }



  void ReadSeedStatus()
  {
    // Read list of seed status (if such list exists)
    // the various status possible are:
    //
    // 'todo'    : not processed yet; previous run crashed
    // 'done'    : already processed
    // 'running' : was being processed and there was a crash
    // 'later'   : previously ran but potentially problematic
    //
    Locker lock(*this);

    cout << Tag() << "Reading the seed status from file " + _seed_status_file + "." << endl;

    if (IsRegularFile(_seed_status_file) && FileSize(_seed_status_file) == 0)
      _seed_status.resize(0);
    else FstreamReadOrContinue(_seed_status_file, _seed_status);

    cout << Tag() << "Found " <<  _seed_status.size() << " seed statuses." << endl;

    // If the previous run was interrupted some of the statuses might be 'running'
    // Here these statuses are updated to 'later'
    for (size_t i = 0; i != _seed_status.size(); i++) {
      if (_seed_status[i] == "running")   // was running but crashed; do at the end
        _seed_status[i] = "later";
      else if (_seed_status[i] != "todo" &&
               _seed_status[i] != "done" &&    // (*seed_status)[i] is unknown; set it to 'later'
               _seed_status[i] != "later")     // prevents issues with corrupted seeds.status file
        _seed_status[i] = "later";
    }

    // make sure all seeds are covered in case seed_status was corrupted and is too short
    for (size_t i = _seed_status.size(); i != _num_seeds; i++)
      _seed_status.push_back("todo");
  }

  void WriteSeedStatus(const size_t num_no_writes_min = 64)
  {
    if (_num_no_writes > num_no_writes_min || _num_timedouts > 0) {
      Locker lock(*this);
      FstreamWrite(_seed_status_file, _seed_status);
      _num_no_writes = _num_timedouts = 0;
    }
    else {
      _num_no_writes++;
    }

  }


public:

  bool GetNextSeedIndex(size_t * i_seed)
  {
    Locker lock(*this);
    if (!_i_seeds_todo.empty()) {
      *i_seed = _i_seeds_todo.back();
      _i_seeds_todo.pop_back();
      return true;
    }
    return false;
  }

  void SetSeedStatusDone(const size_t i_seed)
  {
    {
      Locker lock(*this);
      _seed_status[i_seed] = "done";
    }
    WriteSeedStatus();
  }

  void SetSeedStatusRunning(const size_t i_seed)
  {
    {
      Locker lock(*this);
      _seed_status[i_seed] = "running";
    }
    WriteSeedStatus();
  }

  void SetSeedStatusTimedOut(const size_t i_seed)
  {
    {
      Locker lock(*this);
      _seed_status[i_seed] = "timedout";
      _num_timedouts++;
    }
    WriteSeedStatus();
  }


};






class LRLGCounters : public LockedData
{
public:
  int    n_seeds_done;
  int    n_inserts_walked;
  int    n_inserts_total;
  double T_pathing;
  double T_insert_walking;
  double T_insert_merging;
  vec< pair<int,int> > walks_info;   // inserts_(walked,total) for each nhood

  LRLGCounters( size_t n_seeds ) 
    : n_seeds_done(0), 
      n_inserts_walked(0),
      n_inserts_total(0),
      T_pathing(0.0),
      T_insert_walking(0.0),
      T_insert_merging(0.0)
  {
    walks_info.resize( n_seeds, make_pair( -1, -1 ) );
  }

};




/*****************************************************************************
 *
 *   Processor class for handling multi-threading of the local
 *   assemblies.  For an explanation of this, see system/Worklist.h.
 *
 *****************************************************************************/

#include <sys/resource.h>

class Processor
{
private:
  // Large objects, representing the global assembly.
  // These objects are loaded (by const reference) into each SeedNeighborhood.
  // They should NOT be changed here - if they do not hold the same value for
  // each local assembly, the result will be nonsensical.
  SeedStore & _seeds;

  int _K;
  int _ploidy;
  int _mcn_other;
  vecKmerPath * _paths;
  vecKmerPath * _paths_rc;
  vec<tagged_rpint> * _pathsdb;
  vec<int> * _unipath_lengths;
  vec<int> * _read_lengths;
  PairsManager * _pairs;
  vec<ReadLocationLG> * _unilocs;
  vec<longlong> * _unilocs_index;
  digraphE<fsepdev> * _LG;
  vec<int> * _predicted_CNs;
  vec<Bool> * _branches;
  BaseVecVec * _reads_bases;
  BaseVecVec * _unibases;
  vec< vec<int> > * _unibases_next;
  vec<int> * _to_rc;
  Bool _EVAL;
  Bool _READ_CLOUD_EVAL;
  Bool _LOCAL_DOT;
  vec<String> _LOCAL_DUMP;
  Bool _LOCAL_FASTA;
  Bool _LOCAL_PRIMARY;
  Bool _BRIDGE_BETWEEN_CLOUD_UNIPATHS;
  Bool _LOCAL_LOG;
  Bool _nhood_eval;
  BaseVecVec * _genome;
  VecPlacementVec * _unipath_POGs;
  vec<ReadLocationLG> * _read_POGs;
  LRLGCounters * _counters;
  ofstream * _log_fs;

public:
  Processor(SeedStore & seeds,
            const int K,
            const int ploidy,
            const int mcn_other,
            vecKmerPath * paths,
            vecKmerPath * paths_rc,
            vec<tagged_rpint> * pathsdb,
            vec<int> * unipath_lengths,
            vec<int> * read_lengths,
            PairsManager * pairs,
            vec<ReadLocationLG> * unilocs,
            vec<longlong> * unilocs_index,
            digraphE<fsepdev> * LG,
            vec<int> * predicted_CNs,
            vec<Bool> * branches,
            BaseVecVec * reads_bases,
            BaseVecVec * unibases,
            vec< vec<int> > * unibases_next,
            vec<int> * to_rc,
            Bool EVAL,
            Bool READ_CLOUD_EVAL,
            Bool LOCAL_DOT,
            vec<String> LOCAL_DUMP,
            Bool LOCAL_FASTA,
            Bool LOCAL_PRIMARY,
            Bool BRIDGE_BETWEEN_CLOUD_UNIPATHS,
            Bool LOCAL_LOG,
            Bool nhood_eval,
            BaseVecVec * genome,
            VecPlacementVec * unipath_POGs,
            vec<ReadLocationLG> * read_POGs,
            LRLGCounters * counters,
            ofstream * log_fs)
    : _seeds(seeds),
      _K(K),
      _ploidy(ploidy),
      _mcn_other(mcn_other),
      _paths(paths),
      _paths_rc(paths_rc),
      _pathsdb(pathsdb),
      _unipath_lengths(unipath_lengths),
      _read_lengths(read_lengths),
      _pairs(pairs),
      _unilocs(unilocs),
      _unilocs_index(unilocs_index),
      _LG(LG),
      _predicted_CNs(predicted_CNs),
      _branches(branches),
      _reads_bases(reads_bases),
      _unibases(unibases),
      _unibases_next(unibases_next),
      _to_rc(to_rc),
      _EVAL(EVAL),
      _READ_CLOUD_EVAL(READ_CLOUD_EVAL),
      _LOCAL_DOT(LOCAL_DOT),
      _LOCAL_DUMP(LOCAL_DUMP),
      _LOCAL_FASTA(LOCAL_FASTA),
      _LOCAL_PRIMARY(LOCAL_PRIMARY),
      _BRIDGE_BETWEEN_CLOUD_UNIPATHS(BRIDGE_BETWEEN_CLOUD_UNIPATHS),
      _LOCAL_LOG(LOCAL_LOG),
      _nhood_eval(nhood_eval),
      _genome(genome),
      _unipath_POGs(unipath_POGs),
      _read_POGs(read_POGs),
      _counters(counters),
      _log_fs(log_fs)
  {}

  // operator(): This is the function called in parallel by the Worklist
  // When called, it will assemble the seed neighborhood #s
  void operator() (size_t /*thread_ID*/)
  {
    size_t i_seed = 0;
    while ( _seeds.GetNextSeedIndex( &i_seed ) ) {
      TaskTimer prun_timer;
      prun_timer.Start( );
    
      size_t seed_ID = _seeds.GetSeedID(i_seed);

      String seed_string = Tag("LRLG.P") + "nhood #" + ToString(i_seed) + " seed_ID= " + ToString(seed_ID);
      (*_log_fs) << (seed_string + "  starting local assembly") << endl;
      
      _seeds.SetSeedStatusRunning(i_seed);

      // Create SeedNeighborhood object to perform the local assembly
      SeedNeighborhood nhood(_K, _ploidy, _mcn_other,
                             _paths, _paths_rc, _pathsdb, _unipath_lengths,
                             _read_lengths, _predicted_CNs, _pairs, _unilocs,
                             _unilocs_index, _LOCAL_DUMP);

      nhood.SetSeedID(seed_ID);

      // Find the seed_dir, where we will be placing output files from this nhood.
      // This directory should already exist (DO NOT call Mkdir in parallel code!)
      {
	Locker locker( _seeds );
	String seed_dir = _seeds.GetSeedDir(i_seed);
	nhood._outhead = seed_dir + "/" + ToString(i_seed);
	String log_file = nhood._outhead + ".log";
	String hbv_file = nhood._outhead + ".hbv";
	if (_LOCAL_LOG) nhood.SetLog(log_file);
      }
      
      bool completed = false;

      // Find the neighborhood (unipaths, reads) of this seed
      nhood.MakeUnipathCloud(*_LG, *_predicted_CNs, *_branches);
      if (_EVAL)
	nhood.EvalUnipathCloud(*_genome, *_unipath_POGs);
      
      nhood.MakeReadCloud(*_reads_bases, *_unibases, _LOCAL_PRIMARY);
      if (_EVAL && _READ_CLOUD_EVAL)
	nhood.EvalReadCloud(*_genome, *_unipath_POGs, *_read_POGs);

      // Assemble the nhood's reads into a local unipath graph and HyperKmerPath
      if (_BRIDGE_BETWEEN_CLOUD_UNIPATHS)
	nhood.MakeLocalAssembly( *_unibases, *_unibases_next, *_to_rc, 1 );
      nhood.MakeLocalAssembly( *_unibases, *_unibases_next, *_to_rc, 2 );
      if (_EVAL)
	nhood.EvalLocalAssembly();
      
      // Walk inserts in this neighborhood with the help of fragments;
      // then merge the walked inserts into the local assembly
      completed = nhood.MakeInsertWalks( );
      if (_EVAL && completed)
	nhood.EvalInsertWalks();
      
      // Each locally assembled neighborhood is represented as a
      // HyperBasevector (which is a HyperKmerPath that is independent
      // of kmer numbering) The HyperBasevectors are saved into files
      // in SUBDIR
      prun_timer.Stop( );
      String run_time = ToString( prun_timer.Elapsed( ), 2 ) + "s";
      
      if (completed) {
	String str_done
	  = ToString( seed_string ) + "  local assembly done ("
	  + ToString( run_time ) + "): "
	  + ToString( nhood._n_inserts_walked ) + " inserts walked, out of "
	  + ToString( nhood._n_inserts_total ) + " selected";
        (*_log_fs) << str_done << endl;
	{
	  Locker locker(_seeds);
	  nhood.Write(i_seed, _seeds.getWriter(), _LOCAL_DOT, _LOCAL_FASTA);
 	}
	_seeds.SetSeedStatusDone(i_seed);
      }
      else {
        (*_log_fs) << (seed_string  + "  local assembly TIMED OUT (" + run_time + ")") << endl;
        _seeds.SetSeedStatusTimedOut(i_seed);
      }
      
      // Update global counter variables.
      {
        Locker l(*_counters); // destructor unlocks

        _counters->n_seeds_done++;
        _counters->n_inserts_walked += nhood._n_inserts_walked;
        _counters->n_inserts_total  += nhood._n_inserts_total;
        // There is one time counter for each of the three bottleneck steps in
        // local neighborhood assembly: pathing, insert walking, and insert merging.
        _counters->T_pathing        += nhood._T_pathing;
        _counters->T_insert_walking += nhood._T_insert_walking;
        _counters->T_insert_merging += nhood._T_insert_merging;
	// More granular informations on insert walking.
	_counters->walks_info[i_seed]
	  = make_pair( nhood._n_inserts_walked, nhood._n_inserts_total );
      }
    }
  }
};


/*****************************************************************************
 *
 *   main starts here
 *
 *****************************************************************************/

int main(int argc, char *argv[])
{
  RunTime();


  // Command-line arguments
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault(SUBDIR, "test");
  CommandArgument_Int(K);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Int_OrDefault_Valid_Doc(PLOIDY, -1, "{-1,1,2}",
    "Used in the absence of a ploidy file: 1 for Haploid, 2 for Diploid");
  CommandArgument_String_OrDefault(READS, "reads");
  // Arguments to do chunks of the pipeline (or not)
  CommandArgument_Bool_OrDefault(USE_TRUTH, True);
  CommandArgument_Bool_OrDefault(EVAL, True);
  CommandArgument_Bool_OrDefault_Doc(NHOOD_EVAL, False,
    "[Optional] Run extra eval code on all the neighboroods.");
  CommandArgument_Bool_OrDefault_Doc(READ_CLOUD_EVAL, False,
    "[Optional] Run extra read cloud eval code (only if EVAL=True).\nNB: it needs the output from GenomeReadLocs!\n");
  CommandArgument_String_OrDefault_Doc(SEED_IDS, "",
    "[Optional] Assemble only these IDs from the seed list.");
  CommandArgument_Bool_OrDefault_Doc(OVERRIDE_TIMEOUTS, False,
    "[Optional] If any neighborhoods time out, re-run them afterward, in non-parallel.  Only operative if SEED_IDS=\"\".");

  // Heuristic parameters
  CommandArgument_Bool_OrDefault_Doc(USE_THEO_GRAPHS, False, "Use the file unipath_link_graph.theo, created by TheoreticalUnipathLinkGraph");
  CommandArgument_Int_OrDefault(NHOOD_RADIUS, 20000);
  CommandArgument_Int_OrDefault(NHOOD_RADIUS_INTERNAL, 10000);
  CommandArgument_Int_OrDefault_Doc(NHOOD_GLUEPERFECT_SIZE, 0,
                                 "if >0, try to merge duplicate edge sequence of this size or greater");
  CommandArgument_Int_OrDefault_Doc(MAX_COPY_NUMBER_OTHER, 0, "Control the copy numer of unipaths entering the neighborhood.");
  CommandArgument_String_OrDefault_Doc(LOCAL_DUMP, "",
     "List of zero or more of the following, to generate various files for\n"
     "each local assembly:"
     "\nLOG - make log file n.log"
     "\nFASTA - make local assembly fasta file n.hbv.fasta"
     "\nDOT - make local assembly dot file n.hbv.dot"
     "\nUNIBASES - make local unibase files n.local_unibases.{dot,fasta}"
     "\nPRIMARY - dump primary read close n.primary.fasta"
     "\nUNMERGED - dump unmerged assembly n.unmerged.{dot,fasta}"
     "\nSELECTED_PAIRS - dump selected pairs n.selected_pairs.fasta");
  CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0,
       "Verbosity bits: 1=UNIPATH_CLOUD.");
  CommandArgument_Bool_OrDefault_Doc(USE_ACYCLIC, True,
       "Use acyclic stuff in constructing neighborhoods.");
  CommandArgument_Bool_OrDefault_Doc(BRIDGE_BETWEEN_CLOUD_UNIPATHS, False,
       "Attempt to walk in the global unipath between cloud unipaths.");
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  cout << Tag() << "Beginning LocalizeReadsLG..." << endl;
  cout << Tag() << "Using " << NUM_THREADS << " threads." << endl;
  TaskTimer timer;
  timer.Start();

  // Parse LOCAL_DUMP.

  vec<String> local_dump;
  ParseStringSet( LOCAL_DUMP, local_dump );
  Bool LOCAL_FASTA = Member( local_dump, String("FASTA") );
  Bool LOCAL_LOG = Member( local_dump, String("LOG") );
  Bool LOCAL_DOT = Member( local_dump, String("DOT") );
  Bool DUMP_LOCAL_UNIBASES = Member( local_dump, String("UNIBASES") );
  Bool LOCAL_PRIMARY = Member( local_dump, String("PRIMARY") );
  Bool DUMP_UNMERGED = Member( local_dump, String("UNMERGED") );

  /*****************************************************************************
   *
   *        LOAD INPUT FILES
   *
   ****************************************************************************/

  cout << Tag() << "Loading input files" << endl;

  // Filenames
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String file_head = run_dir + "/" + READS;
  String kK = ".k" + ToString(K);
  String data_dir = PRE + "/" + DATA;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;


  // check the ploidy input
  if (IsRegularFile(data_dir + "/ploidy")) {
    int ploidy_from_file = StringOfFile(data_dir + "/ploidy", 1).Int();
    if (PLOIDY >= 1){
      if (ploidy_from_file != PLOIDY)
        InputErr("PLOIDY set as a command option but different from a value in ploidy file.");
    }else
      PLOIDY = ploidy_from_file;
  }
  if (PLOIDY < 1 || PLOIDY > 2)
    InputErr("Unsupported ploidy value file - must be 1 or 2");

  // Read bases
  vecbasevector reads_bases(file_head + ".fastb");

  // Read paths [_rc, db]
  vecKmerPath paths   (file_head + ".paths" + kK);
  vecKmerPath paths_rc(file_head + ".paths_rc" + kK);
  vec<tagged_rpint> pathsdb;
  BinaryRead3(file_head + ".pathsdb" + kK, pathsdb);

  // Unibases, unipaths [db]
  vecbasevector unibases(file_head + ".unibases" + kK);
  vecKmerPath unipaths(file_head + ".unipaths" + kK);
  vec<tagged_rpint> unipathsdb;
  BinaryRead3(file_head + ".unipathsdb" + kK, unipathsdb);

  vec< vec<int> > unibases_next( unibases.size( ) );
  vec<int> to_rc;
  if (BRIDGE_BETWEEN_CLOUD_UNIPATHS) 
  {    GetNexts( K, unibases, unibases_next );
       UnibaseInvolution( unibases, to_rc, K );    }

  // Unipath link graph
  String graph_file_infix = USE_THEO_GRAPHS ? ".theo" : "";
  digraphE<fsepdev> LG;
  BinaryRead(sub_dir + "/unipath_link_graph" + graph_file_infix + ".cloud" + kK, LG);

  // Unipath copy numbers (predicted)
  VecPdfEntryVec CNs((file_head + ".unipaths.predicted_count" + kK).c_str());

  // Unilocs (i.e., read locations on unipaths)
  vec<ReadLocationLG> unilocs;
  BinaryRead2(file_head + ".unilocs." + ToString(K), unilocs);

  // Calculate read lengths
  vec<int> read_lengths = KmerPathSeqLength( paths, K );

  // Read pairing information
  PairsManager pairs(file_head + ".pairs");
  pairs.makeCache();

  // Variables containing genome information - may be empty!
  // These are filled only if USE_TRUTH=True (and the unipath placements only
  // if the locs file has been created by FindUnipathLocs)
  vecbasevector genome;
  VecPlacementVec unipath_POGs; // unipath Placements On Genome
  vec<ReadLocationLG> read_POGs; // read Placements On Genome

  if (USE_TRUTH) {
    genome.ReadAll(data_dir + "/genome.fastb");
    if (IsRegularFile(file_head + ".unipaths" + kK + ".locs"))
      unipath_POGs.ReadAll(file_head + ".unipaths" + kK + ".locs");
  }

  if (READ_CLOUD_EVAL)
    BinaryRead2(file_head + ".ref.locs", read_POGs);

  // Sanity checks!
  // If any of these asserts fails, the input files are inconsistent.
  size_t n_reads = reads_bases.size();
  size_t n_unipaths = unipaths.size();
  ForceAssertEq(n_reads, paths.size());
  ForceAssertEq(n_reads, paths_rc.size());
  ForceAssertEq(n_reads, pairs.nReads());
  ForceAssertEq(n_reads, read_lengths.size());
  ForceAssertEq(n_unipaths, unibases.size());
  ForceAssertEq(n_unipaths, (size_t) LG.N());
  ForceAssertEq(n_unipaths, CNs.size());
  if (USE_TRUTH) ForceAssertEq(n_unipaths, unipath_POGs.size());

  SeedStore seeds(sub_dir, SEED_IDS,
                  DUMP_LOCAL_UNIBASES || NHOOD_EVAL ||
                  LOCAL_LOG || LOCAL_FASTA || LOCAL_PRIMARY || LOCAL_DOT);

  /*****************************************************************************
   *
   *        DATASET PRE-PROCESSING: Build some useful auxiliary data structures
   *
   ****************************************************************************/

  cout << Tag() << "Pre-processing dataset" << endl;

  // Find unipath lengths
  vec<int> unipath_lengths(n_unipaths);
  for (size_t i = 0; i < n_unipaths; i++)
    unipath_lengths[i] = unipaths[i].KmerCount();


  // Read the CNs and find the predicted copy number (CN) of each unipath
  vec<int> predicted_CNs(n_unipaths, -1);
  for (size_t i = 0; i < n_unipaths; i++)
    GetMostLikelyValue(predicted_CNs[i], CNs[i]);


  // Index of unilocs (indexed by unipath ID, not read ID)
  vec<longlong> unilocs_index(n_unipaths + 1, -1);
  unilocs_index[0] = 0;
  for (size_t i = 0; i < unilocs.size(); i++)
    unilocs_index[ unilocs[i].Contig() + 1 ] = i + 1;
  for (size_t i = 1; i <= n_unipaths; i++)
    if (unilocs_index[i] < 0) unilocs_index[i] = unilocs_index[i-1];

  // Find unipaths which are branches in the graph
  vec<Bool> branches(n_unipaths);
  digraph AG;
  BuildUnipathAdjacencyGraph(paths, paths_rc, pathsdb,
                              unipaths, unipathsdb, AG);
  for (size_t i = 0; i < n_unipaths; i++)
    branches[i] = (AG.From(i).size() > 1 || AG.To(i).size() > 1);

  // Load heuristic parameters into static SeedNeighborhood class variables
  SeedNeighborhood::_NHOOD_RADIUS           = NHOOD_RADIUS;
  SeedNeighborhood::_NHOOD_RADIUS_INTERNAL  = NHOOD_RADIUS_INTERNAL;
  SeedNeighborhood::_NHOOD_GLUEPERFECT_SIZE = NHOOD_GLUEPERFECT_SIZE;
  SeedNeighborhood::_NHOOD_VERBOSITY        = VERBOSITY;
  SeedNeighborhood::_DUMP_UNMERGED          = DUMP_UNMERGED;
  SeedNeighborhood::_DUMP_LOCAL_UNIBASES    = DUMP_LOCAL_UNIBASES;
  SeedNeighborhood::_USE_ACYCLIC            = USE_ACYCLIC;






  /*****************************************************************************
   *
   *        ASSEMBLE SEED NEIGHBORHOODS
   *
   ****************************************************************************/

  cout << endl
       << Tag() << "Setting up parallel runs for seeds."
       << endl << endl;
  double nhood_clock = WallClockTime( );

  // Place each localized assembly into the worklist.
  // As soon as it gets the first item, the worklist will immediately begin
  // processing the function Processor::operator(), which
  // performs each local assembly with the ID passed to it.  The worklist
  // runs its functions in multi-threaded parallel.
  LRLGCounters counters ( seeds.GetNumSeeds( ) );

  {
    const String log_fn = sub_dir + "/LocalizeReadsLG.log";
    cout << Tag() << "Sending log to: " << log_fn << endl;
    ofstream log_fs(log_fn.c_str());
    PrintCommandPretty(log_fs);
  
    Processor processor(seeds,
                        K, PLOIDY, MAX_COPY_NUMBER_OTHER,
                        & paths,
                        & paths_rc,
                        & pathsdb,
                        & unipath_lengths,
                        & read_lengths,
                        & pairs,
                        & unilocs,
                        & unilocs_index,
                        & LG,
                        & predicted_CNs,
                        & branches,
                        & reads_bases,
                        & unibases,
                        & unibases_next,
                        & to_rc,
                        EVAL, READ_CLOUD_EVAL,
                        LOCAL_DOT,
                        local_dump,
                        LOCAL_FASTA,
                        LOCAL_PRIMARY,
                        BRIDGE_BETWEEN_CLOUD_UNIPATHS,
                        LOCAL_LOG,
                        NHOOD_EVAL,
                        & genome,
                        & unipath_POGs,
                        & read_POGs,
                        & counters,
                        & log_fs);
    Worklist<size_t, Processor> worklist(processor, NUM_THREADS - 1, 100000 * 1024);
    worklist.threadDeathMonitor(60,true);

    for (size_t i = 0; i != NUM_THREADS; i++)
      worklist.add(i);
  }
  // The Worklist destructor doesn't return until it has done all its processing
  seeds.close();

  // Stop timer.
  timer.Stop();

  // Report runtime stats.
  cout << "\n"
       << "LOCALIZE READS LG REPORT\n"
       << "\n"
       << "Overall runtimes:\n"
       << "\n"
       << "  Real time:            " << timer.Elapsed( ) << " sec\n"
       << "  CPU time:             " << timer.UserSecs( ) << " sec\n"
       << "  Local assemblies:     " << TimeSince(nhood_clock) << "\n"
       << endl;

  // Nothing to report.
  if ( counters.n_inserts_walked == 0 ) {
    cout << endl << Tag() << "WARNING! No inserts were walked!" << endl;
  } else {

    cout << "Average runtimes of important local assembly steps (wallclock seconds):\n"
	 << "\n"
	 << "  Pathing local reads:  " << (counters.T_pathing / counters.n_seeds_done) << " sec\n"
	 << "  Walking inserts:      " << (counters.T_insert_walking / counters.n_seeds_done) << " sec\n"
	 << "  Merging inserts:      " << (counters.T_insert_merging / counters.n_seeds_done) << " sec\n"
	 << "\n";
    
    // Report success rate for insert walking.
    ReportWalkingRatesCore( cout, False, seeds.GetSeedIDs( ), counters.walks_info );
  }

  // Done.
  cout << Tag() << "LocalizeReadsLG is done" << endl;
  return 0;

}
